#!/usr/bin/env python

"""
Module implementation of a generalized and differentiable Halo Occupation 
Distribution (HOD)for N-body cosmological simulations. 

Add to .bashrc:
export PYTHONPATH="/path/to/GRAND-HOD:$PYTHONPATH"

"""

import numpy as np
import os
import sys
import random
import time
from astropy.table import Table
import astropy.io.fits as pf
from scipy.special import erfc
import h5py
from glob import glob


def n_cen(M_in, design, m_cutoff=4e12): 
    """
    Computes the expected number of central galaxies given a halo mass and 
    the HOD design. 

    Parameters
    ----------

    M_in : float
        Halo mass in solar mass.

    design : dict
        Dictionary containing the five HOD parameters. 
        
    m_cutoff: float, optional
        Ignore halos small than this mass.

    Returns
    -------

    n_cen : float
        Number of centrals expected for the halo within the range (0, 1).
        This number should be interpreted as a probability.

    """
    M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 
                                                      'M1', 
                                                      'sigma', 
                                                      'alpha', 
                                                      'kappa'))
    #M_cut, M1, sigma, alpha, kappa = [design[k] for k in ('M_cut', 
    #                                                      'M1', 
    #                                                      'sigma', 
    #                                                      'alpha', 
    #                                                      'kappa' ) ]

    if M_in < m_cutoff: # this cutoff ignores halos with less than 100 particles
        return 0
    return .5*erfc(np.log(M_cut/M_in)/(2**.5*sigma))

def n_sat(M_in, design): 
    """
    Computes the expected number of satellite galaxies given a halo mass and 
    the HOD design. 

    Parameters
    ----------

    M_in : float
        Halo mass in solar mass.

    design : dict
        Dictionary containing the five HOD parameters. 
    
    Returns
    -------

    n_sat : float
        Expected number of satellite galaxies for the said halo.

    """
    M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 
                                                      'M1', 
                                                      'sigma', 
                                                      'alpha', 
                                                      'kappa'))

    if M_in < kappa*M_cut: # mass cut off
        return 0
    else:
        return ((M_in - kappa*M_cut)/M1)**alpha

def do_rsd(pos_vec, vel_vec, params):
    """
    Perform Redshift Space Distortion (RSD) on a given position and a velocity.
    The z direction position is modified by the z direction velocity using the
    following formula:

        z' = z + v_z/H_0

    z' also wraps around to the other side of the simulation box.

    Parameters
    ----------

    pos_vec : numpy.array
        Particle position in array of length 3, (x, y, z). Each component
        is given in box unit, i.e. between (-0.5, 0.5).

    vel_vec : numpy.array
        Particle velocity in array of length 3, (vx, vy, vz). Each component 
        is given in km/s.

    params : dict
        Dictionary of simulation parameters. 

    Returns
    -------

    pos_vec : numpy.array
        The updated position vector of the particle.

    """
    pos_vec[2] = \
    (pos_vec[2] + vel_vec[2]/params['velz2kms']) % params['Lbox'] 
    return pos_vec

# generate central galaxy given a halo
def gen_cent(halo_ids, halo_pos, halo_vels, halo_vrms, halo_mass, 
             design, decorations, fcent, rsd, params):
    """
    Function that generates central galaxies and its position and velocity 
    given a halo catalog and HOD designs and decorations. The generated 
    galaxies are output to file fcent. 

    Parameters
    ----------

    halo_ids : numpy.array
        Array of halo IDs.

    halo_pos : numpy.array
        Array of halo positions of shape (N, 3) in box units.

    halo_vels : numpy.array
        Array of halo velocities of shape (N, 3) in km/s.

    halo_vrms: numpy.array
        Array of halo particle velocity dispersion in km/s.

    halo_mass : numpy.array
        Array of halo mass in solar mass.

    design : dict
        Dictionary of the five baseline HOD parameters. 

    decorations : dict
        Dictionary of generalized HOD parameters. 

    fcent : file pointer
        Pointer to the central galaxies output file location. 

    rsd : boolean
        Flag of whether to implement RSD. 

    params : dict
        Dictionary of various simulation parameters. 

    Outputs
    -------

    For each halo, if there exists a central, the function outputs the 
    3D position (Mpc), halo ID, and halo mass (Msun) to file.

    """

    # loop over all halos in the catalog
    for i in range(0, len(halo_pos)):
        # undecorated number of centrals
        N_cent = n_cen(halo_mass[i], design)
        # undecorated central pos and vel
        cent_pos = halo_pos[i]
        cent_vel = halo_vels[i]
        # halo velocity dispersion along the LOS
        vrms_los = halo_vrms[i]/np.sqrt(3.0) # km/s

        # if there is no velocity bias, just throw a random number and move on
        if decorations['alpha_c'] == 0:
            np.random.normal(loc = 0, scale = 1) 
        # if there is velocity, generate pecular velocity for the central
        else:
            v_pec = np.random.normal(loc = 0, 
                            scale = abs(decorations['alpha_c'])*vrms_los)
            cent_vel[2] = cent_vel[2] + v_pec 

        # if we do have a central, then store to file
        if np.random.random() < N_cent:
            # modify the pos with z distort
            if rsd:
                cent_pos = do_rsd(cent_pos, cent_vel, params)

            # write the central to file
            newline_cent = np.array([cent_pos[0], 
                                     cent_pos[1], 
                                     cent_pos[2],  
                                     halo_ids[i], 
                                     halo_mass[i]])

            newline_cent.tofile(fcent)


def gen_sats(halo_ids, halo_pos, halo_vels, newpart, halo_mass, 
             halo_pstart, halo_pnum, design, decorations, fsats, rsd, params):

    """
    Function that generates satellite galaxies and their positions and 
    velocities given a halo catalog and HOD designs and decorations. 

    The decorations are implemented using a particle re-ranking procedure
    that preserves the random number thrown for each particle so that 
    the resulting statistic has no induced shot-noise and is thus 
    differentiable.

    The generated galaxies are output to binary file fsats. 

    Parameters
    ----------

    halo_ids : numpy.array
        Array of halo IDs.

    halo_pos : numpy.array
        Array of halo positions of shape (N, 3) in box units.

    halo_vels : numpy.array
        Array of halo velocities of shape (N, 3) in km/s.

    newpart : h5py file pointer
        The pointer to the particle file.

    halo_mass : numpy.array
        Array of halo mass in solar mass.

    halo_pstart : numpy.array
        Array of particle start indices for each halo.

    halo_pnum : numpy.array
        Array of number of particles for halos. 

    design : dict
        Dictionary of the five baseline HOD parameters. 

    decorations : dict
        Dictionary of generalized HOD parameters. 

    fsats : file pointer
        Pointer to the satellite galaxies output file location. 

    rsd : boolean
        Flag of whether to implement RSD. 

    params : dict
        Dictionary of various simulation parameters. 

    Outputs
    -------

    For each halo, the function outputs the satellite galaxies, specifically
    the 3D position (Mpc), halo ID, and halo mass (Msun) to file.


    """

    # standard hod design
    M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 
                                                      'M1', 
                                                      'sigma', 
                                                      'alpha', 
                                                      'kappa'))

    # process the subsample file to pull out the vels and pos
    subsample = newpart['subsamples']
    part_pos = subsample['pos']
    part_vel = subsample['vel']

    # loop through the halos to populate satellites
    for i in range(len(halo_ids)):
        # load the 10% subsample belonging to the halo
        start_ind = halo_pstart[i]
        numparts = halo_pnum[i]
        # if there are no particles in the particle subsample, move on
        if numparts == 0:
            continue
        # extract the particle positions and vels
        ss_pos = part_pos[start_ind: start_ind + numparts]/params['h'] # Mpc
        ss_vels = part_vel[start_ind: start_ind + numparts] # km/s

        # compute the undecorated expected number of satellites
        N_sat = n_sat(halo_mass[i], design)

        # generate a list of random numbers that will track each particle
        random_list = np.random.random(numparts)

        # the undecorated probability of each particle hosting a satellite
        eachprob = float(N_sat)/numparts
        eachprob_array = np.ones(numparts)*eachprob

        temp_indices = np.arange(numparts)
        temp_range = numparts - 1

        # if there is one particle, then we dont need to do any reranking
        if numparts > 1:
            # first ranking parameter, distance to halo center
            if not decorations['s'] == 0:
                # list of relative positions of the particles to center of halo
                pos_rel = ss_pos - halo_pos[i]
                # distance to the center
                dist2_rel = pos_rel[:,0]**2 + pos_rel[:,1]**2 + pos_rel[:,2]**2

                # now we rank the list by relative distance
                sorted_indices = dist2_rel.argsort()[::-1] # furthest to closest
                # rerank the list of random numbers 
                random_list = random_list[sorted_indices]
                eachprob_array = eachprob_array[sorted_indices]
                ss_pos = ss_pos[sorted_indices]
                ss_vels = ss_vels[sorted_indices]
    
                # tilt probability distribution
                eachprob_array = eachprob_array*\
                (1 + decorations['s']*(1 - temp_indices/(temp_range/2.0)))

            # second ranking parameter s_v
            if not decorations['s_v'] == 0:
                # list of peculiar velocities of the particles
                vels_rel = ss_vels - halo_vels[i]
                # speed relative to center
                v2_rel = vels_rel[:,0]**2 + vels_rel[:,1]**2 + vels_rel[:,2]**2

                # now we rank the list by relative velocity
                sorted_indices = v2_rel.argsort()[::-1] # highest to closest
                # rerank the list of random numbers 
                random_list = random_list[sorted_indices]
                eachprob_array = eachprob_array[sorted_indices]
                ss_pos = ss_pos[sorted_indices]
                ss_vels = ss_vels[sorted_indices]

                # tilt probability distribution
                eachprob_array = eachprob_array*\
                (1 + decorations['s_v']*(1 - temp_indices/(temp_range/2.0)))

            # third ranking parameter s_p
            if not decorations['s_p'] == 0:
                # calc relative positions
                r_rel = ss_pos - halo_pos[i] # Mpc
                r0 = np.sqrt(r_rel[:,0]**2+r_rel[:,1]**2+r_rel[:,2]**2) # Mpc
                r_rel_norm = np.zeros(np.shape(r_rel))
                r_rel_norm[:,0] = r_rel[:,0]/r0 # normalized relative positions
                r_rel_norm[:,1] = r_rel[:,1]/r0
                r_rel_norm[:,2] = r_rel[:,2]/r0

                # list of peculiar velocities of the particles
                vels_rel = ss_vels - halo_vels[i] # velocity km/s
                # relative speed to halo center squared
                v_rel2 = vels_rel[:,0]**2 + vels_rel[:,1]**2 + vels_rel[:,2]**2

                # calculate radial and tangential peculiar velocity
                vel_rad = vels_rel*r_rel_norm
                # radial component
                v_rad2 = vel_rad[:,0]**2 + vel_rad[:,1]**2 + vel_rad[:,2]**2
                # tangential component
                v_tan2 = v_rel2 - v_rad2

                # compute the perihelion distance for NFW profile
                m = halo_mass[i]
                c = halo_c[i]
                rs = halo_rs[i]
                r0_kpc = r0*1000 # kpc
                alpha = \
                1.0/(np.log(1+c)-c/(1+c))*2*6.67e-11*m*2e30/r0_kpc/3.086e+19/1e6

                # iterate a few times to solve for rp
                x2 = v_tan2/(v_tan2+v_rad2)
                num_iters = 10 # how many iterations do we want
                for i in range(0, num_iters):
                    oldx = np.sqrt(x2)
                    x2 = v_tan2/(v_tan2 + v_rad2 + alpha*\
                        (np.log(1+oldx*r0_kpc/rs)/oldx - np.log(1+r0_kpc/rs)))

                # final perihelion distance 
                rp2 = r0_kpc**2*x2

                # now we rank the list by perihelion distance to center
                sorted_indices = rp2.argsort()[::-1] # highest to closest
                # rerank the list of random numbers 
                random_list = random_list[sorted_indices]
                eachprob_array = eachprob_array[sorted_indices]
                ss_pos = ss_pos[sorted_indices]
                ss_vels = ss_vels[sorted_indices]

                # tilt probability distribution
                eachprob_array = eachprob_array*\
                (1 + decorations['s_p']*(1 - temp_indices/(temp_range/2.0)))

        # we have finished the reranking routines
        # decide which particle bears galaxy
        newmask = random_list < eachprob_array
        # generate the position and velocity of the galaxies
        sat_pos = ss_pos[newmask]
        sat_vels = ss_vels[newmask]

        # so a lot of the sat_pos are empty, in that case, just pass
        if len(sat_pos) == 0:
            continue

        # rsd, modify the satellite positions by their velocities
        if rsd:
            sat_pos[:,2] = \
            (sat_pos[:,2] + sat_vels[:,2]/params['velz2kms']) % params['Lbox']

        # then output to file
        for j in range(0, len(sat_pos)):
            newline_sat = np.array([sat_pos[j, 0],
                                    sat_pos[j, 1],
                                    sat_pos[j, 2], 
                                    halo_ids[i], 
                                    halo_mass[i]])

            newline_sat.tofile(fsats)


def gen_gals(directory, design, decorations, fcent, fsats, rsd, params):
    """
    Function that compiles halo catalog from directory and implements 
    assembly bias decoration by halo re-ranking. 

    The halo catalogs are then passed on to functions that generate central
    and satellite galaxies. 

    Parameters
    ----------

    directory : string 
        Directory of the halo and particle files. 

    design : dict
        Dictionary of the five baseline HOD parameters. 

    decorations : dict
        Dictionary of generalized HOD parameters. 

    fcent : file pointer
        Pointer to the central galaxies output file location. 

    fsats : file pointer
        Pointer to the satellite galaxies output file location. 

    rsd : boolean
        Flag of whether to implement RSD. 

    params : dict
        Dictionary of various simulation parameters. 


    """

    M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 
                                                      'M1', 
                                                      'sigma', 
                                                      'alpha', 
                                                      'kappa'))

    s, s_v, alpha_c, s_p, A = map(decorations.get, ('s', 
                                                    's_v', 
                                                    'alpha_c', 
                                                    's_p', 
                                                    'A'))

    # loop over all the halos files and pull out the relevant data fields
    files = [h5py.File(fn) for fn in glob(directory+'/halos_0.*.h5')]   
    num_files = len(files)
    for i in range(0, num_files):
        # open the halo files
        newfile = h5py.File(directory+'/halos_0.'+str(i)+'.h5')
        halos = newfile['halos']

        # removing subhalos from the catalog
        masks = [halos['parent_id'] == -1, halos['m']/params['h'] > 4e12] 
        mask = reduce(np.logical_and, masks)
        maskedhalos = halos[mask]
        # extracting the halo properties that we need
        halo_ids = maskedhalos["id"] # halo IDs
        halo_pos = maskedhalos["pos"]/params['h'] # halo positions, Mpc
        halo_vels = maskedhalos["vel"] # halo velocities, km/s
        halo_vrms = maskedhalos["vrms"] # halo velocity dispersions, km/s
        halo_mass = maskedhalos['m']/params['h'] # halo mass, Msun
        halo_r = maskedhalos['r']/params['h'] # virial radius, kpc
        halo_rs = maskedhalos['klypin_rs']/params['h'] # scale radius, kpc
        halo_c = halo_r/halo_rs # halo concentration
        halo_pstart = maskedhalos['subsamp_start'] # starting index of particles
        halo_pnum = maskedhalos['subsamp_len'] # number of particles 

        # if assembly bias parameter is not zero, then we do halo reranking
        if not A == 0:
            # calculate halo concentration median
            datadir = "./data"
            if rsd:
                datadir = datadir+"_rsd"
            pcfilename = datadir+"/cparam_fits.npz"
            # if the file exists, just load it
            if os.path.isfile(pcfilename):
                cfits = np.load(pcfilename)
                pcmed = cfits['cfit']
            # if the file does not exist, generate it
            else: 
                sys.exit("Error: to invoke assembly bias, \
                            you need to run gen_medianc.py first!")

            halo_cmed = np.polyval(pcmed, np.log10(halo_mass))
            # then we can calculate the quantity halo_deltac
            halo_deltac = halo_c - halo_cmed

            # define a ranking parameter
            halo_pseudomass = halo_mass*np.exp(A*(2*(halo_deltac > 0) - 1))

            # create a list that indicates the original order 
            halo_order = np.arange(len(halo_ids))

            # first we sort everything by mass, original mass
            # TODO: do we need to sort everything every time?
            msorted_indices = halo_mass.argsort()[::-1] # descending order
            # now build the halo catalog with the sorted indices
            halo_ids = halo_ids[msorted_indices]
            halo_pos = halo_pos[msorted_indices]
            halo_vels = halo_vels[msorted_indices]
            halo_mass = halo_mass[msorted_indices] 
            halo_pstart = halo_pstart[msorted_indices]
            halo_pnum = halo_pnum[msorted_indices]
            halo_pseudomass = halo_pseudomass[msorted_indices]
            halo_order = halo_order[msorted_indices]

            # now we resort using halo_pseudomass and get the indices
            new_indices = halo_pseudomass.argsort()[::-1] # descending order
            # we use these indices to reorder every array except the mass array
            halo_ids = halo_ids[new_indices]
            halo_pos = halo_pos[new_indices]
            halo_vels = halo_vels[new_indices]
            halo_pstart = halo_pstart[new_indices]
            halo_pnum = halo_pnum[new_indices]
            halo_order = halo_order[new_indices]
            # we dont touch halo mass so it is still sorted by mass

            # revert to the original order 
            original_indices = halo_order.argsort() # ascending order
            halo_ids = halo_ids[original_indices]
            halo_pos = halo_pos[original_indices]
            halo_vels = halo_vels[original_indices]
            halo_pstart = halo_pstart[original_indices]
            halo_pnum = halo_pnum[original_indices]
            # we sort the halo mass too 
            halo_mass = halo_mass[original_indices]

        # for each halo, generate central galaxies and output to file
        gen_cent(halo_ids, halo_pos, halo_vels, halo_vrms, halo_mass, 
                 design, decorations, fcent, rsd, params)

        # open particle file
        newpart = h5py.File(directory+'/particles_0.'+str(i)+'.h5', 'r')

        # for each halo, generate satellites and output to file
        gen_sats(halo_ids, halo_pos, halo_vels, newpart, halo_mass, 
                 halo_pstart, halo_pnum, design, decorations, fsats, 
                 rsd, params)
        newpart.close()


def gen_gal_cat(whichsim, design, decorations, params, 
                whatseed = 0, rsd = True,
                product_dir="/mnt/store2/bigsim_products/emulator_1100box_planck_products/",
                simname = "emulator_1100box_planck_00"):
    """
    Main interface that takes in the simulation number, HOD design and 
    decorations and outputs the resulting central and satellite galaxy
    catalogs to file in binary format. 

    Parameters
    ----------

    whichsim : int
        Simulation number. Ranges between [0, 15] for current Planck 1100 sims.

    design : dict
        Dictionary of the five baseline HOD parameters. 

    decorations : dict
        Dictionary of generalized HOD parameters. 

    params : dict
        Dictionary of various simulation parameters. 

    whatseed : integer, optional
        The initial seed to the random number generator. 

    rsd : boolean, optional
        Flag of whether to implement RSD. 

    product_dir : string, optional
        A string indicating the location of the simulation data. 
        You should not need to change this if you are on Eisenstein group clusters.

    simname : string, optional
        The name of the simulation boxes. Defaulted to 1100 planck boxes.

    """

    # checking for errors
    if not type(whichsim) is int or whichsim < 0:
        print "Error: whichsim has to be a non-negative integer."

    if not type(rsd) is bool:
        print "Error: rsd has to be a boolean."


    M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 
                                                      'M1', 
                                                      'sigma', 
                                                      'alpha', 
                                                      'kappa'))

    s, s_v, alpha_c, s_p, A = map(decorations.get, ('s', 
                                                    's_v', 
                                                    'alpha_c', 
                                                    's_p', 
                                                    'A'))

    # seed the RNG
    np.random.seed(seed = whatseed)

    whichsim = int(whichsim)
    print "Generating galaxy catalog. Realization: ", whichsim

    # directory of the halo and particle files
    directory = product_dir \
    +simname + "-"+str(whichsim)+"_products/"\
    +simname + "-"+str(whichsim)\
    +"_rockstar_halos/z{:.3f}".format(params['z'])

    # create a output directory
    datadir = "./data"
    if rsd:
        datadir = datadir+"_rsd"
    savedir = datadir+"/rockstar_"+str(M_cut)[0:4]+"_"+str(M1)[0:4]+"_"\
    +str(sigma)[0:4]+"_"+str(alpha)[0:4]+"_"+str(kappa)[0:4]+"_decor_"\
    +str(s)+"_"+str(s_v)+"_"+str(alpha_c)+"_"+str(s_p)+"_"+str(A)
    if rsd:
        savedir = savedir+"_rsd"
    # if we are doing repeats, save them in separate directories
    if not whatseed == 0:
        savedir = savedir+"_"+str(whatseed)

    # if this directory does not exist, make it
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    #try:
    #    os.makedirs(savedir)
    #except OSError:
    #    pass

    # binary output galaxy catalog
    print "Building galaxy catalog (binary output)"
    fcent = open(savedir+"/halos_gal_cent_"+str(whichsim),'wb')
    fsats = open(savedir+"/halos_gal_sats_"+str(whichsim),'wb')

    # find the halos, populate them with galaxies and write them to files
    gen_gals(directory, design, decorations, fcent, fsats, rsd, params)

    # close the files in the end
    fcent.close()
    fsats.close()
    print "Galaxy Catalogs Done. Realization: ", whichsim



