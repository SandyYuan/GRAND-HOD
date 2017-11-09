# all sims, testHOD
import numpy as np
import os,sys
import random
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib import rc, rcParams
rcParams.update({'font.size': 20})
from astropy.table import Table
import astropy.io.fits as pf
from astropy.cosmology import WMAP9 as cosmo
from scipy.special import erfc

import h5py
from glob import glob

# number of central galaxies, returns a fraction
def n_cen(M_in, design): 
	M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 'M1', 'sigma', 'alpha', 'kappa'))

	if M_in < 4e12: # this cutoff ignores halos with less than 100 particles
		return 0
	return .5*erfc(np.log(M_cut/M_in)/(2**.5*sigma))

# number of satellite galaxies, returns a float
def n_sat(M_in, design): 
	M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 'M1', 'sigma', 'alpha', 'kappa'))

	if M_in < kappa*M_cut:
		return 0
	else:
		return ((M_in - kappa*M_cut)/M1)**alpha

# an rsd function
def do_rsd(pos_vec, vel_vec, params):
	# pos gives arrays of x, y, z in Mpc
	# vel gives arrays of vx, vz, vz in km/s
	pos_vec[2] = (pos_vec[2] + vel_vec[2]/params['velz2kms']) % params['Lbox'] # wrap around
	return pos_vec

# generate central galaxy given a halo
def gen_cent(halo_ids, halo_pos, halo_vels, halo_vrms, halo_mass, design, decorations, fcent, rsd, params):
	for i in range(0, len(halo_pos)):
		# undecorated number of centrals
		N_cent = n_cen(halo_mass[i], design)
		# undecorated central pos and vel
		cent_pos = halo_pos[i]
		cent_vel = halo_vels[i]

		# generated peculiar velocity for velocity bias
		vrms_los = halo_vrms[i]/np.sqrt(3.0) # km/s
		# if there is no velocity bias, just throw a random number and move on
		if decorations['alpha_c'] == 0:
			np.random.normal(loc = 0, scale = 1) 
		else:
			v_pec = np.random.normal(loc = 0, scale = decorations['alpha_c']*vrms_los) # km / s
			cent_vel[2] = cent_vel[2] + v_pec 

		# if we do have a central, then store to file
		if np.random.random() < N_cent:
			# modify the pos with z distort
			if rsd:
				cent_pos = do_rsd(cent_pos, cent_vel, params)

			# write the central to file
			newline_cent = np.array([cent_pos[0], cent_pos[1], cent_pos[2], halo_ids[i], halo_mass[i]])
			newline_cent.tofile(fcent)


# generate satellites given an h5 file
def gen_sats(halo_ids, halo_pos, newfile, newpart, halo_mass, halo_pstart, halo_pnum, design, decorations, fsats, rsd, params):
	# standard hod design
	M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 'M1', 'sigma', 'alpha', 'kappa'))

	# process the subsample file to pull out the vels and pos
	subsample = newpart['subsamples']
	part_pos = subsample['pos']
	part_vel = subsample['vel']

	# loop through the halos to populate satellites
	for i in range(0, len(halo_ids)):
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
				# list of relative positions of the particles to the center of halo
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
				eachprob_array = eachprob_array*(1 + decorations['s']*(1 - temp_indices/(temp_range/2.0)))

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
				eachprob_array = eachprob_array*(1 + decorations['s_v']*(1 - temp_indices/(temp_range/2.0)))

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
				v_rel2 = vels_rel[:,0]**2 + vels_rel[:,1]**2 + vels_rel[:,2]**2 # speed

				# calculate radial and tangential peculiar velocity
				vel_rad = vels_rel*r_rel_norm
				v_rad2 = vel_rad[:,0]**2 + vel_rad[:,1]**2 + vel_rad[:,2]**2 # speed
				v_tan2 = v_rel2 - v_rad2

				# compute the perihelion distance for NFW profile
				m = halo_mass[i]
				c = halo_c[i]
				rs = halo_rs[i]
				r0_kpc = r0*1000 # kpc
				alpha = 1.0/(np.log(1+c)-c/(1+c))*2*6.67e-11*m*2e30/r0_kpc/3.086e+19/1e6 # (km/s)^2

				# iterate a few times to solve for rp
				x2 = v_tan2/(v_tan2+v_rad2)
				num_iters = 10 # how many iterations do we want
				for i in range(0, num_iters):
					oldx = np.sqrt(x2)
					x2 = v_tan2/(v_tan2 + v_rad2 + alpha*(np.log(1+oldx*r0_kpc/rs)/oldx - np.log(1+r0_kpc/rs)))

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
				eachprob_array = eachprob_array*(1 + decorations['s_p']*(1 - temp_indices/(temp_range/2.0)))

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
			sat_pos[:,2] = (sat_pos[:,2] + sat_vels[:,2]/params['velz2kms']) % params['Lbox']

		# then output to file
		for j in range(0, len(sat_pos)):
			newline_sat = np.array([sat_pos[j, 0],sat_pos[j, 1],sat_pos[j, 2], halo_ids[i], halo_mass[i]])
			newline_sat.tofile(fsats)


# Compile the full halo catalog, and populate each halo with galaxies
def gen_gals(directory, design, decorations, fcent, fsats, rsd, params):
	M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 'M1', 'sigma', 'alpha', 'kappa'))
	s, s_v, alpha_c, s_p, A = map(decorations.get, ('s', 's_v', 'alpha_c', 's_p', 'A'))

	# loop over all the halos files and pull out the relevant data fields
	# now instead of processing each file separately, we are concatenating them into one big array
	# this is easier for ranking and re-ranking

	files = [h5py.File(fn) for fn in glob(directory+'/halos_0.*.h5')]	
	num_files = len(files)
	for i in range(0, num_files):
		# open the halo files
		newfile = h5py.File(directory+'/halos_0.'+str(i)+'.h5')
		halos = newfile['halos']

		# this halos catalog still contain the subhalos, we want to filter that out
		masks = [halos['parent_id'] == -1, halos['m']/params['h'] > 4e12] 
		mask = reduce(np.logical_and, masks)
		maskedhalos = halos[mask]
		halo_ids = maskedhalos["id"]
		halo_pos = maskedhalos["pos"]/params['h'] # Mpc
		halo_vels = maskedhalos["vel"] # km/s
		halo_vrms = maskedhalos["vrms"] # km/s
		halo_mass = maskedhalos['m']/params['h'] # msun
		halo_r = maskedhalos['r']/params['h'] # Mpc
		halo_rs = maskedhalos['klypin_rs']/params['h'] # Mpc
		halo_c = halo_r/halo_rs
		halo_pstart = maskedhalos['subsamp_start']
		halo_pnum = maskedhalos['subsamp_len']

		# if assembly bias parameter is not zero, then we do halo reranking
		if not A == 0:
			# now that we have halo mass array, we could just calculate c median and std.
			halo_cmed = np.polyval(pcmed, np.log10(halo_mass))
			# then we can calculate the quantity halo_deltac
			halo_deltac = halo_c - halo_cmed

			# define a ranking parameter
			halo_pseudomass = halo_mass*np.exp(A*(2*(halo_deltac > 0) - 1))

			# create a list that indicates the original order 
			halo_order = np.array(range(0, len(halo_ids)))

			# first we sort everything by mass, original mass
			msorted_indices = halo_mass.argsort()[::-1] # descending order
			# now build the halo catalog with the sorted indices
			halo_ids = halo_ids[msorted_indices]
			halo_pos = halo_pos[msorted_indices]
			halo_vels = halo_vels[msorted_indices]
			halo_mass = halo_mass[msorted_indices] # we have sorted halo_mass by mass
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

			# revert to the original order so we dont mess up the reseeding procedure
			original_indices = halo_order.argsort() # ascending order
			halo_ids = halo_ids[original_indices]
			halo_pos = halo_pos[original_indices]
			halo_vels = halo_vels[original_indices]
			halo_pstart = halo_pstart[original_indices]
			halo_pnum = halo_pnum[original_indices]
			# we sort the halo mass too 
			halo_mass = halo_mass[original_indices]

		# for each halo, generate central galaxies and output to file
		gen_cent(halo_ids, halo_pos, halo_vels, halo_vrms, halo_mass, design, decorations, fcent, rsd, params)
		# open particle file
		newpart = h5py.File(directory+'/particles_0.'+str(i)+'.h5')
		# for each halo, generate satellites and output to file
		gen_sats(halo_ids, halo_pos, newfile, newpart, halo_mass, halo_pstart, halo_pnum, design, decorations, fsats, rsd, params)


######## generate a central galaxy catalog and a satellite galaxy catalog
# (1) for the central galaxy catalog, we will have, galaxy position, number of galaxies, 
#     halo index, halo mass
# (2) for the satellites catalog, we will have, galaxy position, halo index, halo postion, halo mass
def gen_gal_cat(whichsim, design, decorations, params, whatseed = 0, rsd = True):
	M_cut, M1, sigma, alpha, kappa = map(design.get, ('M_cut', 'M1', 'sigma', 'alpha', 'kappa'))
	s, s_v, alpha_c, s_p, A = map(decorations.get, ('s', 's_v', 'alpha_c', 's_p', 'A'))

	# see the RNG
	np.random.seed(seed = whatseed)

	print "Generating galaxy catalog. Realization: ", whichsim
	directory = "/mnt/store2/bigsim_products/emulator_1100box_planck_products/emulator_1100box_planck_00-"\
	+str(whichsim)+"_products/emulator_1100box_planck_00-"+str(whichsim)+"_rockstar_halos/z0.500"

	# create a save directory
	datadir = "./data"
	if rsd:
		datadir = datadir+"_rsd"
	savedir = datadir+"/rockstar_"\
	+str(M_cut)[0:4]+"_"+str(M1)[0:4]+"_"+str(sigma)[0:4]+"_"+str(alpha)[0:4]+"_"+str(kappa)[0:4]\
	+"_decor_"+str(s)+"_"+str(s_v)+"_"+str(alpha_c)+"_"+str(s_p)+"_"+str(A)
	if rsd:
		savedir = savedir+"_rsd"
	# if we are doing repeats, save them in separate directories
	if not whatseed == 0:
		savedir = savedir+"_"+str(whatseed)

	# if this directory does not exist, make it
	if not os.path.exists(savedir):
		os.makedirs(savedir)

	# binary output galaxy catalog
	print "Building galaxy catalog (binary output)"
	fcent = open(savedir+"/halos_gal_cent_"+str(whichsim),'wb')
	fsats = open(savedir+"/halos_gal_sats_"+str(whichsim),'wb')

	# find the halos, populate them with galaxies and write them to files
	gen_gals(directory, design, decorations, fcent, fsats, rsd, params)

	fcent.close()
	fsats.close()
	print "Galaxy Catalogs Done. Realization: ", whichsim



"""
# 3d scatter plot of all particles within this halo and the halo pos
fig = pl.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.scatter(pos[:,0],pos[:,1],pos[:,2], s=10)
#ax.scatter(halo_pos_all[0][0],halo_pos_all[0][1],halo_pos_all[0][2], s=20, c='r')
pl.show()
"""
