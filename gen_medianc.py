#!/usr/bin/env python
"""
A module to compute a polynomial fit for the median halo concentration as 
a function of halo mass. 

Simply run 'python gen_medianc.py' to execute this script. The script will
output cparam_fits.npz file to your data directory that contains the best fit.

"""

import numpy as np
import os
import sys
import random
from astropy.table import Table
import astropy.io.fits as pf
import h5py
from glob import glob

def load_mc(whichsim, params):
    """
    Compiles the halo mass and concentration for a given simulation box.

    Parameters
    ----------

    whichsim : int
        Simulation number. Ranges between [0, 15] for current Planck 1100 sims.
    params : dict
        Dictionary of various simulation parameters. 

    Outputs
    -------

    mass_axis : array
        An array of halo mass bins. 

    cmedian_axis : array
        An array of median halo concentrations in the same mass bins. 

    """

    # directory of the halo and particle files
    directory = \
    "/mnt/store2/bigsim_products/emulator_1100box_planck_products/"\
    +"emulator_1100box_planck_00-"+str(whichsim)+"_products/"\
    +"emulator_1100box_planck_00-"+str(whichsim)+"_rockstar_halos/z0.500"

    # set up empty arrays
    allms = np.array([])
    allcs = np.array([])
    # loop over all the halos files and pull out the relevant data
    files = [h5py.File(fn) for fn in glob(directory+'/halos_0.*.h5')]   
    num_files = len(files)
    for i in range(0, num_files):
        # open the halo files
        newfile = h5py.File(directory+'/halos_0.'+str(i)+'.h5')
        halos = newfile['halos']

        # filter out the sub halos
        mask = halos['parent_id'] == -1
        maskedhalos = halos[mask]

        # extract concentration 
        halo_rs = maskedhalos["klypin_rs"]
        halo_rvir = maskedhalos["r"]
        halo_c = halo_rvir/halo_rs
        # extract halo mass
        halo_mass = maskedhalos['m']/params['h'] # msun
        # compile them in data array
        allms = np.concatenate((allms, halo_mass))
        allcs = np.concatenate((allcs, halo_c))

    # compile a list of mass and a list of median concentrations
    allms_log = np.log10(allms)
    # bin the mass in log bins 
    numbins = 100
    massbins = np.linspace(np.log10(4e12), np.max(allms_log), num=numbins+1)
    mass_axis = 0.5*(massbins[:-1]+massbins[1:]) # log10 (msun)
    cmedian_axis = np.zeros(numbins)
    # go through the list to compute the median concentration
    for i in range(0, numbins):
        newlo = massbins[i]
        newhi = massbins[i+1]

        # calculate the median concentration
        masks = [allms_log > newlo, allms_log < newhi]
        mask = reduce(np.logical_and, masks)
        mask = np.array(mask)

        newcs = allcs[mask]
        newc_med = np.median(newcs)

        # calculate the standard deviation in c
        newcstd = np.sqrt(np.mean(newcs**2) - np.mean(newcs)**2)

        cmedian_axis[i] = newc_med

    return mass_axis, cmedian_axis

def avg_c(params, rsd = True):
    """
    Function to calculate the polynomial fit to the median concentration 
    and save it to a file.

    Parameters
    ---------- 

    params : dict
        Dictionary of various simulation parameters. 

    rsd : boolean, optional
        Flag of whether to implement RSD. 

    """
    # compile the mass and median concentration from all sims.
    all_mass = 0
    all_cmedians = 0
    all_cmedians2 = 0
    for i in range(0, params['num_sims']):
        print "Loading data from simulation box ", i
        mass_axis, cmedian_axis = load_mc(i, params)
        all_mass = all_mass + mass_axis
        all_cmedians = all_cmedians + cmedian_axis
        all_cmedians2 = all_cmedians2 + cmedian_axis**2

    # calculate the average median c and mass. 
    c_med = all_cmedians/params['num_sims']
    mass = all_mass/params['num_sims']
    # calculate the uncertainty on median c
    var_cmedians = all_cmedians2/params['num_sims'] - c_med**2
    c_medsig = np.sqrt(var_cmedians)
    # remove nans
    masks = [~np.isnan(c_med), ~np.isnan(c_medsig)]
    totmask = reduce(np.logical_and, masks)
    mass1 = mass[totmask]
    c_med1 = c_med[totmask]
    c_medsig1 = c_medsig[totmask]
    # fit the mass-median c relationship with polynomial
    weights = 1.0/(1.25*c_medsig1[1:]/np.sqrt(params['num_sims']))
    p_c = np.polyfit(mass1[1:], c_med1[1:], 3, w = weights)
    # save the fit
    datadir = "./data"
    if rsd:
        datadir  = datadir + "_rsd"
    fitfile = datadir+"/cparam_fits"
    np.savez(fitfile, cfit = p_c)

# constants
params = {}
params['z'] = 0.5
params['h'] = 0.6726
params['Nslab'] = 3
params['Lbox'] = 1100/params['h'] # Mpc, box size
params['Mpart'] = 3.88537e+10/params['h'] # Msun, mass of each particle
params['velz2kms'] = 9.690310687246482e+04/params['Lbox'] # H(z)/(1+Z), km/s/Mpc
params['maxdist'] = 30 # Mpc
params['num_sims'] = 16

rsd = True 

avg_c(params, rsd)
