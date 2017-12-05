#!/usr/bin/env python

import numpy as np

# which simulation
whichsim = 0
# set up where the output files are
savedir = '/path/to/output/files'

# read in the galaxy catalog
fcent = np.fromfile(savedir+"/halos_gal_cent_"+str(whichsim))
fsats = np.fromfile(savedir+"/halos_gal_sats_"+str(whichsim))
# reshape the file data
fcent = np.array(np.reshape(fcent, (-1, 5)))
fsats = np.array(np.reshape(fsats, (-1, 5)))

# load the galaxy catalogs
pos_cent = fcent[:,0:3]           # central galaxy positions, Mpc
halo_indx_cent = fcent[:,3]       # central galaxy halo indices
halo_mass_cent = fcent[:,4]       # central galaxy halo mass, Msun

pos_sats = fsats[:,0:3]           # satellite galaxy positions, Mpc
halo_indx_sats = fsats[:,3]       # satellite galaxy halo indicies
halo_mass_sats = fsats[:,4]       # satellite galaxy halo mass, Msun
