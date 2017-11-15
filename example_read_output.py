#!/usr/bin/env python

import numpy as np

whichsim = 0
savedir = '/path/to/output/files'

# read in the galaxy catalog
fcent = np.fromfile(savedir+"/halos_gal_cent_"+str(whichsim))
fsats = np.fromfile(savedir+"/halos_gal_sats_"+str(whichsim))
# reshape the file data
fcent = np.array(np.reshape(fcent, (-1, 5)))
fsats = np.array(np.reshape(fsats, (-1, 5)))

pos_cent = fcent[:,0:3]
halo_indx_cent = fcent[:,3]
halo_mass_cent = fcent[:,4]

pos_sats = fsats[:,0:3]
halo_indx_sats = fsats[:,3]
halo_mass_sats = fsats[:,4]
