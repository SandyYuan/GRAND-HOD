#!/usr/bin/env python

from GRAND_HOD import gen_gal_catalog_rockstar as galcat
from GRAND_HOD.gen_medianc import avg_c

import numpy as np

# do you want to invoke rsd?
rsd = True
# constants
h = 0.6726
params = { 'z': 0.5,
           'h': h,
           'Nslab': 3,                            # number of data files per simulation box
           'Lbox': 1100/h,                        # Mpc, box size
           'Mpart': 3.88537e+10/h,                # Msun, mass of each particle
           'num_sims': 16,
           'rsd': rsd}

# parameter for converting velocity in km/s to position in Mpc
params['velz2kms'] = 9.690310687246482e+04/params['Lbox']   # H(z)/(1+Z), km/s/Mpc

# baseline HOD  (Zheng+2009, Kwan+2015)
M_cut = 10**13.35 
log_Mcut = np.log10(M_cut)
M1 = 10**13.8
log_M1 = np.log10(M1)
sigma = 0.85
alpha = 1.0
kappa = 1.0

# generalized HOD prescription 
design = {'M_cut': M_cut, 'M1': M1, 'sigma': sigma, 'alpha': alpha, 'kappa': kappa}
decorations = {'s': 0., 's_v': 0., 'alpha_c': 0., 's_p': 0., 'A': 0.}

# which simulation box are we using?
whichsim = 0

# compute the median halo concentration fit, you dont need to run this if you dont 
# plan on invoking assembly bias decoration. 
#avg_c(params, rsd)

# generate galaxy catalogs
galcat.gen_gal_cat(whichsim, design, decorations, params, rsd = rsd)
