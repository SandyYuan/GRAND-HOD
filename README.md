# GRAND-HOD: GeneRalized ANd Differentiable Halo Occupation Distribution

## Authors:
Sihan (Sandy) Yuan, Daniel Eisenstein & Lehman Garrison

## Introduction:
An Halo Occupation Distribution (HOD) code in Python that is differentiable and incorporates various generalizations to the standard HOD. The code takes a generalized HOD prescription as input and outputs the corresponding mock galaxy catalogs in binary files.

The code is currently written specifically for the Abacus simulations, but the main functionalities can be also easily adapted for other halo catalogs with the appropriate properties. We will update the code soon to support halo catalogs from other cosmological simulations. 

## Usage:

The code does not currently have dependencies other than basic Python packages. 

To install, simply download the .py files to the directory you want the mock catalogs to live in. If you are not on the Eisenstein Group computer clusters at CfA, you may need to change the `directory` variable to point to the location of the simulation data. Keep the simulation box tag `str(whichsim)` as an unknown variable. 

### Input:
The main interface of the code is the function `gen_gal_cat()`, which takes the following inputs:
- `whichsim` : integer. The index of the simulation box. For the Abacus 1100/h Mpc simulation with Planck cosmology, this number ranges from 0 to 15.  
- `design` : dictionary. The baseline HOD parameters. The dictionary requires the following five parameters:
  - `M_cut` : The cut-off mass for the halo to host in a central galaxy. Given in solar mass.
  - `M1` : The scale mass for the number of satellite galaxies. Given in solar mass.
  - `sigma` : Parameter that modulates the shape of the number of central galaxy.
  - `alpha` : The power law index of the number of satellite galaxies. 
  - `kappa` : Parameter that affects the cut-off mass for satellite galaxies. 
  A detailed discussion and best-fit values of these parameters can be found in Zheng+2009.
- `decorations` : dictionary. The HOD generalization parameters. The dictionary contains the following five parameters:
  - `s` : float. The satellite profile modulation parameter. Modulates how the radial distribution of satellite galaxies within halos deviate from the radial profile of the halo. Positive value favors satellite galaxies to populate the outskirts of the halo whereas negative value favors satellite galaxy to live near the center of the halo. 
  - `s_v` : float. The satellite velocity bias parameter. Modulates how the satellite galaxy peculiar velocity deviates from that of the local dark matter particle. Positive value favors high peculiar velocity satellite galaxies and vice versa. Note that our implementation preserves the Newton's second law of the satellite galaxies.
  - `alpha_c` : float. The central velocity bias parameter. Modulates the peculiar velocity of the central galaxy. The larger the absolute value of the parameter, the larger the peculiar velocity of the central. The sign of the value does not matter. 
  - `s_p` : float. The perihelion distance modulation parameter. A positive value favors satellite galaxies to have larger distances to the halo center upon their closest approach to the center and vice versa. This can be regarded as a "fancier" satellite profile modulation parameter. 
  - `A` : float. The assembly bias parameter. Introduces the effect of assembly bias. A positive value favors higher concentration halos to host galaxies whereas a negative value favors lower concentration halos to host galaxies. If you are invoking assembly bias decoration, i.e. a non-zero `A` parameter, you need to run gen_medianc.py first. 
  A detailed discussion of these parameters can be found in Yuan et al. in prep. To turn off any of the five decorations, just set the corresponding parameter to 0.
- `params` : dictionary. Simulation parameters. Following are the required parameters:
  - `z` : float. Redshift of the simulation snapshot. With the current directory, `z = 0.5`.
  - `h` : float. Hubble constant. For Planck cosmology, we use `h = 0.6726`.
  - `Nslab` : float. The number of halo data files each simulation box produces. For current Abacus 1100 boxes, `Nslab = 3`.
  - `Lbox` : float. The size of the simulation box in Mpc. For current Abacus 1100 boxes, `Lbox = 1100/h`.
  - `Mpart` : float. The particle mass in solar mass. For current Abacus 1100 boxes, `Mpart = 3.88e10/h`.
  - `velz2kms` : float. The conversion from simulation velocity to velocity in km/s and vice versa. The value can be calculated by H(z)/(1+z). 
  - `num_sims` : integer. The number of simulation boxes. For current Abacus 1100 boxes, `num_sims = 16`.
- `whatseed` : integer (optional). The seed to the random number generator. Default value is 0.
- `rsd` : boolean (optional). The redshift space distortion flag. Shifts the LOS locations of galaxies. Default is True. 
  
An example code to read these data files is
```
import gen_gal_catlog_rockstar as galcat
from gen_medianc import avg_c

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

# rsd?
rsd = True
params['rsd'] = rsd

# HOD, Zheng+2009, Kwan+2015
M_cut = 10**13.35 # these constants are taken at the middle of the design, Kwan+15
log_Mcut = np.log10(M_cut)
M1 = 10**13.8
log_M1 = np.log10(M1)
sigma = 0.85
alpha = 1
kappa = 1

# HOD prescription 
design = {'M_cut': M_cut, 'M1': M1, 'sigma': sigma, 'alpha': alpha, 'kappa': kappa}
decorations = {'s': 0, 's_v': 0, 'alpha_c': 0, 's_p': 0, 'A': 0}

# which simulation box are we using?
whichsim = 0

# compute the median halo concentration fit, you dont need to run this if you dont 
# plan on invoking assembly bias decoration. 
avg_c(params, rsd)

# generate galaxy catalogs
galcat.gen_gal_cat(whichsim, design, decorations, params, rsd = rsd)

```

### Output:
The script outputs two files, a central galaxy catalog and a satellite galaxy catalog, in binary format. 
The outputs are stored in a folder `./data_rsd` if the `rsd` flag is True and `./data` if the `rsd` flag is False. For each generalized HOD prescription, the resulting data files are further located in a folder with the name `rockstar_<HOD prescription>_decor_<decorations>`. 

The central and satellite galaxy catalog files are named `halos_gal_cent_<whichsim>` and `halos_gal_sats_<whichsim>`, respectively. Each file contains the 3D position, halo ID and halo mass of the galaxies. 

An example code to read these data files is
```
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
```
where `savedir` is the directory where the data files live in, and `whichsim` is the index of the simulation box. 


## Citation:
If you use this code, please cite [Yuan, Eisenstein & Garrison (2017)](http://adsabs.harvard.edu/doi/10.1093/mnras/stx2032).

## Help 
If you need help with the code, please contact me (Sandy) at sihan.yuan@cfa.harvard.edu. 
