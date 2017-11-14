# GRAND-HOD: GeneRalized ANd Differentiable Halo Occupation Distribution

## Authors:
Sihan (Sandy) Yuan, Daniel Eisenstein & Lehman Garrison

## Introduction:
An Halo Occupation Distribution (HOD) code in Python that is differentiable and incorporates various generalizations to the standard HOD. The code takes a generalized HOD prescription as input and outputs the corresponding mock galaxy catalogs in binary files.

The code is currently written specifically for the Abacus simulations, but the main functionalities can be also easily adapted for other halo catalogs with the appropriate properties. We will update the code soon to support halo catalogs from other cosmological simulations. 

## Usage:

The code does not currently have dependencies other than basic Python packages. 

To install, simply download the script to the directory you want the mock catalogs to live in. If you are not on the Eisenstein Group computer clusters at CfA, you may need to change the `directory` variable to point to the location of the simulation data. Keep the simulation box tag `str(whichsim)` as an unknown variable. 

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
- 'decorations' : dictionary. The HOD generalization parameters. The dictionary may contain the following five parameters:
 - 

### Output:

## Citation:

## Help 
If you need help with the code, please contact me (Sandy) at sihan.yuan@cfa.harvard.edu. 
