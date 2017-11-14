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
- `decorations` : dictionary. The HOD generalization parameters. The dictionary may contain the following five parameters:
  - `s` : float. The satellite profile modulation parameter. Modulates how the radial distribution of satellite galaxies within halos deviate from the radial profile of the halo. Positive value favors satellite galaxies to populate the outskirts of the halo whereas negative value favors satellite galaxy to live near the center of the halo. 
  - `s_v` : float. The satellite velocity bias parameter. Modulates how the satellite galaxy peculiar velocity deviates from that of the local dark matter particle. Positive value favors high peculiar velocity satellite galaxies and vice versa. Note that our implementation preserves the Newton's second law of the satellite galaxies.
  - `alpha_c` : float. The central velocity bias parameter. Modulates the peculiar velocity of the central galaxy. The larger the absolute value of the parameter, the larger the peculiar velocity of the central. The sign of the value does not matter. 
  - `s_p` : float. The perihelion distance modulation parameter. A positive value favors satellite galaxies to have larger distances to the halo center upon their closest approach to the center and vice versa. This can be regarded as a "fancier" satellite profile modulation parameter. 
  - `A` : float. The assembly bias parameter. Introduces the effect of assembly bias. A positive value favors higher concentration halos to host galaxies whereas a negative value favors lower concentration halos to host galaxies. 
  A detailed discussion of these parameters can be found in Yuan et al. in prep. 
 

### Output:

## Citation:

## Help 
If you need help with the code, please contact me (Sandy) at sihan.yuan@cfa.harvard.edu. 
