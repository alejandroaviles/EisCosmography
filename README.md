# EisCosmography

Alejandro Aviles 
(ABACUS-CINVESTAV)
aviles@ciencias.unam.mx
avilescervantes@gmail.com

In collaboration with Jaime Klapp and Orlando Luongo

CosmoMC module for the Eis cosmography method. + Simulated supernovae catalogues + bias statistics. 

The file cosmography.f90 has the module to the code CosmoMC to fit the statefinders parameters in different methods in cosmography. As described in arXiv:160X.XXXX (or see EisCosmography.pdf in this repository). This module does not depend on any other CosmoMC module or externals, thus it should work on any CosmoMC version.

The dirs /*SimulatedData have all the simulated catalogues used in the paper arXiv:1606.XXXX. Further statistics than those presented in the paper can be found in dir /statistics

A long table that extends Table 1 in the paper, for the 100 dispersed catalogs fits, can be found in file longtable_statefinders.pdf, as well as an equivalent table for the eis parameters (longtable_eis.pdf)
