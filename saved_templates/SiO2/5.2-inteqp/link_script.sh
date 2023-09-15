#!/bin/bash

# From wfn folder. Coarse wavefunctions. 
ln -sf ../2.1-wfn/WFN.h5 ./WFN_co.h5

# From bands folder. Fine wavefunctions. 
ln -sf ../1.5-bands/WFN.h5 ./WFN_fi.h5

# From sigma folder. Quasiparticle energy corrections. 
ln -sf ../5.1-sigma/eqp1.dat ./eqp_co.dat 
