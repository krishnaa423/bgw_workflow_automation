#!/bin/bash 

# Copy for sigma folder

# From wfn folder. 
ln -sf ../2.1-wfn/vxc.dat ./vxc.dat
ln -sf ../2.1-wfn/RHO ./RHO
ln -sf ../2.1-wfn/WFN.h5 ./WFN_inner.h5


# From epsilon folder. 
ln -sf ../4-epsilon/eps0mat.h5 ./eps0mat.h5
ln -sf ../4-epsilon/epsmat.h5 ./epsmat.h5
