#!/bin/bash 


# Copy for absorption folder

# From wfn folder. 
ln -sf ../2.1-wfn/WFN.h5 ./WFN_co.h5
ln -sf ../2.1-wfn/WFN.h5 ./WFN_fi.h5
ln -sf ../2.2-wfnq/WFNq.h5 ./WFNq_fi.h5


# From epsilon folder. 
ln -sf ../4-epsilon/eps0mat.h5 ./eps0mat.h5
ln -sf ../4-epsilon/epsmat.h5 ./epsmat.h5

# From sigma folder. 
ln -sf ../5.1-sigma/eqp1.dat ./eqp_co.dat

# From kernel folder. 
ln -sf ../6-kernel/bsemat.h5 ./bsemat.h5
