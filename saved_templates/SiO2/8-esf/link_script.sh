#!/bin/bash

# dft_force. 
ln -sf ../1.1-scf/pwscf.in ./pwscf.in
ln -sf ../1.1-scf/pwscf.out ./pwscf.out
ln -sf ../2.1-wfn/WFN.h5 ./WFN.h5

# kgrid. 
ln -sf ../2.1-wfn/kgrid.out ./kgrid.out
#This overrides the above step and creates the full kgrid order using kmesh.pl. 
kmesh.pl 2 2 2 > kgrid.out          

# eqp (eqp).
ln -sf ../7-absorption/eqp.dat ./eqp.dat

# elph (g).
ln -sf ../1.3-epw/epw.in ./epw.in
file_idx=0
file_name=""
for file in $(find ../1.1-scf/tmp -name SiO2_elph_*)
do
    ((file_idx++))
    file_name="elph_"$file_idx".h5"
    ln -sf $file ./$file_name
done

# kernel (K).
ln -sf ../7-absorption/hbse.h5 ./hbse.h5

# eigenvectors and eigenvalues(A).
ln -sf ../7-absorption/eigenvectors.h5 ./eigenvectors.h5
ln -sf ../7-absorption/eigenvalues.dat ./eigenvalues.dat
