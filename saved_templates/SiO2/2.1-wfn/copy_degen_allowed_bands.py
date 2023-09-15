
import pyparsing as pp
import re
import os
import numpy as np

# Debugging. 
# os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/2.1-wfn')


# Proposed max bands, val bands and conduction bands. 
eps_bands = 1000
sigma_bands = 1000
sigma_min = 38
sigma_max = 70
coarse_val_bands = 10
fine_val_bands = 10
coarse_cond_bands = 22
fine_cond_bands = 22

# File names. 
epsilon = '../4-epsilon/epsilon.inp'
sigma = '../5.1-sigma/sigma.inp'
kernel = '../6-kernel/kernel.inp'
absorption = '../7-absorption/absorption.inp'


# Functions. 
def replace_eps(bnd):
    with open(epsilon, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_bands\s+\d+', f'number_bands {bnd}', txt)
    with open(epsilon, 'w') as f: txt = f.write(new_txt)

def replace_sigma(bnd, min, max):
    
    # Sigma bands. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_bands\s+\d+', f'number_bands {bnd}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

    # min. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'band_index_min\s+\d+', f'band_index_min {min}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

    # max. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'band_index_max\s+\d+', f'band_index_max {max}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

def replace_kernel(vbc, cbc):
    # vbc
    with open(kernel, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands\s+\d+', f'number_val_bands {vbc}', txt)
    with open(kernel, 'w') as f: txt = f.write(new_txt)

    # cbc
    with open(kernel, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands\s+\d+', f'number_cond_bands {cbc}', txt)
    with open(kernel, 'w') as f: txt = f.write(new_txt)

def replace_absorption(vbf, vbc, cbf, cbc):
    # vbc
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands_coarse\s+\d+', f'number_val_bands_coarse {vbc}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # vbf
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands_fine\s+\d+', f'number_val_bands_fine {vbf}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # cbc
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands_coarse\s+\d+', f'number_cond_bands_coarse {cbc}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # cbf
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands_fine\s+\d+', f'number_cond_bands_fine {cbf}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)


os.system('degeneracy_check.x WFN > degen.log')
inner_array = []
outer_array = []

with open('degen.log', 'r') as f:

    lines = f.readlines()

    lines = [line for line in lines if not (line[0].isalpha() or line[0] == '=')]

    with open('filtered_degen.log', 'w') as g:
        g.writelines(lines[1:])


with open('filtered_degen.log', 'r') as f:

    lines = f.readlines()
    for line in lines:
        try:
            val = float(line)
            inner_array.append(val)
        except Exception as e:
            outer_array.append(inner_array)
            inner_array = []
            continue

    outer_array.append(inner_array)

print(f'Outer array length: {len(outer_array)}')

max_allowed_eps_bands = int([bnd for bnd in outer_array[0] if bnd <= eps_bands][-1])
max_allowed_sigma_bands = int([bnd for bnd in outer_array[0] if bnd <= sigma_bands][-1])
allowed_sigma_min = int([bnd for bnd in outer_array[0] if bnd >= sigma_min][0])
allowed_sigma_max = int([bnd for bnd in outer_array[0] if bnd <= sigma_max][-1])
max_allowed_val_bands_fine = int([bnd for bnd in outer_array[1] if bnd <= fine_val_bands][-1])
max_allowed_val_bands_coarse = int([bnd for bnd in outer_array[1] if bnd <= coarse_val_bands][-1])
max_allowed_cond_bands_fine = int([bnd for bnd in outer_array[2] if bnd <= fine_cond_bands][-1])
max_allowed_cond_bands_coarse = int([bnd for bnd in outer_array[2] if bnd <= coarse_cond_bands][-1])

# replace files. 
replace_eps(max_allowed_eps_bands)
replace_sigma(max_allowed_sigma_bands, allowed_sigma_min, allowed_sigma_max)
replace_kernel(max_allowed_val_bands_coarse, max_allowed_cond_bands_coarse)
replace_absorption(max_allowed_val_bands_fine, max_allowed_val_bands_coarse, max_allowed_cond_bands_fine, max_allowed_cond_bands_coarse)
