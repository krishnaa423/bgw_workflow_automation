import argparse
import pyparsing as pp
import re
import os

parser = argparse.ArgumentParser()

parser.add_argument('--coarse-kgrid-list')         # Only MP grid types 2x2x2 4x4x4 ...
parser.add_argument('--coarse-cond-list')         # 2 4 5 10
parser.add_argument('--coarse-val-list')         # 2 4 5 10

parser.add_argument('--fine-kgrid-list')         # Only MP grid types 2x2x2 4x4x4 ...
parser.add_argument('--fine-cond-list')         # 2 4 5 10
parser.add_argument('--fine-val-list')         # 2 4 5 10

# Get arguments. 
coarse_kgrid_list = []
coarse_cond_list = []
coarse_val_list = []

fine_kgrid_list = []
fine_cond_list = []
fine_val_list = []

# They are done independently. 
pattern_kgrid = ...+ pp.Literal('K_POINTS') + 'automatic' + ...+ pp.Char(pp.srange('[A-Z]'))[1, ...]
pattern_ecutwfn = ...+ pp.Literal('ecutwfc') + '=' + ...+ pp.Char(pp.srange('[a-zA-Z]')) + pp.Regex(r'.*', re.DOTALL)



# Functions. 
def save_initial():
    pass

def update_coarse_kgrid(kgrid):
    pass

def update_coarse_bands(cond_bnd, val_bnd):
    pass

def update_fine_kgrid(kgrid):
    pass

def update_fine_bands(cond_bnd, val_bnd):
    pass

def restore_defaults():
    pass

def make_initial_dirs():
    pass

def run_cmd(cmd):
    pass

save_initial()

make_initial_dirs()

# Loop over kgrid and save.
for kgrid in coarse_kgrid_list:

    os.chdir('../')  # to 7.1/absorption folder.
    update_coarse_kgrid(kgrid)          # updates jobscript_kgrid.run in 2.1-wfn. 

    # Run 2.1-wfn.
    os.chdir('../2.1-wfn')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    run_cmd('./jobscript_parabands.run')

    # Run 2.2-wfnq.
    os.chdir('../2.2-wfnq')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')

    
    # Run 4-epsilon. 
    os.chdir('../4-epsilon')
    run_cmd('./jobscript_epsilon.run')


    # Run 5.1-sigma. 
    os.chdir('../5.1-sigma')
    run_cmd('./jobscript_sigma.run')

    # Run 6-kernel. 
    os.chdir('../6-kernel')
    run_cmd('./jobscript_kernel.run')

    # Run 7.1-absorption.
    os.chdir('../7.1-absorption')
    run_cmd('./jobscript_absorption.run')

    # Copy results. 
    os.system(f'cp eigenvalues.dat ./convergence/coarse_kgrid/eigenvalues_{kgrid}.dat')
    os.system(f'cp absorption_eh.dat ./convergence/coarse_kgrid/absorption_eh_{kgrid}.dat')

    # Switch back to convergence directory. 
    os.chdir('./convergence')

# Loop over fine k-grid and save.
for kgrid in fine_kgrid_list:

    os.chdir('../')  # to 7.1/absorption folder.
    update_coarse_kgrid(kgrid)          # updates jobscript_kgrid.run in 3.1-wfn_fi. 

    # Run 3.1-wfn_fi.
    os.chdir('../3.1-wfn_fi')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    run_cmd('./jobscript_parabands.run')

    # Run 3.2-wfnq_fi.
    os.chdir('../3.2-wfnq_fi')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')

    # Run 7.1-absorption.
    os.chdir('../7.1-absorption')
    run_cmd('./jobscript_absorption.run')

    # Copy results. 
    os.system(f'cp eigenvalues.dat ./convergence/fine_kgrid/eigenvalues_{kgrid}.dat')
    os.system(f'cp absorption_eh.dat ./convergence/fine_kgrid/absorption_eh_{kgrid}.dat')

    # Switch back to convergence directory. 
    os.chdir('./convergence')

# Loop over coarse bands and save.
for cond_bnd, val_band in zip(coarse_cond_list, coarse_val_list):

    os.chdir('../')  # to 7.1/absorption folder.
    update_coarse_bands(cond_bnd, val_band)          # updates kernel.inp.

    # Run 6-kernel. 
    os.chdir('../6-kernel')
    run_cmd('./jobscript_kernel.run')

    # Run 7.1-absorption.
    os.chdir('../7.1-absorption')
    run_cmd('./jobscript_absorption.run')

    # Copy results. 
    os.system(f'cp eigenvalues.dat ./convergence/coarse_bands/eigenvalues_{cond_bnd}_{val_band}.dat')
    os.system(f'cp absorption_eh.dat ./convergence/coarse_bands/absorption_eh_{cond_bnd}_{val_band}.dat')

    # Switch back to convergence directory. 
    os.chdir('./convergence')


# Loop over fine bands and save.
for cond_bnd, val_band in zip(fine_cond_list, fine_val_list):

    os.chdir('../')  # to 7.1/absorption folder.
    update_fine_bands(cond_bnd, val_band)          # updates absorption.inp.

    # Run 7.1-absorption.
    os.chdir('../7.1-absorption')
    run_cmd('./jobscript_absorption.run')

    # Copy results. 
    os.system(f'cp eigenvalues.dat ./convergence/fine_bands/eigenvalues_{cond_bnd}_{val_band}.dat')
    os.system(f'cp absorption_eh.dat ./convergence/fine_bands/absorption_eh_{cond_bnd}_{val_band}.dat')

    # Switch back to convergence directory. 
    os.chdir('./convergence')

restore_defaults()

# Print time taken. 