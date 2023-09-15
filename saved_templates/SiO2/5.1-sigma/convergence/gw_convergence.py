import argparse
import pyparsing as pp
import re
import os

parser = argparse.ArgumentParser()

parser.add_argument('--sigma-bands-list')         # 1000 500 200
parser.add_argument('--screened-cut-list')        # 40 30 20 10
parser.add_argument('--epsilon-bands-list')        # 1000 500 200
parser.add_argument('--kgrid-list')                 # 2x2x2 4x4x4. Only MP type grids. 
parser.add_argument('--infinite-bands-assump')        # 1000
parser.add_argument('--infinite-cutoff-assump')        # 40

# Get arguments. 
sigma_bands_init = 1
screened_cut_init = 1
epsilon_bands_init = 1
kgrid_init = 1
sigma_bands_list = []
screened_cut_list = []
epsilon_bands_list = []
kgrid_list = []

# They are done independently. 

pattern_kgrid = ...+ pp.Literal('K_POINTS') + 'automatic' + ...+ pp.Char(pp.srange('[A-Z]'))[1, ...]
pattern_ecutwfn = ...+ pp.Literal('ecutwfc') + '=' + ...+ pp.Char(pp.srange('[a-zA-Z]')) + pp.Regex(r'.*', re.DOTALL)

# Functions. 
def save_initial():
    pass

def update_sigma_bands(nbnd):
    pass

def update_screened_cutoff(cutoff):
    pass

def update_epsilon_bands(nbnd):
    pass

def update_kgrid(kgrid):
    pass

def run_cmd(ecutwfc):
    pass

def make_initial_dirs():
    pass

def restore_defaults():
    pass

save_initial()

make_initial_dirs()

# Loop over sigma_bands and save.
for nbnd in sigma_bands_list:

    os.chdir('../')  # to 5.1-sigma folder.
    update_sigma_bands(nbnd) 

    # Run epsilon. 
    os.chdir('../4-epsilon')
    run_cmd('./jobscript_epsilon.run')          # Batch and wait with slurm. 

    # Run sigma. 
    os.chdir('../5.1-sigma')
    run_cmd('./jobscript_sigma.run')

    # Copy results. 
    os.system(f'cp ch_convergence.dat ./convergence/sigma_bands/ch_convergence_{nbnd}.dat')
    os.system(f'cp eqp1.dat ./convergence/sigma_bands/eqp1_{nbnd}.dat')
    os.system(f'cp sigma_hp.log ./convergence/sigma_bands/sigma_hp_{nbnd}.log')

    # Switch back to convergence directory. 
    os.chdir('./convergence')

# Loop over screened_cutoff and save. 
for screened_cut in screened_cut_list:

    os.chdir('../')  # to 5.1-sigma folder.
    update_screened_cutoff(screened_cut) 

    # Run epsilon. 
    os.chdir('../4-epsilon')
    run_cmd('./jobscript_epsilon.run')          # Batch and wait with slurm. 

    # Run sigma. 
    os.chdir('../5.1-sigma')
    run_cmd('./jobscript_sigma.run')

    # Copy results. 
    os.system(f'cp ch_convergence.dat ./convergence/screened_cut/ch_convergence_{nbnd}.dat')
    os.system(f'cp eqp1.dat ./convergence/screened_cut/eqp1_{nbnd}.dat')
    os.system(f'cp sigma_hp.log ./convergence/screened_cut/sigma_hp_{nbnd}.log')

    # Switch back to convergence directory. 
    os.chdir('./convergence')


# Loop over epsilon_bands and save.
for nbnd in epsilon_bands_list:

    os.chdir('../')  # to 5.1-sigma folder.
    update_epsilon_bands(nbnd) 

    # Run epsilon. 
    os.chdir('../4-epsilon')
    run_cmd('./jobscript_epsilon.run')          # Batch and wait with slurm. 

    # Run sigma. 
    os.chdir('../5.1-sigma')
    run_cmd('./jobscript_sigma.run')

    # Copy results. 
    os.system(f'cp ch_convergence.dat ./convergence/epsilon_bands/ch_convergence_{nbnd}.dat')
    os.system(f'cp eqp1.dat ./convergence/epsilon_bands/eqp1_{nbnd}.dat')
    os.system(f'cp sigma_hp.log ./convergence/epsilon_bands/sigma_hp_{nbnd}.log')

    # Switch back to convergence directory. 
    os.chdir('./convergence')


# Loop over kgrids and save.
for kgrid in kgrid_list:

    os.chdir('../')  # to 5.1-sigma folder.
    update_kgrid(kgrid)         # Will update the kgrid joscripts with the given grid. 

    # Run 2.1-wfn. 
    os.chdir('../2.1-wfn')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    run_cmd('./jobscript_parabands.run')

    # Run 2.2-wfnq. 
    os.chdir('../2.1-wfn')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')

    # Run epsilon. 
    os.chdir('../4-epsilon')
    run_cmd('./jobscript_epsilon.run')          # Batch and wait with slurm. 

    # Run sigma. 
    os.chdir('../5.1-sigma')
    run_cmd('./jobscript_sigma.run')

    # Copy results. 
    os.system(f'cp ch_convergence.dat ./convergence/kgrid/ch_convergence_{kgrid}.dat')
    os.system(f'cp eqp1.dat ./convergence/kgrid/eqp1_{kgrid}.dat')
    os.system(f'cp sigma_hp.log ./convergence/kgrid/sigma_hp_{kgrid}.log')

    # Switch back to convergence directory. 
    os.chdir('./convergence')


restore_defaults()

# Print time taken. 