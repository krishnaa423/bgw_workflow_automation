import os
from ase import Atoms
from ase.io import espresso, read
import numpy as np
import pyparsing as pp
import io 
import json 
import re


# For debugging. 
# os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/Frontera/SiO2_v3_esf')

# region: Global variables.

# region: Units.
bohr2A = 0.529177249  # bohr 2 A conversion. 
ry2eV = 13.6057039763 # Rydberg to eV conversion. 
# endregion.

# region: Create structure object from CIF file. 
struct = read('SiO2_2x1.cif')
# endregion.


# region: Create input dict. Read from file if needed. Write to file if needed. 
input= dict()
input = {
    # region: structure.
    'struct_formula': 'SiO2',
    'num_electrons_uc': 96,
    'num_val_bands_uc': 48,
    'cell': struct.cell.cellpar().tolist(), # a, b, c, \alpha, \beta, \gamma.
    'positions': struct.get_positions().tolist(),          # Positions in A. 
    'pbc': struct.get_pbc().tolist(),
    # endregion. 

    # region: scf
    'pwscf_wfc_cutoff': 60.0,
    'pwscf_kpts': (2, 2, 2),
    'pw_pseudopotentials': {
        'Si': 'Si_ONCV_PBE_sr.upf',
        'O': 'O_ONCV_PBE_sr.upf'
    },
    'pseudo_folder': 'qe_pseudo',
    # endregion. 

    # region: ph
    'ph_bands': 70,         # vb is at 48. 70 gives 22 valence bands. 
    'ph_kpts': (2, 2, 2),
    'ph_nirk': 8,
    'ph_scf_tol': '1.0d-6',
    # endregion. 

    # region bands.
    'bands_pwbands': 100,
    'bands_parabands': 100,
    'bands_plot_kpath_string': 
'''K_POINTS crystal
        201
      0.0000000000    0.0000000000    0.0000000000    1.0
      0.0000000000    0.0089285714    0.0000000000    1.0
      0.0000000000    0.0178571429    0.0000000000    1.0
      0.0000000000    0.0267857143    0.0000000000    1.0
      0.0000000000    0.0357142857    0.0000000000    1.0
      0.0000000000    0.0446428571    0.0000000000    1.0
      0.0000000000    0.0535714286    0.0000000000    1.0
      0.0000000000    0.0625000000    0.0000000000    1.0
      0.0000000000    0.0714285714    0.0000000000    1.0
      0.0000000000    0.0803571429    0.0000000000    1.0
      0.0000000000    0.0892857143    0.0000000000    1.0
      0.0000000000    0.0982142857    0.0000000000    1.0
      0.0000000000    0.1071428571    0.0000000000    1.0
      0.0000000000    0.1160714286    0.0000000000    1.0
      0.0000000000    0.1250000000    0.0000000000    1.0
      0.0000000000    0.1339285714    0.0000000000    1.0
      0.0000000000    0.1428571429    0.0000000000    1.0
      0.0000000000    0.1517857143    0.0000000000    1.0
      0.0000000000    0.1607142857    0.0000000000    1.0
      0.0000000000    0.1696428571    0.0000000000    1.0
      0.0000000000    0.1785714286    0.0000000000    1.0
      0.0000000000    0.1875000000    0.0000000000    1.0
      0.0000000000    0.1964285714    0.0000000000    1.0
      0.0000000000    0.2053571429    0.0000000000    1.0
      0.0000000000    0.2142857143    0.0000000000    1.0
      0.0000000000    0.2232142857    0.0000000000    1.0
      0.0000000000    0.2321428571    0.0000000000    1.0
      0.0000000000    0.2410714286    0.0000000000    1.0
      0.0000000000    0.2500000000    0.0000000000    1.0
      0.0000000000    0.2589285714    0.0000000000    1.0
      0.0000000000    0.2678571429    0.0000000000    1.0
      0.0000000000    0.2767857143    0.0000000000    1.0
      0.0000000000    0.2857142857    0.0000000000    1.0
      0.0000000000    0.2946428571    0.0000000000    1.0
      0.0000000000    0.3035714286    0.0000000000    1.0
      0.0000000000    0.3125000000    0.0000000000    1.0
      0.0000000000    0.3214285714    0.0000000000    1.0
      0.0000000000    0.3303571429    0.0000000000    1.0
      0.0000000000    0.3392857143    0.0000000000    1.0
      0.0000000000    0.3482142857    0.0000000000    1.0
      0.0000000000    0.3571428571    0.0000000000    1.0
      0.0000000000    0.3660714286    0.0000000000    1.0
      0.0000000000    0.3750000000    0.0000000000    1.0
      0.0000000000    0.3839285714    0.0000000000    1.0
      0.0000000000    0.3928571429    0.0000000000    1.0
      0.0000000000    0.4017857143    0.0000000000    1.0
      0.0000000000    0.4107142857    0.0000000000    1.0
      0.0000000000    0.4196428571    0.0000000000    1.0
      0.0000000000    0.4285714286    0.0000000000    1.0
      0.0000000000    0.4375000000    0.0000000000    1.0
      0.0000000000    0.4464285714    0.0000000000    1.0
      0.0000000000    0.4553571429    0.0000000000    1.0
      0.0000000000    0.4642857143    0.0000000000    1.0
      0.0000000000    0.4732142857    0.0000000000    1.0
      0.0000000000    0.4821428571    0.0000000000    1.0
      0.0000000000    0.4910714286    0.0000000000    1.0
      0.0000000000    0.5000000000    0.0000000000    1.0
      0.0104166667    0.4895833333    0.0000000000    1.0
      0.0208333333    0.4791666667    0.0000000000    1.0
      0.0312500000    0.4687500000    0.0000000000    1.0
      0.0416666667    0.4583333333    0.0000000000    1.0
      0.0520833333    0.4479166667    0.0000000000    1.0
      0.0625000000    0.4375000000    0.0000000000    1.0
      0.0729166667    0.4270833333    0.0000000000    1.0
      0.0833333333    0.4166666667    0.0000000000    1.0
      0.0937500000    0.4062500000    0.0000000000    1.0
      0.1041666667    0.3958333333    0.0000000000    1.0
      0.1145833333    0.3854166667    0.0000000000    1.0
      0.1250000000    0.3750000000    0.0000000000    1.0
      0.1354166667    0.3645833333    0.0000000000    1.0
      0.1458333333    0.3541666667    0.0000000000    1.0
      0.1562500000    0.3437500000    0.0000000000    1.0
      0.1666666667    0.3333333333    0.0000000000    1.0
      0.1770833333    0.3229166667    0.0000000000    1.0
      0.1875000000    0.3125000000    0.0000000000    1.0
      0.1979166667    0.3020833333    0.0000000000    1.0
      0.2083333333    0.2916666667    0.0000000000    1.0
      0.2187500000    0.2812500000    0.0000000000    1.0
      0.2291666667    0.2708333333    0.0000000000    1.0
      0.2395833333    0.2604166667    0.0000000000    1.0
      0.2500000000    0.2500000000    0.0000000000    1.0
      0.2604166667    0.2395833333    0.0000000000    1.0
      0.2708333333    0.2291666667    0.0000000000    1.0
      0.2812500000    0.2187500000    0.0000000000    1.0
      0.2916666667    0.2083333333    0.0000000000    1.0
      0.3020833333    0.1979166667    0.0000000000    1.0
      0.3125000000    0.1875000000    0.0000000000    1.0
      0.3229166667    0.1770833333    0.0000000000    1.0
      0.3333333333    0.1666666667    0.0000000000    1.0
      0.3437500000    0.1562500000    0.0000000000    1.0
      0.3541666667    0.1458333333    0.0000000000    1.0
      0.3645833333    0.1354166667    0.0000000000    1.0
      0.3750000000    0.1250000000    0.0000000000    1.0
      0.3854166667    0.1145833333    0.0000000000    1.0
      0.3958333333    0.1041666667    0.0000000000    1.0
      0.4062500000    0.0937500000    0.0000000000    1.0
      0.4166666667    0.0833333333    0.0000000000    1.0
      0.4270833333    0.0729166667    0.0000000000    1.0
      0.4375000000    0.0625000000    0.0000000000    1.0
      0.4479166667    0.0520833333    0.0000000000    1.0
      0.4583333333    0.0416666667    0.0000000000    1.0
      0.4687500000    0.0312500000    0.0000000000    1.0
      0.4791666667    0.0208333333    0.0000000000    1.0
      0.4895833333    0.0104166667    0.0000000000    1.0
      0.5000000000    0.0000000000    0.0000000000    1.0
      0.4903846154    0.0000000000    0.0096153846    1.0
      0.4807692308    0.0000000000    0.0192307692    1.0
      0.4711538462    0.0000000000    0.0288461538    1.0
      0.4615384615    0.0000000000    0.0384615385    1.0
      0.4519230769    0.0000000000    0.0480769231    1.0
      0.4423076923    0.0000000000    0.0576923077    1.0
      0.4326923077    0.0000000000    0.0673076923    1.0
      0.4230769231    0.0000000000    0.0769230769    1.0
      0.4134615385    0.0000000000    0.0865384615    1.0
      0.4038461538    0.0000000000    0.0961538462    1.0
      0.3942307692    0.0000000000    0.1057692308    1.0
      0.3846153846    0.0000000000    0.1153846154    1.0
      0.3750000000    0.0000000000    0.1250000000    1.0
      0.3653846154    0.0000000000    0.1346153846    1.0
      0.3557692308    0.0000000000    0.1442307692    1.0
      0.3461538462    0.0000000000    0.1538461538    1.0
      0.3365384615    0.0000000000    0.1634615385    1.0
      0.3269230769    0.0000000000    0.1730769231    1.0
      0.3173076923    0.0000000000    0.1826923077    1.0
      0.3076923077    0.0000000000    0.1923076923    1.0
      0.2980769231    0.0000000000    0.2019230769    1.0
      0.2884615385    0.0000000000    0.2115384615    1.0
      0.2788461538    0.0000000000    0.2211538462    1.0
      0.2692307692    0.0000000000    0.2307692308    1.0
      0.2596153846    0.0000000000    0.2403846154    1.0
      0.2500000000    0.0000000000    0.2500000000    1.0
      0.2403846154    0.0000000000    0.2596153846    1.0
      0.2307692308    0.0000000000    0.2692307692    1.0
      0.2211538462    0.0000000000    0.2788461538    1.0
      0.2115384615    0.0000000000    0.2884615385    1.0
      0.2019230769    0.0000000000    0.2980769231    1.0
      0.1923076923    0.0000000000    0.3076923077    1.0
      0.1826923077    0.0000000000    0.3173076923    1.0
      0.1730769231    0.0000000000    0.3269230769    1.0
      0.1634615385    0.0000000000    0.3365384615    1.0
      0.1538461538    0.0000000000    0.3461538462    1.0
      0.1442307692    0.0000000000    0.3557692308    1.0
      0.1346153846    0.0000000000    0.3653846154    1.0
      0.1250000000    0.0000000000    0.3750000000    1.0
      0.1153846154    0.0000000000    0.3846153846    1.0
      0.1057692308    0.0000000000    0.3942307692    1.0
      0.0961538462    0.0000000000    0.4038461538    1.0
      0.0865384615    0.0000000000    0.4134615385    1.0
      0.0769230769    0.0000000000    0.4230769231    1.0
      0.0673076923    0.0000000000    0.4326923077    1.0
      0.0576923077    0.0000000000    0.4423076923    1.0
      0.0480769231    0.0000000000    0.4519230769    1.0
      0.0384615385    0.0000000000    0.4615384615    1.0
      0.0288461538    0.0000000000    0.4711538462    1.0
      0.0192307692    0.0000000000    0.4807692308    1.0
      0.0096153846    0.0000000000    0.4903846154    1.0
      0.0000000000    0.0000000000    0.5000000000    1.0
      0.0000000000    0.0000000000    0.4886363636    1.0
      0.0000000000    0.0000000000    0.4772727273    1.0
      0.0000000000    0.0000000000    0.4659090909    1.0
      0.0000000000    0.0000000000    0.4545454545    1.0
      0.0000000000    0.0000000000    0.4431818182    1.0
      0.0000000000    0.0000000000    0.4318181818    1.0
      0.0000000000    0.0000000000    0.4204545455    1.0
      0.0000000000    0.0000000000    0.4090909091    1.0
      0.0000000000    0.0000000000    0.3977272727    1.0
      0.0000000000    0.0000000000    0.3863636364    1.0
      0.0000000000    0.0000000000    0.3750000000    1.0
      0.0000000000    0.0000000000    0.3636363636    1.0
      0.0000000000    0.0000000000    0.3522727273    1.0
      0.0000000000    0.0000000000    0.3409090909    1.0
      0.0000000000    0.0000000000    0.3295454545    1.0
      0.0000000000    0.0000000000    0.3181818182    1.0
      0.0000000000    0.0000000000    0.3068181818    1.0
      0.0000000000    0.0000000000    0.2954545455    1.0
      0.0000000000    0.0000000000    0.2840909091    1.0
      0.0000000000    0.0000000000    0.2727272727    1.0
      0.0000000000    0.0000000000    0.2613636364    1.0
      0.0000000000    0.0000000000    0.2500000000    1.0
      0.0000000000    0.0000000000    0.2386363636    1.0
      0.0000000000    0.0000000000    0.2272727273    1.0
      0.0000000000    0.0000000000    0.2159090909    1.0
      0.0000000000    0.0000000000    0.2045454545    1.0
      0.0000000000    0.0000000000    0.1931818182    1.0
      0.0000000000    0.0000000000    0.1818181818    1.0
      0.0000000000    0.0000000000    0.1704545455    1.0
      0.0000000000    0.0000000000    0.1590909091    1.0
      0.0000000000    0.0000000000    0.1477272727    1.0
      0.0000000000    0.0000000000    0.1363636364    1.0
      0.0000000000    0.0000000000    0.1250000000    1.0
      0.0000000000    0.0000000000    0.1136363636    1.0
      0.0000000000    0.0000000000    0.1022727273    1.0
      0.0000000000    0.0000000000    0.0909090909    1.0
      0.0000000000    0.0000000000    0.0795454545    1.0
      0.0000000000    0.0000000000    0.0681818182    1.0
      0.0000000000    0.0000000000    0.0568181818    1.0
      0.0000000000    0.0000000000    0.0454545455    1.0
      0.0000000000    0.0000000000    0.0340909091    1.0
      0.0000000000    0.0000000000    0.0227272727    1.0
      0.0000000000    0.0000000000    0.0113636364    1.0
      0.0000000000    0.0000000000    0.0000000000    1.0
''',
    # endregion. 

    # region: dos
    'dos_bands': 100,
    'dos_kpts': (12, 12, 12),
    'dos_e_min': -5.0,
    'dos_e_max': 15.0,
    # endregion

    # region: wfn
    'wfn_pw_bands': 100,
    'wfn_parabands': 1200,
    'wfnq_pw_bands': 100,
    'wfnq_parabands': 100,
    'wfn_kpts': (2, 2, 2),
    'wfn_qshift': (0.000, 0.000, 0.001),
    'wfn_fi_kpts': (2, 2, 2),
    'wfn_vxc_diag_min': 1,
    'wfn_vxc_diag_max': 100,
    'wfn_vxc_offdiag_min': 0,
    'wfn_vxc_offdiag_max': 0,
    # endregion. 

    # region: epsilon. 
    'epsilon_cutoff': 20.0,
    'epsilon_bands': 1000,
    # endregion.

    # region: sigma.
    'sigma_bands': 1000,
    'sigma_bands_min': 38,          # 11 valence bands.     
    'sigma_bands_max': 70,          # 22 conduction bands. 
    # endregion. 

    # region: kernel. 
    'coarse_cond_bands': 22,
    'coarse_val_bands': 10,
    # endregion. 

    # region: absorption. 
    'fine_cond_bands': 22,
    'fine_val_bands': 10,
    'max_eigenvectors':5,
    'plotxct_state': 1,
    # endregion. 

    # region: esf. 
    'esf_xct_state': 0,
    'esf_tol': 9e-2,
    'esf_max_steps': 5,
    # endregion. 

    # region: SLURM. 
    # Below options are written with Frontera in mind. Will have to change for other clusters. 
    'sched_type': 'slurm', # Can be 'WSL' or 'slurm'
    'slurm_account_flag': '#SBATCH --account=PHY20032',
    'slurm_partition_flag': '#SBATCH --partition',
    'slurm_nodes_flag': '#SBATCH --nodes',
    'slurm_ntasks_flag': '#SBATCH --ntasks',
    'slurm_time_flag': '#SBATCH --time',
    'slurm_mail_user_flag': '#SBATCH --mail-user=krishnaa.vadivel@yale.edu',
    'slurm_mail_type_flag': '#SBATCH --mail-type=ALL',
    'slurm_mpi_type': 'ibrun',  # Can be ibrun in Frontera and srun in perlmutter. 
    # endregion. 
}

# Add key values that depend on other keys. 
# region: SLURM.
input['slurm_job_flag'] = f'#SBATCH --job-name={input["struct_formula"]}'
# endregion. 


# Read from file flag. Overrides the above input dict. 
def write_input_json():
    with open('input.json', 'w') as f:
        f.write(json.dumps(input))  

read_input_from_file = False
fname = 'input.json'        # Name of file where input parameters are recorded. 
if read_input_from_file == True:
    # Used to override above input dict with one from file.
    with open(fname, 'r') as f: 
        input = json.loads(f.read())

else:
    # Else write the above created dict to file. 
    write_input_json()
# endregion. 
# endregion


# region: Global functions.
def mkdir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def touch(file_name):
    if not os.path.exists(file_name):
        with open(file_name, 'w') as f:
            pass
# endregion


# Generate folders. 

# region: 1.1-scf
# Switch to current directory.
dir_name = '1.1-scf'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script and other init.
# Make pseudo dir.
mkdir('pseudo')
os.chdir('pseudo')
pseudo_dir = f'../../{input["pseudo_folder"]}'
for val in input['pw_pseudopotentials'].values():
    pseudo_file = pseudo_dir + '/' + val
    os.system(f'cp {pseudo_file} {val}')
os.chdir('../')
# endregion

# region: Input files.
# pwscf.in 
fname = 'pwscf.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'scf',
                'tprnfor': True,
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                # 'nosym':  True              
            },
            'electrons': {
                'diagonalization': 'paro'
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['pwscf_kpts']         # MP grid.
    )
# endregion

# region: Jobscripts.
# jobscript_pwscf.run
fname = 'jobscript_pwscf.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


module load impi
{input['slurm_mpi_type']} pw.x < pwscf.in &> pwscf.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

pw.x < pwscf.in &> pwscf.out
'''
        )
# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.1-scf/* ./')

# endregion

# Go back to parent directory. 
os.chdir('../')  
# endregion

# region: 1.2-ph
# Switch to current directory.
dir_name = '1.2-ph'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands'
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['ph_bands'],
                # 'nosym': True       
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['ph_kpts']         # MP grid will be replaced by what we get after wfn. 
    )

# ph.in
fname = 'ph.in'
with open(fname, 'w') as f:
    f.write(
f'''--
&inputph
  prefix   = '{input['struct_formula']}'
  outdir = './tmp'
  fildvscf = 'dvscf'
  ldisp    = .true.
  fildyn   = '{input['struct_formula']}.dyn'
  nq1={input['ph_kpts'][0]},
  nq2={input['ph_kpts'][1]},
  nq3={input['ph_kpts'][2]},
  tr2_ph   =  {input['ph_scf_tol']}
 /
 '''
    )

# q2r.in
fname = 'q2r.in'
with open(fname, 'w') as f:
    f.write(
f'''&INPUT
  fildyn = '{input['struct_formula']}.dyn'
  zasr = 'crystal'
  flfrc = '{input['struct_formula']}.fc'
/
'''
    )

# matdyn.in
fname = 'matdyn.in'
with open(fname, 'w') as f:
    f.write(
f'''&INPUT
  flfrc = '{input['struct_formula']}.fc'
  flfrq = '{input['struct_formula']}.freq'
  fleig = '{input['struct_formula']}.eig'
  flvec = '{input['struct_formula']}.modes'
/
'''
    )

# matdyn_bands.in
fname = 'matdyn_bands.in'
append_string = '\n'.join(input['bands_plot_kpath_string'].splitlines()[1:])
with open(fname, 'w') as f:
    f.write(
f'''&INPUT
  flfrc = '{input['struct_formula']}.fc'
  flfrq = '{input['struct_formula']}.freq'
  fleig = '{input['struct_formula']}.eig'
  flvec = '{input['struct_formula']}.modes'
  q_in_band_form = .true.
/
{append_string}
'''
    )

# create_epw_save.py. For epw save folder. 
epw_save_name = 'create_epw_save.py'
with open(epw_save_name, 'w') as f:
    f.write(
f'''
#
# Post-processing script QE --> EPW
# 14/07/2015 - Samuel Ponce
#

from builtins import input
import numpy as np
import os

# Enter the number of irr. q-points
# prefix = input('Enter the prefix used for PH calculations (e.g. diam)\\n')
prefix = '{input['struct_formula']}'
outdir = './tmp'

# Enter the number of irr. q-points
nqpt = {input['ph_nirk']}

# try:
#   nqpt = int(nqpt)
# except ValueError:
#   raise Exception('The value you enter is not an integer!')

if not os.path.exists('save'):
        os.mkdir('save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp ' + outdir + '/_ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r ' + outdir + '/_ph0/'+prefix+'.phsave save/')

    # # KV: Maybe. 
    # os.system('cp '+prefix+'.dyn'+' save/'+prefix+'.dyn_q'+label)
  else:
    os.system('cp ' + outdir + '/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('rm ' + outdir + '/_ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )

'''
    )

# copy_kgrid_pwbands_ph_matdyn.py
fname = 'copy_kgrid_matdyn.py'
copy_kgrid_matdyn_text = \
'''
import pyparsing as pp
import re


# Copy to pwbands.in
pattern = ... + pp.Literal('K_POINTS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('./pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: 
    subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: 
    f.write(new_text)


# Copy to ph.in and matdyn.in. 
pattern = ... + pp.Literal('/\\n')  + pp.Regex(r'.*', re.DOTALL)

ph_parsed_text = pattern.parse_file('ph.in').as_list()
matdyn_parsed_text = pattern.parse_file('matdyn.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: subst_text = ''.join(f.readlines()[1:])

ph_new_text = ''.join([
    ph_parsed_text[0],
    ph_parsed_text[1],
    subst_text
])

matdyn_new_text = ''.join([
    matdyn_parsed_text[0],
    matdyn_parsed_text[1],
    subst_text
])

# with open('ph.in', 'w') as f: f.write(ph_new_text)        # Commented out. Not writing to ph.in for now. 
with open('matdyn.in', 'w') as f: f.write(matdyn_new_text)

# Copy number of irreducible k-points to create_epw_save.py. 
pattern = ...+ pp.Literal('nqpt') + '=' + pp.Word(pp.nums) + pp.Regex(r'.*', re.DOTALL)
parsed_text = pattern.parse_file('./create_epw_save.py')

nirk = len(open('kgrid.out', 'r').readlines()) - 2
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    parsed_text[2],
    str(nirk),
    parsed_text[4]
])
open('create_epw_save.py', 'w').write(new_text)
'''
with open(fname, 'w') as f: f.write(copy_kgrid_matdyn_text)

# endregion

# region: Jobscripts.
# jobscript_kgrid.run
fname = 'jobscript_kgrid.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
{input['slurm_mpi_type']} kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid_matdyn.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid_matdyn.py
'''
        )

# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in &> pwbands.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in &> pwbands.out
'''
        )

# jobscript_ph.run
fname = 'jobscript_ph.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
rm -r ./save
{input['slurm_mpi_type']} ph.x -nk {input['ph_kpts'][0]*input['ph_kpts'][1]*input['ph_kpts'][2]} < ph.in &> ph.out
python create_epw_save.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
rm -r ./save
mpirun -n {input['ph_kpts'][0]*input['ph_kpts'][1]*input['ph_kpts'][2]} ph.x -nk {input['ph_kpts'][0]*input['ph_kpts'][1]*input['ph_kpts'][2]} < ph.in &> ph.out
python create_epw_save.py
'''
        )  

# jobscript_q2r.run
fname = 'jobscript_q2r.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} q2r.x < q2r.in &> q2r.out
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
q2r.x < q2r.in &> q2r.out
'''
        )  

# jobscript_matdyn.run
fname = 'jobscript_matdyn.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} matdyn.x < matdyn.in &> matdyn.out
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
matdyn.x < matdyn.in > matdyn.out
'''
        )  

# jobscript_matdyn_bands.run
fname = 'jobscript_matdyn_bands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} matdyn.x < matdyn_bands.in &> matdyn_bands.out
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
matdyn.x < matdyn_bands.in &> matdyn_bands.out
'''
        )  

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.2-ph/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../') 
# endregion

# region: 1.3-epw
# Switch to current directory.
dir_name = '1.3-epw'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# epw.in
fname = 'epw.in'
with open(fname, 'w') as f:
    f.write(
f'''--
&inputepw
  prefix      = '{input['struct_formula']}'
  outdir      = './tmp'

  elph        = .true.
  epbwrite    = .true.    ! Writes the coarse matrix elements in the Bloch representation. 
  epbread     = .false.		! Default is false. 

  prtgkk      = .true.

  nbndsub     =  {input['ph_bands']}        ! Default is zero. 

  wannierize  = .true.
  num_iter    = 300
  iprint      = 2
  proj(1)     = 'random'

  iverbosity  = 1

  dvscf_dir   = '../1.2-ph/save'
  
  nk1         = {input['ph_kpts'][0]}
  nk2         = {input['ph_kpts'][1]}
  nk3         = {input['ph_kpts'][2]}

  nq1         = {input['ph_kpts'][0]}
  nq2         = {input['ph_kpts'][1]}
  nq3         = {input['ph_kpts'][2]}

  nqf1         = {input['ph_kpts'][0]}
  nqf2         = {input['ph_kpts'][1]}
  nqf3         = {input['ph_kpts'][2]}
  nkf1         = {input['ph_kpts'][0]}
  nkf2         = {input['ph_kpts'][1]}
  nkf3         = {input['ph_kpts'][2]}
 /
'''
    )

# endregion

# region: Jobscripts.
# jobscript_epw.run
n_pools = input['ph_kpts'][0]*input['ph_kpts'][1]*input['ph_kpts'][2]
fname = 'jobscript_epw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}={int(n_pools)}
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} epw.x -nk {n_pools} < epw.in &> epw.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
mpirun -n {n_pools} epw.x -nk {n_pools} < epw.in &> epw.out
'''
        )  

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.3-epw/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 1.5-bands
dir_name = '1.5-bands'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['bands_pwbands'],
                # 'nosym': True            
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],  
        kpts=input['pwscf_kpts']         # Need to overwrite k-points. 
    )
# pwbands.in: copy k-points to file. 
pattern = ... + pp.Literal('K_POINTS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)
parsed_text = pattern.parse_file(fname).as_list()
new_text = ''.join([
    parsed_text[0],
    input['bands_plot_kpath_string'],
    parsed_text[3],
    parsed_text[4]
])
with open(fname, 'w') as f: f.write(new_text)

# bands.in
fname = 'bands.in'
with open(fname, 'w') as f:
    f.write(
f'''&BANDS
outdir='./tmp'
prefix='{input['struct_formula']}'
filband='{input['struct_formula']}_dft_bandstructure.dat'
/
'''
    )

# pw2bgw.in
fname = 'pw2bgw.in'
with open(fname, 'w') as f:
    f.write(
f'''&input_pw2bgw
   outdir = './tmp/'
   prefix = '{input["struct_formula"]}'
   real_or_complex = 2
   wfng_flag = .true.
   wfng_file = 'WFN'
   wfng_kgrid = .true.
   wfng_nk1 = 0
   wfng_nk2 = 0
   wfng_nk3 = 0
   wfng_dk1 = 0.0
   wfng_dk2 = 0.0
   wfng_dk3 = 0.0
   vscg_flag = .true.
   vscg_file = 'VSC'
   vkbg_flag = .true.
   vkbg_file = 'VKB'
/
'''
    )

# parabands.inp
fname = 'parabands.inp'
with open(fname, 'w') as f:
    f.write(
f'''
input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input['bands_parabands']}
'''
    )

# endregion

# region: Jobscripts.
# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in &> pwbands.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in &> pwbands.out
'''
        ) 


# jobscript_bands.run
fname = 'jobscript_bands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} bands.x < bands.in &> bands.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
bands.x < bands.in &> bands.out
'''
        ) 

# jobscript_pw2bgw.run
fname = 'jobscript_pw2bgw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
module load espresso 
{input['slurm_mpi_type']} pw2bgw.x < pw2bgw.in &> pw2bgw.out
cp ./tmp/WFN ./WFN
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands, run below. Else, comment out.
#wfn2hdf.x BIN WFN WFN.h5
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/WFN ./WFN
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands, run below. Else, comment out.
#wfn2hdf.x BIN WFN WFN.h5
'''
        )

# jobscript_parabands.run
fname = 'jobscript_parabands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} parabands.cplx.x &> parabands.out
mv WFN_generated WFN.h5
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
parabands.cplx.x &> parabands.out
mv WFN_generated WFN.h5
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.5-bands/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')  
# endregion

# region: 1.6-dos
dir_name = '1.6-dos'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['dos_bands']              
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],  # Need to write k-points. 
        kpts=input['dos_kpts']
    )

# dos.in
fname = 'dos.in'
with open(fname, 'w') as f:
    f.write(
f'''&DOS
  outdir='./tmp/',
  prefix='{input['struct_formula']}',
  fildos='{input['struct_formula']}_dos.dat',
  emin={input['dos_e_min']},
  emax={input['dos_e_max']}
/
'''
    )
# endregion

# region: Jobscripts.
# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in &> pwbands.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in &> pwbands.out
'''
        ) 

# jobscript_dos.run 
fname = 'jobscript_dos.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} dos.x -pd .true. < dos.in &> dos.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
dos.x -pd .true. < dos.in &> dos.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.6-dos/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 1.7-pdos
dir_name = '1.7-pdos'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pdos.in
fname = 'pdos.in'
with open(fname, 'w') as f:
    f.write(
f'''&PROJWFC
  outdir= './tmp/',
  prefix= '{input['struct_formula']}',
  Emin = {input['dos_e_min']},
  Emax = {input['dos_e_max']},
  filpdos= '{input['struct_formula']}_pdos.dat'
/
'''
    )

# endregion

# region: Jobscripts.
# jobscript_pdos.run
fname = 'jobscript_pdos.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} projwfc.x -pd .true. < pdos.in &> pdos.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
projwfc.x -pd .true. < pdos.in &> pdos.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.7-pdos/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 1.8-pp
dir_name = '1.8-pp'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pp_charge_density.in
fname = 'pp_charge_density.in'
with open(fname, 'w') as f:
    f.write(
f'''&INPUTPP
  prefix= '{input['struct_formula']}',
  outdir= './tmp/',
!  filplot='{input['struct_formula']}_charge_density.dat',
  plotnum=0  
/

&PLOT
  iflag=3,  ! 3D plot.
  output_format=5, ! 3s .xsf
  filout='{input['struct_formula']}_charge_density.xsf'
/
'''
    )

# endregion

# region: Jobscripts.
# jobscript_pp_charge_density.run
fname = 'jobscript_pp_charge_density.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


./link_script.sh
module load impi
{input['slurm_mpi_type']} pp.x < pp_charge_density.in &> pp_charge_density.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pp.x < pp_charge_density.in &> pp_charge_density.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/1.8-pp/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 2.1-wfn
dir_name = '2.1-wfn'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['wfn_pw_bands'],
                # 'nosym': True             
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['wfn_kpts']
    )

# pw2bgw.in
fname = 'pw2bgw.in'
with open(fname, 'w') as f:
    f.write(
f'''&input_pw2bgw
   prefix = '{input['struct_formula']}'
   outdir = './tmp/'
   real_or_complex = 2
   wfng_flag = .true.
   wfng_file = 'WFN'
   wfng_kgrid = .true.
   wfng_nk1 = {input['wfn_kpts'][0]}
   wfng_nk2 = {input['wfn_kpts'][1]}
   wfng_nk3 = {input['wfn_kpts'][2]}
   wfng_dk1 = 0.0
   wfng_dk2 = 0.0
   wfng_dk3 = 0.0
   rhog_flag = .true.
   rhog_file = 'RHO'
   vxcg_flag = .true.
   vxcg_file = 'VXC'
   vxc_flag = .true.
   vxc_file = 'vxc.dat'
   vxc_diag_nmin = {input['wfn_vxc_diag_min']}
   vxc_diag_nmax = {input['wfn_vxc_diag_max']}
   vxc_offdiag_nmin = {input['wfn_vxc_offdiag_min']}
   vxc_offdiag_nmax = {input['wfn_vxc_offdiag_max']}
   vscg_flag = .true.
   vscg_file = 'VSC'
   vkbg_flag = .true.
   vkbg_file = 'VKB'
/
'''
    )

# parabands.inp
fname = 'parabands.inp'
with open(fname, 'w') as f:
    f.write(
f'''
input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input['wfn_parabands']}
'''
    )

# copy_kgrid.py
fname = 'copy_kgrid.py'
copy_kgrid_text = \
'''
import pyparsing as pp
import re

pattern = ... + pp.Literal('K_POINTS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('./pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: 
    subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: 
    f.write(new_text)
'''
with open(fname, 'w') as f: f.write(copy_kgrid_text)

# region: copy_degen_allowed_bands.py
fname = 'copy_degen_allowed_bands.py'
copy_degen_allowed_text = \
f'''
import pyparsing as pp
import re
import os
import numpy as np

# Debugging. 
# os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/2.1-wfn')


# Proposed max bands, val bands and conduction bands. 
eps_bands = {input['epsilon_bands']}
sigma_bands = {input['sigma_bands']}
sigma_min = {input['sigma_bands_min']}
sigma_max = {input['sigma_bands_max']}
coarse_val_bands = {input['coarse_val_bands']}
fine_val_bands = {input['fine_val_bands']}
coarse_cond_bands = {input['coarse_cond_bands']}
fine_cond_bands = {input['fine_cond_bands']}

# File names. 
epsilon = '../4-epsilon/epsilon.inp'
sigma = '../5.1-sigma/sigma.inp'
kernel = '../6-kernel/kernel.inp'
absorption = '../7-absorption/absorption.inp'


# Functions. 
def replace_eps(bnd):
    with open(epsilon, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_bands\s+\d+', f'number_bands {{bnd}}', txt)
    with open(epsilon, 'w') as f: txt = f.write(new_txt)

def replace_sigma(bnd, min, max):
    
    # Sigma bands. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_bands\s+\d+', f'number_bands {{bnd}}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

    # min. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'band_index_min\s+\d+', f'band_index_min {{min}}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

    # max. 
    with open(sigma, 'r') as f: txt = f.read()
    new_txt = re.sub(r'band_index_max\s+\d+', f'band_index_max {{max}}', txt)
    with open(sigma, 'w') as f: txt = f.write(new_txt)

def replace_kernel(vbc, cbc):
    # vbc
    with open(kernel, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands\s+\d+', f'number_val_bands {{vbc}}', txt)
    with open(kernel, 'w') as f: txt = f.write(new_txt)

    # cbc
    with open(kernel, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands\s+\d+', f'number_cond_bands {{cbc}}', txt)
    with open(kernel, 'w') as f: txt = f.write(new_txt)

def replace_absorption(vbf, vbc, cbf, cbc):
    # vbc
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands_coarse\s+\d+', f'number_val_bands_coarse {{vbc}}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # vbf
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_val_bands_fine\s+\d+', f'number_val_bands_fine {{vbf}}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # cbc
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands_coarse\s+\d+', f'number_cond_bands_coarse {{cbc}}', txt)
    with open(absorption, 'w') as f: txt = f.write(new_txt)

    # cbf
    with open(absorption, 'r') as f: txt = f.read()
    new_txt = re.sub(r'number_cond_bands_fine\s+\d+', f'number_cond_bands_fine {{cbf}}', txt)
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

print(f'Outer array length: {{len(outer_array)}}')

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
'''
with open(fname, 'w') as f: f.write(copy_degen_allowed_text)
# endregion
# endregion

# region: Jobscripts.
# jobscript_kgrid.run
fname = 'jobscript_kgrid.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
{input['slurm_mpi_type']} kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
'''
        )

# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in &> pwbands.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
mpirun -n 8 pw.x < pwbands.in &> pwbands.out
'''
        ) 

# jobscript_pw2bgw.run
fname = 'jobscript_pw2bgw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw2bgw.x -pd .true. < pw2bgw.in &> pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
# wfn2hdf.x BIN WFN WFN.h5
# Below if no parabands. Will update the degen allowed bands in files. 
#python copy_degen_allowed_bands.py    
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# Below if no parabands. Comment out if using parabands to extend bands. 
#wfn2hdf.x BIN WFN WFN.h5
# Below if no parabands. Will update the degen allowed bands in files. 
#python copy_degen_allowed_bands.py          
'''
        ) 

# jobscript_parabands.run
fname = 'jobscript_parabands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} parabands.cplx.x &> parabands.out
mv WFN_generated WFN.h5
hdf2wfn.x BIN WFN.h5 WFN
python copy_degen_allowed_bands.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
parabands.cplx.x &> parabands.out
mv WFN_generated WFN.h5
hdf2wfn.x BIN WFN.h5 WFN
python copy_degen_allowed_bands.py
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/2.1-wfn/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 2.2-wfnq
dir_name = '2.2-wfnq'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['wfn_pw_bands']             
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['wfn_kpts']
    )

# pw2bgw.in
fname = 'pw2bgw.in'
with open(fname, 'w') as f:
    f.write(
f'''&input_pw2bgw
   prefix = '{input['struct_formula']}'
   outdir = './tmp/'
   real_or_complex = 2
   wfng_flag = .true.
   wfng_file = 'WFN'
   wfng_kgrid = .true.
   wfng_nk1 = {input['wfn_kpts'][0]}
   wfng_nk2 = {input['wfn_kpts'][1]}
   wfng_nk3 = {input['wfn_kpts'][2]}
   wfng_dk1 = {input['wfn_qshift'][0]}
   wfng_dk2 = {input['wfn_qshift'][1]}
   wfng_dk3 = {input['wfn_qshift'][2]}
   rhog_flag = .true.
   rhog_file = 'RHO'
   vxcg_flag = .true.
   vxcg_file = 'VXC'
   vxc_flag = .true.
   vxc_file = 'vxc.dat'
   vxc_diag_nmin = {input['wfn_vxc_diag_min']}
   vxc_diag_nmax = {input['wfn_vxc_diag_min']}
   vxc_offdiag_nmin = {input['wfn_vxc_offdiag_min']}
   vxc_offdiag_nmax = {input['wfn_vxc_offdiag_max']}
   vscg_flag = .true.
   vscg_file = 'VSC'
   vkbg_flag = .true.
   vkbg_file = 'VKB'
/
'''
    )

# parabands.inp
fname = 'parabands.inp'
with open(fname, 'w') as f:
    f.write(
f'''
input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input['wfnq_parabands']}
'''
    )

# copy_kgrid.py
fname = 'copy_kgrid.py'
copy_kgrid_text = \
'''
import pyparsing as pp
import re

pattern = ... + pp.Literal('K_POINTS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('./pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: 
    subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: 
    f.write(new_text)
'''
with open(fname, 'w') as f: f.write(copy_kgrid_text)

# endregion

# region: Jobscripts.
# jobscript_kgrid.run
fname = 'jobscript_kgrid.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift {input['wfn_qshift'][0]} {input['wfn_qshift'][1]} {input['wfn_qshift'][2]} ./tmp/{input['struct_formula']}.xml kgrid.inp
{input['slurm_mpi_type']} kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
data-file2kgrid.py --kgrid {input['wfn_kpts'][0]} {input['wfn_kpts'][1]} {input['wfn_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift {input['wfn_qshift'][0]} {input['wfn_qshift'][1]} {input['wfn_qshift'][2]} ./tmp/{input['struct_formula']}.xml kgrid.inp
kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
'''
        )

# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in > pwbands.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in > pwbands.out
'''
        ) 

# jobscript_pw2bgw.run
fname = 'jobscript_pw2bgw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw2bgw.x -pd .true. < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFNq.h5
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFNq.h5
'''
        ) 

# jobscript_parabands.run
fname = 'jobscript_parabands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} parabands.cplx.x &> parabands.out
mv WFN_generated WFNq.h5
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
parabands.cplx.x &> parabands.out
mv WFN_generated WFNq.h5
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/2.2-wfnq/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 3.1-wfn_fi
dir_name = '3.1-wfn_fi'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['wfn_pw_bands']             
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['wfn_fi_kpts']
    )

# pw2bgw.in
fname = 'pw2bgw.in'
with open(fname, 'w') as f:
    f.write(
f'''&input_pw2bgw
   prefix = '{input['struct_formula']}'
   outdir = './tmp/'
   real_or_complex = 2
   wfng_flag = .true.
   wfng_file = 'WFN'
   wfng_kgrid = .true.
   wfng_nk1 = {input['wfn_fi_kpts'][0]}
   wfng_nk2 = {input['wfn_fi_kpts'][1]}
   wfng_nk3 = {input['wfn_fi_kpts'][2]}
   wfng_dk1 = 0.0
   wfng_dk2 = 0.0
   wfng_dk3 = 0.0
   rhog_flag = .true.
   rhog_file = 'RHO'
   vxcg_flag = .true.
   vxcg_file = 'VXC'
   vxc_flag = .true.
   vxc_file = 'vxc.dat'
   vxc_diag_nmin = {input['wfn_vxc_diag_min']}
   vxc_diag_nmax = {input['wfn_vxc_diag_min']}
   vxc_offdiag_nmin = {input['wfn_vxc_offdiag_min']}
   vxc_offdiag_nmax = {input['wfn_vxc_offdiag_max']}
   vscg_flag = .true.
   vscg_file = 'VSC'
   vkbg_flag = .true.
   vkbg_file = 'VKB'
/
'''
    )

# parabands.inp
fname = 'parabands.inp'
with open(fname, 'w') as f:
    f.write(
f'''
input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input['wfnq_parabands']}
'''
    )

# copy_kgrid.py
fname = 'copy_kgrid.py'
copy_kgrid_text = \
'''
import pyparsing as pp
import re

pattern = ... + pp.Literal('K_POINTS') + ... + pp.ZeroOrMore(pp.srange('A-Z')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: f.write(new_text)
'''
with open(fname, 'w') as f: f.write(copy_kgrid_text)

# endregion

# region: Jobscripts.
# jobscript_kgrid.run
fname = 'jobscript_kgrid.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
module load espresso 
data-file2kgrid.py --kgrid {input['wfn_fi_kpts'][0]} {input['wfn_fi_kpts'][1]} {input['wfn_fi_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
{input['slurm_mpi_type']} kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
data-file2kgrid.py --kgrid {input['wfn_fi_kpts'][0]} {input['wfn_fi_kpts'][1]} {input['wfn_fi_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift 0.0 0.0 0.0 ./tmp/{input['struct_formula']}.xml kgrid.inp
kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
'''
        )

# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in > pwbands.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in > pwbands.out
'''
        ) 

# jobscript_pw2bgw.run
fname = 'jobscript_pw2bgw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFN_fi.h5
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFN_fi.h5
'''
        ) 

# jobscript_parabands.run
fname = 'jobscript_parabands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} parabands.cplx.x &> parabands.out
mv WFN_generated WFN_fi.h5
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
parabands.cplx.x &> parabands.out
mv WFN_generated WFN_fi.h5
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/3.1-wfn_fi/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 3.2-wfnq_fi
dir_name = '3.2-wfnq_fi'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

ln -sf ../1.1-scf/tmp ./
ln -sf ../1.1-scf/pseudo ./
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# pwbands.in
fname = 'pwbands.in'
with open(fname, 'w') as f:
    espresso.write_espresso_in(
        f,
        struct,
        input_data={
            'control': {
                'outdir': './tmp',
                'prefix': f'{input["struct_formula"]}',
                'pseudo_dir': './pseudo',
                'calculation': 'bands' 
            },
            'system': {
                'ecutwfc': input['pwscf_wfc_cutoff'],
                'nbnd': input['wfn_pw_bands']             
            },
            'electrons': {
            }
        },
        pseudopotentials=input['pw_pseudopotentials'],
        kpts=input['wfn_fi_kpts']
    )

# pw2bgw.in
fname = 'pw2bgw.in'
with open(fname, 'w') as f:
    f.write(
f'''&input_pw2bgw
   prefix = '{input['struct_formula']}'
   outdir = './tmp/'
   real_or_complex = 2
   wfng_flag = .true.
   wfng_file = 'WFN'
   wfng_kgrid = .true.
   wfng_nk1 = {input['wfn_fi_kpts'][0]}
   wfng_nk2 = {input['wfn_fi_kpts'][1]}
   wfng_nk3 = {input['wfn_fi_kpts'][2]}
   wfng_dk1 = {input['wfn_qshift'][0]}
   wfng_dk2 = {input['wfn_qshift'][1]}
   wfng_dk3 = {input['wfn_qshift'][2]}
   rhog_flag = .true.
   rhog_file = 'RHO'
   vxcg_flag = .true.
   vxcg_file = 'VXC'
   vxc_flag = .true.
   vxc_file = 'vxc.dat'
   vxc_diag_nmin = {input['wfn_vxc_diag_min']}
   vxc_diag_nmax = {input['wfn_vxc_diag_min']}
   vxc_offdiag_nmin = {input['wfn_vxc_offdiag_min']}
   vxc_offdiag_nmax = {input['wfn_vxc_offdiag_max']}
   vscg_flag = .true.
   vscg_file = 'VSC'
   vkbg_flag = .true.
   vkbg_file = 'VKB'
/
'''
    )

# parabands.inp
fname = 'parabands.inp'
with open(fname, 'w') as f:
    f.write(
f'''
input_wfn_file WFN
output_wfn_file WFN_generated
vkb_file VKB
vsc_file VSC

number_bands {input['wfnq_parabands']}
'''
    )

# copy_kgrid.py
fname = 'copy_kgrid.py'
copy_kgrid_text = \
'''
import pyparsing as pp
import re

pattern = ... + pp.Literal('K_POINTS') + ... + pp.ZeroOrMore(pp.srange('A-Z')) + pp.Regex(r'.*', re.DOTALL)

parsed_text = pattern.parse_file('pwbands.in').as_list()

subst_text = ''
with open('kgrid.out', 'r') as f: subst_text = f.read()

new_text = ''.join([
    parsed_text[0],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('pwbands.in', 'w') as f: f.write(new_text)
'''
with open(fname, 'w') as f: f.write(copy_kgrid_text)

# endregion

# region: Jobscripts.
# jobscript_kgrid.run
fname = 'jobscript_kgrid.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
module load espresso 
data-file2kgrid.py --kgrid {input['wfn_fi_kpts'][0]} {input['wfn_fi_kpts'][1]} {input['wfn_fi_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift {input['wfn_qshift'][0]} {input['wfn_qshift'][1]} {input['wfn_qshift'][2]} ./tmp/{input['struct_formula']}.xml kgrid.inp
{input['slurm_mpi_type']} kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
data-file2kgrid.py --kgrid {input['wfn_fi_kpts'][0]} {input['wfn_fi_kpts'][1]} {input['wfn_fi_kpts'][2]} --kshift 0.0 0.0 0.0 --qshift {input['wfn_qshift'][0]} {input['wfn_qshift'][1]} {input['wfn_qshift'][2]} ./tmp/{input['struct_formula']}.xml kgrid.inp
kgrid.x kgrid.inp kgrid.out kgrid.log
python copy_kgrid.py
'''
        )

# jobscript_pwbands.run
fname = 'jobscript_pwbands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw.x < pwbands.in > pwbands.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
f'''#!/bin/bash

./link_script.sh
pw.x < pwbands.in > pwbands.out
'''
        ) 

# jobscript_pw2bgw.run
fname = 'jobscript_pw2bgw.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type']  == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFNq_fi.h5
'''
        )
    elif input['sched_type']  == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
pw2bgw.x < pw2bgw.in > pw2bgw.out
cp ./tmp/RHO ./RHO
cp ./tmp/WFN ./WFN
cp ./tmp/vxc.dat ./vxc.dat
cp ./tmp/VXC ./VXC
cp ./tmp/VSC ./VSC
cp ./tmp/VKB ./VKB

# If no parabands. Comment out if using parabands to extend bands. 
wfn2hdf.x BIN WFN WFNq_fi.h5
'''
        ) 

# jobscript_parabands.run
fname = 'jobscript_parabands.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
module load impi
{input['slurm_mpi_type']} parabands.cplx.x &> parabands.out
mv WFN_generated WFNq_fi.h5
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
parabands.cplx.x &> parabands.out
mv WFN_generated WFNq_fi.h5
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/3.2-wfnq_fi/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 4-epsilon
dir_name = '4-epsilon'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash -l


# Copy for epsilon folder
ln -sf ../2.1-wfn/WFN.h5 ./WFN.h5
ln -sf ../2.2-wfnq/WFNq.h5 ./WFNq.h5
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# epsilon.inp
fname = 'epsilon.inp'
with open(fname, 'w') as f:
    f.write(
f'''
epsilon_cutoff {input['epsilon_cutoff']}
number_bands {input['epsilon_bands']}
use_wfn_hdf5

begin qpoints
end
'''
    )

# copy_qgrid.py
fname = 'copy_qgrid.py'
copy_qgrid_text = \
f'''
import pyparsing as pp
import re
import numpy as np

# Create the pattern. 
pattern = ... + pp.Literal('begin qpoints\\n') + ... + pp.Literal('end') + pp.Regex(r'.*', re.DOTALL)

# Parse epsilon.inp.
parsed_text = pattern.parse_file('epsilon.inp').as_list()

# Create text to substitute. 
subst_text = ''
data = np.loadtxt('../2.1-wfn/kgrid.out', skiprows=2)
shape = data.shape
if len(shape) == 1:
    data = data.reshape(1, shape[0])
max_rows = data.shape[0]
new_data = np.zeros((max_rows, 5), dtype='f8')
new_data[:, 0:4] = data
new_data[0, 0] = {input['wfn_qshift'][0]}  # For the qshifted point.
new_data[0, 1] = {input['wfn_qshift'][1]}   # For the qshifted point.
new_data[0, 2] = {input['wfn_qshift'][2]}   # For the qshifted point.
new_data[0, 3] = 1.0  # For the qshifted point.
new_data[0, 4] = 1    # For the qshifted point.

with open('epsilon_qpoints.csv', 'w') as f:
    for row in new_data:
        f.write(f'{{row[0]:10.6f}} {{row[1]:10.6f}} {{row[2]:10.6f}} {{row[3]:10.6f}} {{int(row[4])}}\\n')

with open('epsilon_qpoints.csv', 'r') as f: subst_text = f.read()

# Create the new text. 
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('epsilon.inp', 'w') as f: f.write(new_text)
'''

with open(fname, 'w') as f: f.write(copy_qgrid_text)

# endregion

# region: Jobscripts.
# jobscript_epsilon.run
fname = 'jobscript_epsilon.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


./link_script.sh
python copy_qgrid.py
{input['slurm_mpi_type']} epsilon.cplx.x &> epsilon.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
python copy_qgrid.py
mpirun -n 8 epsilon.cplx.x &> epsilon.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/4-epsilon/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 5.1-sigma
dir_name = '5.1-sigma'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

# Copy for sigma folder

# From wfn folder. 
ln -sf ../2.1-wfn/vxc.dat ./vxc.dat
ln -sf ../2.1-wfn/RHO ./RHO
ln -sf ../2.1-wfn/WFN.h5 ./WFN_inner.h5


# From epsilon folder. 
ln -sf ../4-epsilon/eps0mat.h5 ./eps0mat.h5
ln -sf ../4-epsilon/epsmat.h5 ./epsmat.h5
'''
    )
os.system(f'chmod u+x {fname}')

# endregion

# region: Input files
# sigma.inp
fname = 'sigma.inp'
with open(fname, 'w') as f:
    f.write(
f'''
number_bands {input['sigma_bands']}

band_index_min {input['sigma_bands_min']}
band_index_max {input['sigma_bands_max']}

screening_semiconductor

use_wfn_hdf5

begin kpoints
end
'''
    )

# copy_kpts_sigma.py
fname = 'copy_kpts_sigma.py'
copy_kgrid_text = \
f'''
import pyparsing as pp
import re
import numpy as np

# Create the pattern. 
pattern = ... + pp.Literal('begin kpoints') + ... + pp.Literal('end') + pp.Regex(r'.*', re.DOTALL)

# Parse epsilon.inp.
parsed_text = pattern.parse_file('sigma.inp').as_list()

# Create text to substitute. 
subst_text = ''
data = np.loadtxt('../2.1-wfn/kgrid.out', skiprows=2)
shape = data.shape
if len(shape) == 1:
    data = data.reshape(1, shape[0])
np.savetxt('sigma_kpoints.csv', data, delimiter=' ', fmt='%10.6f')
with open('sigma_kpoints.csv', 'r') as f: subst_text = f.read()

# Create the new text. 
new_text = ''.join([
    parsed_text[0],
    parsed_text[1],
    '\\n',
    subst_text,
    parsed_text[3],
    parsed_text[4]
])

with open('sigma.inp', 'w') as f: f.write(new_text)
'''

with open(fname, 'w') as f: f.write(copy_kgrid_text)

# endregion

# region: Jobscripts.
# jobscript_sigma.run
jobscript_name = 'jobscript_sigma.run' 
touch(jobscript_name)
with open(jobscript_name, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
python copy_kpts_sigma.py
{input['slurm_mpi_type']} sigma.cplx.x &> sigma.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
python copy_kpts_sigma.py
mpirun -n 8 sigma.cplx.x &> sigma.out
'''
        )

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/5.1-sigma/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 5.2-inteqp
dir_name = '5.2-inteqp'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash

# From wfn folder. Coarse wavefunctions. 
ln -sf ../2.1-wfn/WFN.h5 ./WFN_co.h5

# From bands folder. Fine wavefunctions. 
ln -sf ../1.5-bands/WFN.h5 ./WFN_fi.h5

# From sigma folder. Quasiparticle energy corrections. 
ln -sf ../5.1-sigma/eqp1.dat ./eqp_co.dat 
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# inteqp.inp
fname = 'inteqp.inp'
with open(fname, 'w') as f:
    f.write(
f'''
number_val_bands_coarse {input['coarse_val_bands']}
number_val_bands_fine {input['fine_val_bands']}
number_cond_bands_coarse {input['coarse_cond_bands']}
number_cond_bands_fine 6 {input['fine_cond_bands']}

use_wfn_hdf5

use_symmetries_coarse_grid
no_symmetries_fine_grid
'''
    )


# endregion

# region: Jobscripts.
# jobscript_inteqp.run
fname = 'jobscript_inteqp.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=10
{input['slurm_ntasks_flag']}=560
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
{input['slurm_mpi_type']} inteqp.cplx.x &> inteqp.out
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
inteqp.cplx.x &> inteqp.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/5.2-inteqp/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 6-kernel ..
dir_name = '6-kernel'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 

# Copy for kernel folder

# From wfn folder. 
ln -sf ../2.1-wfn/WFN.h5 ./WFN_co.h5


# From epsilon folder. 
ln -sf ../4-epsilon/eps0mat.h5 ./eps0mat.h5
ln -sf ../4-epsilon/epsmat.h5 ./epsmat.h5
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# kernel.inp
fname = 'kernel.inp'
with open(fname, 'w') as f:
    f.write(
f'''
number_val_bands {input['coarse_val_bands']}
number_cond_bands {input['coarse_cond_bands']}

use_wfn_hdf5

screening_semiconductor
'''
    )

# endregion

# region: Jobscripts.
# jobscript_kernel.run
fname = 'jobscript_kernel.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
{input['slurm_mpi_type']} kernel.cplx.x &> kernel.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
mpirun -n 4 kernel.cplx.x &> kernel.out
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/6-kernel/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 7-absorption
dir_name = '7-absorption'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
'''#!/bin/bash 


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
'''
    )
os.system(f'chmod u+x {fname}')
# endregion

# region: Input files.
# absorption.inp
fname = 'absorption.inp'
with open(fname, 'w') as f:
    f.write(
f'''
diagonalization

use_wfn_hdf5

number_val_bands_coarse {input['coarse_val_bands']}
number_val_bands_fine {input['fine_val_bands']}
number_cond_bands_coarse {input['coarse_cond_bands']}
number_cond_bands_fine {input['fine_cond_bands']}

screening_semiconductor

use_velocity

gaussian_broadening
energy_resolution 0.15

dump_bse_hamiltonian

eqp_co_corrections
write_eigenvectors {input['max_eigenvectors']}
'''
    )

# plotxct.inp
fname = 'plotxct.inp'
with open(fname, 'w') as f:
    f.write(
f'''
hole_position {input['wfn_fi_kpts'][0]/2} {input['wfn_fi_kpts'][1]/2} {input['wfn_fi_kpts'][2]/2} 
plot_spin 1
plot_state {input['plotxct_state']}
q_shift {input['wfn_qshift'][0]} {input['wfn_qshift'][1]} {input['wfn_qshift'][2]}
supercell_size {input['wfn_fi_kpts'][0]} {input['wfn_fi_kpts'][1]} {input['wfn_fi_kpts'][2]} 
verbosity 2
use_wfn_hdf5
'''
    )

# endregion

# region: Jobscripts.
# jobscript_absorption.run
fname = 'jobscript_absorption.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
{input['slurm_mpi_type']} absorption.cplx.x &> absorption.out
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
mpirun -n 4 absorption.cplx.x &> absorption.out
'''
        ) 


# jobscript_plotxct.run
fname = 'jobscript_plotxct.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=20
{input['slurm_ntasks_flag']}=1120
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}


./link_script.sh
{input['slurm_mpi_type']} plotxct.cplx.x &> plotxct.out
volume.py ../1.1-scf/pwscf.in espresso *.a3Dr a3dr xct_1.xsf xsf false abs2 true
rm *a3Dr
rm slurm*
'''
        )

    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
plotxct.cplx.x &> plotxct.out
volume.py ../1.1-scf/pwscf.in espresso *.a3Dr a3dr xct_1.xsf xsf false abs2 true
rm *a3Dr
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/7-absorption/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')
# endregion

# region: 8-esf
dir_name = '8-esf'
mkdir(dir_name)
os.chdir(dir_name)

# region: Link script.
fname = 'link_script.sh'
touch(fname)
with open(fname, 'w') as f:
    f.write(
f'''#!/bin/bash

# dft_force. 
ln -sf ../1.1-scf/pwscf.in ./pwscf.in
ln -sf ../1.1-scf/pwscf.out ./pwscf.out
ln -sf ../2.1-wfn/WFN.h5 ./WFN.h5

# kgrid. 
ln -sf ../2.1-wfn/kgrid.out ./kgrid.out
#This overrides the above step and creates the full kgrid order using kmesh.pl. 
kmesh.pl {input['ph_kpts'][0]} {input['ph_kpts'][1]} {input['ph_kpts'][2]} > kgrid.out          

# eqp (eqp).
ln -sf ../7-absorption/eqp.dat ./eqp.dat

# elph (g).
ln -sf ../1.3-epw/epw.in ./epw.in
file_idx=0
file_name=""
for file in $(find ../1.1-scf/tmp -name {input['struct_formula']}_elph_*)
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
'''
    )
os.system(f'chmod u+x {fname}') 

# endregion

# region: Input files.
# endregion

# region: Jobscripts.
# jobscript_esf.run
fname = 'jobscript_esf.run' 
touch(fname)
with open(fname, 'w') as f:
    if input['sched_type'] == 'slurm':
        f.write(
f'''#!/bin/bash -l
{input['slurm_account_flag']}
{input['slurm_partition_flag']}=development
{input['slurm_job_flag']}
{input['slurm_nodes_flag']}=1
{input['slurm_ntasks_flag']}=56
{input['slurm_time_flag']}=01:30:00
{input['slurm_mail_user_flag']}
{input['slurm_mail_type_flag']}

./link_script.sh
python esf.py &> esf.log
rm slurm*
'''
        )
    elif input['sched_type'] == 'WSL':
        f.write(
'''#!/bin/bash

./link_script.sh
python esf.py &> esf.log
'''
        ) 

# endregion

# region: Utilities.
# Copy relevant directories.
os.system('cp -r ../utilities/8-esf/* ./')
# endregion

# Go back to parent directory. 
os.chdir('../')   
# endregion