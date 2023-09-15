import os 
import time
from petsc4py import PETSc 
import numpy as np
import h5py 
import re 
import json 
import subprocess


# Debugging. 
# os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/8-esf')

# region: Global variables. 
sched_type = 'slurm'   # WSL or slurm. 
is_time_logging = True   # Want to record how much time each one takes? 
# endregion

# region: Global classes. 
class Logger:

    def __init__(self):
        self.start = 0.0
        self.stop = 0.0

    def time_start(self):
        self.start = time.time()

    def time_stop(self):
        self.stop = time.time()

    def time_elapsed(self):
        return self.stop - self.start
local_logger = Logger()
global_logger = Logger()
# endregion

# region: Global functions. 
def print_flush(msg):
    print(msg, flush=True)

def run_cmd(cmd):
    global sched_type
    global local_logger
    global global_logger

def print_flush(msg):
    print(msg, flush=True)

def run_cmd(cmd):
    global sched_type
    global local_logger
    global global_logger

    # Starting a command. 
    print_flush(f'Starting : {cmd}')
    ps_result = -1

    # Run the command. 
    if is_time_logging: local_logger.time_start()

    if sched_type == 'WSL':
        ps_result = subprocess.run(f'{cmd}', shell=True)
    elif sched_type == 'slurm':
        ps_result = subprocess.run(f'sbatch --wait {cmd}', shell=True)

    if is_time_logging: local_logger.time_stop()

    # See if it finished properly. 
    if ps_result.returncode == 0:  # No error. 
        # Done with a command. 
        if is_time_logging:
            print_flush(f'Done with: {cmd}. Took {local_logger.time_elapsed()} seconds.\n\n')
        else:
            print_flush(f'Done with: {cmd}.\n\n')
    else:
        # Error. 
        if is_time_logging:
            print_flush(f'Error finishing: {cmd}. Exited with code {ps_result.returncode}. Took time {local_logger.time_elapsed()}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        else:
            print_flush(f'Error finishing: {cmd}. Exited with code {ps_result.returncode}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        os._exit(ps_result.returncode)     # Exit with the error code. 
# endregion

# region: Start global timer.
global_logger.time_start()
# endregion


# region: ESD variables. 
# TAO solver settings. 
tao = PETSc.TAO().create()
tao_solver_type = PETSc.TAO.Type.LMVM

max_tol = 9e-3  # Force/mass where force is in eV/A and mass is in amu. 
max_steps = 4
my_eps = 1e-4
iter_offset = 0    # This is the variable to continue from an interrupted iteration if needed. 

# Get number of atoms in the system. 
with open('../input.json', 'r') as f: input = json.loads(f.read())
initial_positions = input['positions']
nat = len(initial_positions)
flattened_positions = []
for loc in range(nat):
    for coord in range(3): 
        flattened_positions.append(initial_positions[loc][coord])
print_flush(f'Number of atoms read from input.json: {nat}')
print_flush(f'Flattened positions: {flattened_positions}')

# excited state varibles. 
xctd_pos = np.zeros((max_steps, nat*3), dtype='f8')
xctd_energy = np.zeros((max_steps,), dtype='f8')        # Total energy for the excited state, as a sum of dft and gw-bse. 
gwbse_energy = np.zeros((max_steps,), dtype='f8')       # Only the gwbse energy. 
xctd_force = np.zeros((max_steps, nat*3), dtype='f8')
xctd_amass = np.zeros((max_steps, nat), dtype='f8')

# Query variables. 
query_pos = np.zeros((nat*3,), dtype='f8') # Size will be 3*nat. 
query_e = 0.0 # Size will be 1.
query_gwbse_e = 0.0 # Size will be 1.
query_force = np.zeros((nat*3,), dtype='f8') # Size will be 3*nat. 
query_amass = np.zeros((nat,), dtype='f8') # Size will be nat. 



# Create PETSC force over mass vector.  
fm = PETSc.Vec().create()
fm.setType(PETSc.Vec.Type.SEQ)
fm.setSizes(3*nat)
fm.setUp()
fm.setValues(range(3*nat), np.zeros((3*nat,), dtype='f8'))
fm.assemble()
print_flush(f'Initial fm: {fm[:]}')

# Create PETSC position vector. 
pos = PETSc.Vec().create()
pos.setType(PETSc.Vec.Type.SEQ)
pos.setSizes(3*nat)
pos.setUp()
pos.setValues(range(3*nat), np.array(flattened_positions, dtype='f8'))
pos.assemble()
print_flush(f'Initial pos: {pos[:]}')
# endregion. 


# region: ESD functions.

def is_within_my_eps(x, y):
    # Calc max of diff. 
    diff_max = np.max(np.abs(x - y))

    # Check if it is within error bound. 
    if (diff_max > my_eps):
        return False
    else:
        return True


def set_read_json_flag(flag):
    print_flush(f'Calling set_read_json_flag with {flag}')

    def repl(match):
        return f'read_input_from_file = {flag}'

    flag_pattern = r'read_input_from_file\s+=\s+(?P<flag>True|False)' 

    with open('../create_workflow.py', 'r') as f: text = f.read()
    new_text = re.sub(flag_pattern, repl, text)
    with open('../create_workflow.py', 'w') as f: f.write(new_text)

def update_workflow_pos(pos):
    print_flush(f'Calling update_workflow_pos with {pos[:]}')
    
    # Reshape position. 
    pos_np = pos[:]
    nat = int(pos_np.shape[0]/3)
    pos_np = pos_np.reshape(nat, 3)
    pos_list = [row.tolist() for row in pos_np]    # Creates a list out of it. 
    print_flush(f'Position list is {pos_list}')

    # Update input.json.
    with open('../input.json', 'r') as f: input_txt = f.read() 
    input = json.loads(input_txt)                                                # Reads the input as json object. 
    input['positions'] = pos_list
    with open('../input.json', 'w')as f: f.write(json.dumps(input))

def remove_prev_step():
    print_flush('Removing directories.')
    inodes = os.listdir('./')

    for inode in inodes:

        # Remove only directory names that start with a number.
        if os.path.isdir(inode) and inode[0].isnumeric() and int(inode[0]) < 8: # Anything less that esf folder.   
            os.system(f'rm -r {inode}')
            print_flush(f'Removed directory {inode}')

def get_e_and_esf(pos):
    '''
    Function forces getting new energy and force after updating files with new position. 
    '''
    
    print_flush(f'Calling get_e_and_esf with {pos[:]}')

    # Write position to file. 
    update_workflow_pos(pos)

    # Do everything and run one iteration/step. 
    os.chdir('../')
    remove_prev_step()
    run_cmd('./jobscript_create_workflow.run')  # Create files if needed. 
    run_cmd('./jobscript_run_workflow.run')     # Run iteration. 
    os.chdir('./8-esf') # Get back to esf directory.


    # Read the values. 
    get_esf_info()

def return_e(pos):
    print_flush(f'Calling return_e with {pos[:]}')
    print_flush(f'Saved query_pos is :{query_pos}')
    
    # Run calculation if necessary. Then return exciton energy. 
    if not is_within_my_eps(query_pos, pos[:]):
        get_e_and_esf(pos)

    return query_e

def return_esf(pos):
    print_flush(f'Calling return_esf with {pos[:]}')
    print_flush(f'Saved query_pos is :{query_pos}')

    # Run calculation if necessary. Then return ESF. 
    if not is_within_my_eps(query_pos, pos[:]):
        get_e_and_esf(pos)

    return query_force


# region: PETSC objective, gradient, and monitor functions. 
def tao_objective(tao, pos):
    # return energy. 
    global query_pos
    global query_e
    global query_force

    print_flush(f'Calling tao_objective with {pos[:]}')

    query_e = return_e(pos)

    print_flush(f'tao_objective returned {query_e}')
    # query_pos = pos[:]      # To override read pos. 

    return query_e 

def tao_gradient(tao, pos, fm):
    # return force. 
    global query_pos
    global query_e
    global query_force

    print_flush(f'Calling tao_gradient with pos: {pos[:]} and fm: {fm[:]}')

    fm[:] = -return_esf(pos)/np.repeat(query_amass, 3)
    fm.assemble()

    print_flush(f'tao_gradient returned {fm[:]}')
    # query_pos = pos[:]      # To override read pos. 

def tao_objective_and_gradient(tao, pos, fm):
    # return energy and force. 
    global query_pos
    global query_e
    global query_force

    print_flush(f'Calling tao_objective_and_gradient and gradient with {pos[:]}, and fm: {fm[:]}')

    query_e = return_e(pos)

    fm[:] = -return_esf(pos)/np.repeat(query_amass, 3)
    fm.assemble()

    print_flush(f'tao_objective_and_gradient energy: {query_e}, fm: {fm[:]}')

    return query_e 

def tao_monitor(tao):
    global xctd_pos
    global xctd_energy
    global gwbse_energy
    global xctd_force
    global xctd_amass

    print_flush('Calling tao monitor.')

    # Get current iteration data. 
    current_iter = tao.its + iter_offset
    current_position = pos[:]
    current_energy = tao.objective
    current_force = -fm[:]*np.repeat(query_amass, 3)


    print_flush(f'TAO:')
    print_flush(f'Iteration number: {tao.getIterationNumber()}')
    print_flush(f'Converged reason: {tao.getConvergedReason()}')
    # print_flush(f'Function value: {tao.getFunctionValue()}')
    # print_flush(f'Objective value: {tao.getObjectiveValue()}')
    # print_flush(f'Gradient: {tao.getGradient()[:]}')
    # print_flush(f'Gradient norm: {tao.getGradientNorm()}')
    # print_flush(f'Max func evaluations: {tao.getMaximumFunctionEvaluations()}')
    # print_flush(f'Max iterations: {tao.getMaximumIterations()}')
    # print_flush(f'Solution: {tao.getSolution()[:]}')
    # print_flush(f'Solution Norm: {tao.getSolutionNorm()}')
    # print_flush(f'Solution status: {tao.getSolutionStatus()}')
    # print_flush(f'Tolerances: {tao.getTolerances()}')

    # Update current iteration position, energy and force. 
    xctd_pos[current_iter, :] = current_position     # Position in A. 
    xctd_energy[current_iter] = current_energy       # Energy in eV. 
    gwbse_energy[current_iter] = query_gwbse_e       # Energy in eV. 
    xctd_force[current_iter, :] = current_force   # Get the force in eV/A.
    xctd_amass[current_iter, :] =  query_amass      # Get the mass in amu. 

    # Update hdf5 file too. 
    write_esd_to_h5()

# endregion


def get_esf_info():
    global query_pos
    global query_e
    global query_gwbse_e
    global query_force
    global query_amass
    global nat

    print_flush(f'Calling get_esf_info')

    with h5py.File('esf.h5', 'r') as esf:
        query_pos = esf['position'][:].flatten()
        query_e = float(esf['energy'][()])
        query_gwbse_e = float(esf['gwbse_energy'][()])
        query_force = esf['force'][:].flatten()
        query_amass = esf['atomic_masses'][:]

        nat = int(query_amass.shape[0])

        print_flush(f'Read query_pos: {query_pos}')
        print_flush(f'Read query_e: {query_e}')
        print_flush(f'Read query_gwbse_e: {query_gwbse_e}')
        print_flush(f'Read query_force: {query_force}')
        print_flush(f'Read query_amass: {query_amass}')
        print_flush(f'Read nat: {nat}')

def setup_optim():
    global tao       # Modify the global object. 
    global query_pos
    global query_e
    global query_force

    print_flush('Calling setup_optim')

    # Set TAO settings. 
    tao.setType(tao_solver_type)
    tao.setObjectiveGradient(tao_objective_and_gradient, fm)
    tao.setMonitor(tao_monitor)
    tao.setMaximumIterations(max_steps)
    tao.setTolerances(gatol=max_tol) 

def read_from_progress():
    # TODO: Set iter_offset, and read other variables calculated so far. 
    pass

def write_esd_to_h5():
    print_flush('Calling write_esd_to_h5')

    # Write output hdf5 file. 
    with h5py.File('esd.h5', 'w') as f:
        ds_pos = f.create_dataset('positions', (max_steps, nat*3), dtype='f8')
        ds_e = f.create_dataset('energies', (max_steps,), dtype='f8')
        ds_gwbse_e = f.create_dataset('gwbse_energies', (max_steps,), dtype='f8')
        ds_force = f.create_dataset('forces', (max_steps, nat*3), dtype='f8')
        ds_amass = f.create_dataset('atomic_masses', (max_steps, nat), dtype='f8')
        ds_iterations = f.create_dataset('done_iterations', (1,), dtype='i4')

        ds_pos[:] = xctd_pos[:].reshape(max_steps, nat*3)
        ds_e[:] = xctd_energy[:]
        ds_gwbse_e[:] = gwbse_energy[:]
        ds_force[:] = xctd_force[:].reshape(max_steps, nat*3)
        ds_amass[:] = xctd_amass[:]
        ds_iterations[0] = tao.its

def write_esd_to_xsf():
    # TODO: write to XSF
    print_flush('Calling write_esd_to_xsf')
    pass

# endregion. 



# region: Main code. 

# Set flag. 
set_read_json_flag(True)    # To update input positions using json file. 

# TAO solve. 
read_from_progress()
setup_optim()       # Sets all the optimization parameters. 
tao.solve(pos)
print(f'TAO final position output: {pos[:]}')

# Write outputs. HDF5 and XSF files. 
write_esd_to_h5()
write_esd_to_xsf()

# Clear flag. 
set_read_json_flag(False)   # To set to default: not update positions using json file. 

# endregion. 

# region: Stop global timer.
global_logger.time_stop()
print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 
# endregion. 