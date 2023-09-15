import pyparsing as pp 
import argparse 
import re
import os 
import time 

parser = argparse.ArgumentParser()


# Parser add argument types. 
parser.add_argument('--kgrids', nargs='+', type=str, help='List of kgrids. Example: 2x2x2, 4x4x4')
parser.add_argument('--ecutwfns', nargs='+', type=float, help='List of wfn cutoffs. 40 50 60.')

# Parse the input arguments. 
args = parser.parse_args()


# Functions. 
def write_kgrid(kgrid, file):
    pattern_kgrid = ...+ pp.Literal('K_POINTS automatic\n') + ...+ pp.Literal(pp.srange('[a-zA-Z]')) + pp.Regex('.*', re.DOTALL)

def write_ecutwfn(ecut, file):
    pattern_ecutwfn = ...+ pp.Literal('ecutwfn') + pp.Literal('=') + ...+ pp.Literal(pp.srange('[a-zA-Z/]')) + pp.Regex('.*', re.DOTALL)

def get_kgrid(file):
    pass

def get_ecutwfn(file):
    pass

def copy_output_kgrid_run(kgrid):
    
    dir_name = f'kgrid_{kgrid}'

    # Create new directory. 
    os.mkdir(f'./{dir_name}')

    # Copy the output to new directory. 
    os.system(f'cp ../pwscf.out ./{dir_name}')

def copy_output_ecutwfn_run(ecut):
    dir_name = f'ecut_{ecut}'

    # Create new directory. 
    os.mkdir(f'./{dir_name}')

    # Copy the output to new directory. 
    os.system(f'cp ../pwscf.out ./{dir_name}')


# Variables. 
pwscf_input_file = '../pwscf.in'
pwscf_output_file = '../pwscf.out'
cached_kgrid = get_kgrid(pwscf_input_file)
cached_ecutwfn = get_ecutwfn(pwscf_input_file)


# region: Workflow code. 

# region: Global variables. 
sched_type = 'WSL'   # WSL or slurm. 
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

    # Starting a command. 
    print_flush(f'Starting : {cmd}')
    exit_code = -1

    # Run the command. 
    if is_time_logging: local_logger.time_start()

    if sched_type == 'WSL':
        exit_code = os.system(f'{cmd}')
    elif sched_type == 'slurm':
        exit_code = os.system(f'sbatch --wait {cmd}')

    if is_time_logging: local_logger.time_stop()

    # See if it finished properly. 
    if exit_code == 0:  # No error. 
        # Done with a command. 
        if is_time_logging:
            print_flush(f'Done with: {cmd}. Took {local_logger.time_elapsed()} seconds.\n\n')
        else:
            print_flush(f'Done with: {cmd}.\n\n')
    else:
        # Error. 
        if is_time_logging:
            print_flush(f'Error finishing: {cmd}. Exited with code {exit_code}. Took time {local_logger.time_elapsed()}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        else:
            print_flush(f'Error finishing: {cmd}. Exited with code {exit_code}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the dft convergence run is: {global_logger.time_elapsed()} seconds.') 

        os._exit(exit_code)     # Exit with the error code. 
# endregion

# region: Start global timer.
global_logger.time_start()
# endregion

# endregion

# region: Iterate through kgrids.
for kgrid_str in args.kgrids:
    kgrid_int = map(int, kgrid_str.split('x'))

    # Write the new kgrid. 
    write_kgrid(kgrid_int, pwscf_input_file)

    # region: 1.1-scf
    os.chdir('../')
    run_cmd('./jobscript_pwscf.run')
    os.chdir('./convergence')
    write_kgrid(cached_kgrid, pwscf_input_file)         # Write the old one back. 
    copy_output_kgrid_run(kgrid_str)
    # endregion

# endregion. 

# region: Iterate through ecutwfns.
for ecutwfn_float in args.ecutwfns:

    # Write the new ecutwfn. 
    write_ecutwfn(ecutwfn_float, pwscf_input_file)

    # region: 1.1-scf
    os.chdir('../')
    run_cmd('./jobscript_pwscf.run')
    os.chdir('./convergence')
    write_kgrid(cached_ecutwfn, pwscf_input_file)         # Write the old one back. 
    copy_output_ecutwfn_run(ecutwfn_float)
    # endregion

# endregion. 



# region: Stop global timer.
global_logger.time_stop()
print_flush(f'Total time elapsed for the dft convergence run is: {global_logger.time_elapsed()} seconds.') 
# endregion. 