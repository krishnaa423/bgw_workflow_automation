import os 
import time 
import re

# for debugging. 
os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/tests/dwfn')

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
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        os._exit(exit_code)     # Exit with the error code. 
# endregion

# region: Start global timer.
global_logger.time_start()
# endregion

# region: Main code.

# region: Variables.
position = 1.15
delta = 0.01 
pos_pattern = r'positions\s=\s(?P<pos>\[.*?\])'
# endregion. 

# region: Functions. 
def pos_repl(match):

    new_string = \
    f'''positions = [
        (0.0, 0.0, 0.0),
        ({position}, 0.0, 0.0)
    ]'''
    return new_string
# endregion. 

# region: Run 

# region: Unstepped position.
# Read create workflow.py. 
position = position
with open('../../create_workflow.py', 'r') as f: wflow = f.read()
# Substitute new position in it. 
new_wflow = re.sub(pos_pattern, pos_repl, wflow, flags=re.DOTALL)
with open('../../create_workflow.py', 'w') as f: f.write(new_wflow)
# Run. 
os.chdir('../../')
run_cmd('python remove_workflow.py')
run_cmd('python create_workflow.py')
run_cmd('python run_workflow.py > run_workflow.log')
os.chdir('./tests/dwfn/')
# Copy. 
os.system(f'cp ../../8-esf/WFN.h5 ./WFN_{position}.h5')
os.system(f'cp ../../8-esf/elph_1.h5 ./elph_1_{position}.h5')
os.system(f'cp ../../8-esf/kgrid.out ./kgrid.out')
# endregion. 


# region: Stepped position. 
# Read create workflow.py. 
position = position + delta
with open('../../create_workflow.py', 'r') as f: wflow = f.read()
# Substitute new position in it. 
new_wflow = re.sub(pos_pattern, pos_repl, wflow, flags=re.DOTALL)
with open('../../create_workflow.py', 'w') as f: f.write(new_wflow)
# Run. 
os.chdir('../../')
run_cmd('python remove_workflow.py')
run_cmd('python create_workflow.py')
run_cmd('python run_workflow.py > run_workflow.log')
os.chdir('./tests/dwfn/')
# Copy. 
os.system(f'cp ../../8-esf/WFN.h5 ./WFN_{position}.h5')
# endregion.

# endregion.

# endregion. 

# region: Stop global timer.
global_logger.time_stop()
print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 
# endregion. 

