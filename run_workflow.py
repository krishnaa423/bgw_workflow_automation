import os 
import time 
import subprocess


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

# region: 1.1-scf
os.chdir('./1.1-scf')
# run_cmd('./jobscript_pwscf.run')
os.chdir('../')
# endregion

# region: 1.2-ph
os.chdir('./1.2-ph')
# run_cmd('./jobscript_kgrid.run')
# run_cmd('./jobscript_pwbands.run')
# run_cmd('./jobscript_ph.run')
# run_cmd('./jobscript_q2r.run')
# run_cmd('./jobscript_matdyn.run')
# run_cmd('./jobscript_matdyn_bands.run')
os.chdir('../')
# endregion

# region: 1.3-epw
os.chdir('./1.3-epw')
# run_cmd('./jobscript_epw.run')
os.chdir('../')
# endregion

# region: 1.5-bands
# os.chdir('./1.5-bands')
# run_cmd('./jobscript_pwbands.run')
# run_cmd('./jobscript_pw2bgw.run')
# run_cmd('./jobscript_parabands.run')
# run_cmd('./jobscript_bands.run')
# os.chdir('../')
# endregion

# region: 1.6-dos
# os.chdir('./1.6-dos')
# run_cmd('./jobscript_dos.run')
# os.chdir('../')
# endregion

# region: 1.7-pdos
# os.chdir('./1.7-dos')
# run_cmd('./jobscript_dos.run')
# os.chdir('../')
# endregion

# region: 1.8-pp
# os.chdir('./1.8-pp')
# run_cmd('./jobscript_pp_charge_density.run')
# os.chdir('../')
# endregion

# region: 2.1-wfn
os.chdir('./2.1-wfn')
# run_cmd('./jobscript_kgrid.run')
# run_cmd('./jobscript_pwbands.run')
# run_cmd('./jobscript_pw2bgw.run')
run_cmd('./jobscript_parabands.run')
os.chdir('../')
# endregion


# region: 2.2-wfnq
os.chdir('./2.2-wfnq')
run_cmd('./jobscript_kgrid.run')
run_cmd('./jobscript_pwbands.run')
run_cmd('./jobscript_pw2bgw.run')
# run_cmd('./jobscript_parabands.run')
os.chdir('../')
# endregion


# region: 3.1-wfn_fi
# os.chdir('./3.1-wfn_fi')
# run_cmd('./jobscript_kgrid.run')
# run_cmd('./jobscript_pwbands.run')
# run_cmd('./jobscript_pw2bgw.run')
# run_cmd('./jobscript_parabands.run')
# os.chdir('../')
# endregion


# region: 3.2-wfnq_fi
# os.chdir('./3.2-wfnq_fi')
# run_cmd('./jobscript_kgrid.run')
# run_cmd('./jobscript_pwbands.run')
# run_cmd('./jobscript_pw2bgw.run')
# os.chdir('../')
# endregion


# region: 4-epsilon
os.chdir('./4-epsilon')
run_cmd('./jobscript_epsilon.run')
os.chdir('../')
# endregion


# region: 5.1-sigma
os.chdir('./5.1-sigma')
run_cmd('./jobscript_sigma.run')
os.chdir('../')
# endregion

# region: 5.2-inteqp
# os.chdir('./5.2-inteqp')
# run_cmd('./jobscript_inteqp.run')
# os.chdir('../')
# endregion


# region: 6-kernel
os.chdir('./6-kernel')
run_cmd('./jobscript_kernel.run')
os.chdir('../')
# endregion


# region: 7-absorption
os.chdir('./7-absorption')
run_cmd('./jobscript_absorption.run')
# run_cmd('./jobscript_plotxct.run')          # If you want to plot the exciton too. 
os.chdir('../')
# endregion


# region: 8-esf
os.chdir('./8-esf')
run_cmd('./jobscript_esf.run')
os.chdir('../')
# endregion

# region: Stop global timer.
global_logger.time_stop()
print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 
# endregion. 