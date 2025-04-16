#!/usr/bin/python

# Basics
import os
import subprocess
import sys

import numpy as np

# Miscellaneous
import time
import datetime
from glob import glob
import getpass

# CHARMM input, output and slurm submission files
runfile = "run.sh"  # Bash script to run the setup
setupfile = "create_run_setup.sh"  # Bash script to initialize the setup
inpfile = "gpu_step5_production.inp"  # Input file for charmm
outfile = "gpu_step5_production.out"  # Output file for charmm
run_checker_file = "run_setup.psf"  # if this file exists then don't run setup

# CHARMM dcd file name tag
dcdftag = "step5_*.dcd"
dcdsplt = [["step5_", -1], [".", 0]]

# Input file line to replace
# first: string to identify line; second: new line accepting .format(new_step)
repline_step = ["set cnt", "set cnt {:d}"]

# Input file line to replace if equi run complete
# first: file to look at, if exist, apply replace line;
# second: string to identify line; third: new line
repline_equi = ["step5_0.pdb", "set equi", f"set equi 0"]

# Error line (or part of line) to look in 'outfile'.
# If found restart run at latest completed step.
# Else, don't restart
# If error line is None, ignore this check
# errline = "Energy is NaN"
errline = None

# Sleep time between status checks
slptime = 30

# Stop flag file, if in directory, stop observation
stpfile = "stop"


def log_2_file(message):

    print("logging", message)

    log_path = os.path.join(os.getcwd(), "observer_log.txt")
    if os.path.exists(log_path):
        with open(log_path, "a") as logfile:
            logfile.write(message + "\n")
        logfile.close()
        return

    with open(log_path, "w") as logfile:
        logfile.write(message + "\n")
    logfile.close()
    return


# Function to apply changes to the input file
def update_inpfile(new_step):
    # Read input file
    with open(inpfile, 'r') as f:
        inplines = f.read()

    # Check equi run
    if os.path.exists(repline_equi[0]):
        for line in inplines.split("\n"):
            if repline_equi[1] in line:
                inplines = inplines.replace(line, repline_equi[2])

    # Replace step number
    for line in inplines.split("\n"):
        if repline_step[0] in line:
            inplines = inplines.replace(line, repline_step[1].format(new_step))

    # Write modified input file
    with open(inpfile, 'w') as f:
        f.write(inplines)

    return


# Function to submit simulation
def submit_setupfile():
    # Submit run script file
    task = subprocess.run(['sbatch', setupfile], capture_output=True)

    # Read slurm id
    tskid = int(task.stdout.split()[-1])

    return tskid


# Function to submit simulation
def submit_runfile():
    # Submit run script file
    task = subprocess.run(['sbatch', runfile], capture_output=True)

    # Read slurm id
    tskid = int(task.stdout.split()[-1])

    return tskid


# Function to find new step number
def find_new_step():
    # Get all dcd files
    dcdfiles = glob(dcdftag)

    # Check steps from dcd files, if not defined new step is 1
    if len(dcdfiles):

        dcdidcs = np.zeros_like(dcdfiles, dtype=int)
        for idcd, dcdfile in enumerate(dcdfiles):
            idx = dcdfile[:]
            for isplt in dcdsplt:
                idx = idx.split(isplt[0])[isplt[1]]
            dcdidcs[idcd] = int(idx)

        # Get highest dcd index as new step
        new_step = np.max(dcdidcs)

    else:

        new_step = 1

    return new_step


# Function to check current job status
def check_status(tskid):
    # Get task id list
    user = getpass.getuser()
    tsklist = subprocess.run(['squeue', '-u', user], capture_output=True)
    idslist = [
        int(ids.split()[0])
        for ids in tsklist.stdout.decode().split('\n')[1:]
        if len(ids.strip())]

    # Return if job is still running
    if tskid in idslist:
        return tskid
    else:
        log_2_file("")
        log_2_file(str(datetime.datetime.now()) + f": Job {tskid:d} done.")

    if not os.path.exists(os.path.join(os.getcwd(), outfile)):
        return tskid

    # Read output file
    with open(outfile, 'r') as f:
        outlines = f.readlines()

    # If not, check for normal termination
    for line in outlines[::-1]:
        if " NORMAL TERMINATION" in line: # keep space to avoid going here through keywords ABNORMAL TERMINATION
            log_2_file(f": Job {tskid:d} finished.")
            return 0
    log_2_file(f": Job {tskid:d} broken. Check restart")

    # If not, check for error line in output file

    # # Check for error line
    # if errline is None:
    #     errline_found = True
    # else:
    #     errline_found = False
    #     for line in outlines[::-1]:
    #         if errline in line:
    #             errline_found = True
    #             break
    #
    # # If error line, found start resubmission, else return
    # if not errline_found:
    #     log_2_file("Error not found! Job is done")
    #     return 0
    # else:
    #     log_2_file("Error found! Restart job")

    # Find new step
    new_step = find_new_step()
    log_2_file(f"  New step: {new_step:d}")

    # Update input file with new step
    update_inpfile(new_step)

    # Resubmit simulation
    new_tskid = submit_runfile()
    log_2_file(f"  New task id: {new_tskid:d}")

    # Return new task id
    return new_tskid


# Function to run observation loop
def observe_loop(tskid):
    # Infinite loop until job is finished
    finished = False
    current_tskid = tskid
    while not finished:

        # Wait time
        time.sleep(slptime)

        # Check job status
        current_tskid = check_status(current_tskid)

        # Check for stop flag
        if os.path.exists(stpfile):
            log_2_file("Stop flag found!")
            finished = True

        # If tskid is zero, job is finished
        if current_tskid == 0:
            log_2_file("Job is done!(?)")
            finished = True
    return


# Start job and observation
if __name__ == "__main__":

    # if setup didn't run yet, submit setup first
    # and wait observe till completed then run the setup
    if not os.path.exists(os.path.join(os.getcwd(), run_checker_file)):
        tskid = submit_setupfile()
        log_2_file(f"Setup and run task id: {tskid:d}")
        observe_loop(tskid)

    # Find new step
    new_step = find_new_step()
    log_2_file(f"Start at step {new_step:d}")

    # Update input file with new step
    update_inpfile(new_step)

    # Submit initial simulation setup
    tskid = submit_runfile()
    log_2_file(f"Initial task id: {tskid:d}")

    # Run observation
    observe_loop(tskid)
