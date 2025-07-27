#!/usr/bin/python

# Basics
import os
import subprocess
import numpy as np

# Miscellaneous
import time
import datetime
from glob import glob
import getpass

log_file = "observer.log"
if not os.path.exists(log_file):
    with open(log_file, "w") as f:
        f.write("Created Logger")
        f.close()

# CHARMM input, output and slurm submission files
runfile = "run.sh"
setupfile = "create_run_setup.sh"
inpfile = "npt2_step5.inp"
outfile = "npt2_step5.out"

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


def log(message):
    with open(log_file, "a") as f:
        f.write(message)
        f.close()

# Function to apply changes to the input file
def update_inpfile(new_step_):
    # Read input file
    with open(inpfile, 'r') as f:
        inplines = f.read()
        f.close()

    # Check equi run
    if os.path.exists(repline_equi[0]):
        for line in inplines.split("\n"):
            if repline_equi[1] in line:
                inplines = inplines.replace(line, repline_equi[2])

    # Replace step number
    for line in inplines.split("\n"):
        if repline_step[0] in line:
            inplines = inplines.replace(line, repline_step[1].format(new_step_))

    # Write modified input file
    with open(inpfile, 'w') as f:
        f.write(inplines)
        f.close()

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
        log(f"{datetime.datetime.now()} : Job {tskid} still running. Sleeping")
        return tskid
    else:
        log("")
        log(f"{datetime.datetime.now()} : Job {tskid:d} done.")

    if not os.path.exists(os.path.join(os.getcwd(), outfile)):
        return tskid

    # Read output file
    with open(outfile, 'r') as f:
        outlines = f.readlines()

    # If not, check for normal termination
    for line in outlines[::-1]:
        if "NORMAL TERMINATION" in line:
            log(f": Job {tskid:d} finished.")
            return 0
    log(f": Job {tskid:d} broken. Check restart")

    # If not, check for error line in output file

    # Check for error line
    if errline is None: # If this setting is chosen, then there is no option besides NORMAL TERMINATION to allow exit
        errline_found = True
    else:
        errline_found = False
        for line in outlines[::-1]:
            if errline in line:
                errline_found = True
                break

    # If error line, found start resubmission, else return
    if not errline_found:
        log("Error not found! Job is done")
        return 0
    else:
        log(f"{datetime.datetime.now()} : Error found! Restart job")

    # Find new step
    new_step = find_new_step()
    log(f"  New step: {new_step:d}")

    # Update input file with new step
    update_inpfile(new_step)

    # Resubmit simulation
    new_tskid = submit_runfile()
    log(f"  New task id: {new_tskid:d}")

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
            log("Stop flag found! Most likely TERMINATION due to CANCELATION")
            finished = True

        # If tskid is zero, job is finished
        if current_tskid == 0:
            print("Job is done!(?)")
            finished = True
    return


# Start job and observation
if __name__ == "__main__":
    # Find new step
    new_step = find_new_step()
    log(f"{datetime.datetime.now()} : Start at step {new_step:d}")

    # Update input file with new step
    update_inpfile(new_step)

    # Submit initial simulation setup
    tskid = submit_setupfile()
    log(f"{datetime.datetime.now()} : Initial task id: {tskid:d}")

    # Run observation
    observe_loop(tskid)