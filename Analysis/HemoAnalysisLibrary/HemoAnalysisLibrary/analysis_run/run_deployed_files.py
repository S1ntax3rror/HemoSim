import os
import subprocess


def submit_analysis(file_path):
    # Submit run script file
    task = subprocess.run(['sbatch', file_path], capture_output=True)

    # Read slurm id
    tskid = int(task.stdout.split()[-1])
    return tskid


def get_subdirs(path):
    files = os.listdir()
    dirs = []
    for file in files:
        if os.path.isdir(os.path.join(path, file)):
            dirs.append(os.path.join(path, file))
    return dirs


if __name__ == "__main__":
    os.chdir("../..")
    cwd = os.getcwd()
    sub_dirs = get_subdirs(cwd)

    run_analysis = "run_analysis.sh"

    for directory in sub_dirs:
        submit_analysis(os.path.join(directory, run_analysis))
