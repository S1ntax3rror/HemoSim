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
    depth = 2
    analysis_file_name = "fe_CO_direct_analysis_with_angle.py"
    cwd = os.getcwd()
    sub_dirs = get_subdirs(cwd)


    for directory in sub_dirs:
        analysis_path = cwd

        submit_analysis(analysis_path)
