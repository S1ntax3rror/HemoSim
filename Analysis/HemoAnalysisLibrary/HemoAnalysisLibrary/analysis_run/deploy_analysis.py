import os
import shutil


def get_subdirs(path):
    files = os.listdir()
    dirs = []
    for file in files:
        if os.path.isdir(os.path.join(path, file)):
            dirs.append(os.path.join(path, file))
    return dirs


if __name__ == "__main__":
    cwd = os.getcwd()
    analysis_script = "analysis_runner.py"
    analysis_script_path = os.path.join(cwd, analysis_script)
    analysis_exec = "run_analysis.sh"
    analysis_exec_path = os.path.join(cwd, analysis_exec)
    os.chdir("..")
    cwd_backed = os.getcwd()
    directories = get_subdirs(cwd_backed)
    for directory in directories:
        if os.path.isdir(directory) and "perm" in directory:
            print(f"copying files to {directory}")
            shutil.copy(analysis_script_path, os.path.join(directory, analysis_script))
            shutil.copy(analysis_exec_path, os.path.join(directory, analysis_exec))
