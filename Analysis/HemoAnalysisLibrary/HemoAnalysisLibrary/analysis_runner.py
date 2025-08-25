import os
import subprocess


def get_subdirs(path):
    files = os.listdir()
    dirs = []
    for file in files:
        if os.path.isdir(os.path.join(path, file)):
            dirs.append(os.path.join(path, file))
    return dirs


def start_multiprocessing(file_name):
    cwd = os.getcwd()
    subdirs = get_subdirs(cwd)

    for subdir in subdirs:
        print(f"Staring to run {file_name}")
        subprocess.Popen(["python3", os.path.join(subdir, file_name)], cwd=subdir)


if __name__ == "__main__":
    start_multiprocessing("fe_CO_direct_analysis_with_angle.py")
