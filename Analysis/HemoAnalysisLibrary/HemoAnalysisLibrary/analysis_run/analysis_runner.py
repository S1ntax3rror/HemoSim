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

    procs = []
    for subdir in subdirs:
        print(f"Staring to run {os.path.join(subdir, file_name)}")
        p = subprocess.Popen(["python3", os.path.join(subdir, file_name)], cwd=subdir)
        procs.append(p)

    for p in procs:
        p.wait()


if __name__ == "__main__":
    start_multiprocessing("../fe_CO_direct_analysis_with_angle.py")
