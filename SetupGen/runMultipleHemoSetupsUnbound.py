import subprocess
import os
import sys
from multiprocessing import Pool


def replace_line_in_file(file, old_line, new_line):
    with open(file, "r") as f:
        lines = f.readlines()
        lines_copy = lines.copy()
        for i, line in enumerate(lines):
            if old_line in line:
                lines_copy[i] = new_line
        f.close()
    with open(file, "w") as f:
        f.writelines(lines_copy)
        f.close()


def observe_run(path):
    result = subprocess.run(['python3', 'Observer.py'], capture_output=True, text=True, cwd=path)
    text = result.stdout
    text_err = result.stderr

    with open("outfile_sh.out", "a") as file:
        file.write(text)
        file.write(text_err)


if __name__ == "__main__":
    active_dir = "/cluster/home/kaeserj/HemoSimulations/a_Permutations/Simulation_container/branches_of_main"
    print("active directory:", active_dir[-20:])

    perm_dir_paths = os.listdir(active_dir)
    clean_perm_dirs = []
    for perm_dir in perm_dir_paths:
        if os.path.isdir(os.path.join(active_dir, perm_dir)):
            clean_perm_dirs.append(os.path.join(active_dir, perm_dir))

    sh_file_paths = []

    for clean_dir in clean_perm_dirs:
        paths = os.listdir(clean_dir)
        for path_ in paths:
            run_path = os.path.join(clean_dir, path_) + "/"
            if os.path.isdir(run_path) and not os.path.exists(os.path.join(run_path, "step5_999.pdb")):
                sh_file_paths.append(run_path)
                print("run dictionary:", clean_dir[-25])


    with Pool(processes=len(sh_file_paths)) as pool:
        pool.map(observe_run, sh_file_paths)
