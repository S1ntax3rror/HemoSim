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
    args = sys.argv
    if len(args) > 1 and args[1] == "setup":
        result = subprocess.run(['python3', 'Observer.py', 'setup'], capture_output=True, text=True, cwd=path)
    else:
        result = subprocess.run(['python3', 'Observer.py'], capture_output=True, text=True, cwd=path)
    text = result.stdout
    text_err = result.stderr

    with open("outfile_sh.out", "w") as file:
        file.write(text)
        file.write(text_err)

if __name__ == "__main__":
    cwd = os.getcwd()
    active_dir = cwd + "/" + "all_perm_unbound_with_ghosts/"
    print("active directory:", active_dir)

    paths = os.listdir(active_dir)

    sh_file_paths = []

    for path_ in paths:
        run_path = os.path.join(active_dir, path_) + "/"
        if os.path.isdir(run_path) and not os.path.exists(os.path.join(run_path, "step5_999.pdb")):
            sh_file_paths.append(run_path)
            print("run dictionary:", path_)

    with Pool(processes=len(sh_file_paths)) as pool:
        pool.map(observe_run, sh_file_paths)
