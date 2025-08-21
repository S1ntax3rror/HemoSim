import os
import shutil
from distutils.dir_util import copy_tree

# shutil.copy(crd_path, os.path.join(new_setup_path, "init_setup.crd"))
# copy_tree(template_setup_folder, new_setup_path)

# permutations = ["1110", "1011", "0111"]
permutations = ["1110", "1101", "0111", "1011"]
# permutations = ["1110", "0111"]
perm_prefix = "/cluster/home/kaeserj/HemoSimulations/a_Permutations/Simulation_container/branches_of_main/10ns"
clean_setup_path = "/cluster/home/kaeserj/HemoSimulations/a_Permutations/V7_ShortRun_NVT/clean"
# specify which frame to take in 100ps steps (10 = 1ns)
# structures = [[10, 40], [23], [10, 47]]
structures = [[0, 124], [18], [95, 506], [294]]

# structures = [[0, 124], [18], [95, 506], [294]]


def gen_simulation_setup(location, clean_setup, run_files):
    copy_tree(clean_setup, location)
    for file in run_files:
        shutil.copy(file["data"], file["location"])


cwd = os.getcwd()

for perm, struct_list in zip(permutations, structures):
    for struct in struct_list:
        perm_path = os.path.join(perm_prefix, "permutation"+perm)
        struct_perm_path = os.path.join(cwd, f"perm{perm}struct{struct}")
        os.mkdir(struct_perm_path)

        for i in range(100):
            location_path = os.path.join(struct_perm_path, f"run{i}")
            os.mkdir(location_path)

            run_setup_str = {"data": os.path.join(perm_path, "run_setup.str"),
                             "location": os.path.join(location_path, "run_setup.str")}
            run_setup_psf = {"data": os.path.join(perm_path, "run_setup.psf"),
                             "location": os.path.join(location_path, "run_setup.psf")}
            run_setup_crd = {"data": os.path.join(perm_path, f"step5_{struct}.crd"),
                             "location": os.path.join(location_path, "run_setup.crd")}

            run_file_list = [run_setup_psf, run_setup_str, run_setup_crd]

            print(f"Generating setup for structure {struct}. Run {i}")
            # print("Generating setup with parameters", location_path, "\n", clean_setup_path, "\n", run_file_list, "\n\n")

            gen_simulation_setup(location_path, clean_setup_path, run_file_list)
