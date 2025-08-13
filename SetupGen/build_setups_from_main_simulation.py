import os
import shutil
from distutils.dir_util import copy_tree
"""
This file needs to be placed into the directory where the setup folders should be generated

Requires:
Path to simulation folder
Prefix/suffix of crd files
Path to PSF from sim folder
Path to template setup folder
"""

debug = False
main_sim_folder_path = "/cluster/home/kaeserj/HemoSimulations/a_Permutations/all_perm_unbound_with_ghosts/permutation1111"
prefix_crd = "step5_"
suffix_crd = ".crd"
template_setup_folder = "/cluster/home/kaeserj/HemoSimulations/a_Permutations/V3_Kai_setup"
psf_short_path = "run_setup.psf"
i_files = [10, 20, 30, 80, 130, 180]

if debug:
    i_files = [10]

psf_path = os.path.join(main_sim_folder_path, psf_short_path)

print("Using ", main_sim_folder_path, "structure source")
print("Using ", template_setup_folder, "template source")

for i in i_files:
    print("Starting creation of setup folder ", i)
    # setup paths and mkdir
    crd_path = os.path.join(main_sim_folder_path, prefix_crd + str(i) + suffix_crd)
    new_setup_path = "hemo_setup_" + str(i) + "ns"

    if not os.path.exists(new_setup_path):
        os.makedirs(new_setup_path)
    else:
        shutil.rmtree(new_setup_path)
        os.makedirs(new_setup_path)

    # copy template setup into new setup
    copy_tree(template_setup_folder, new_setup_path)

    # replace crd and psf files
    if os.path.exists(os.path.join(new_setup_path, "init_setup.crd")):
        os.remove(os.path.join(new_setup_path, "init_setup.crd"))
    if os.path.exists(os.path.join(new_setup_path, "init_setup.psf")):
        os.remove(os.path.join(new_setup_path, "init_setup.psf"))

    shutil.copy(crd_path, os.path.join(new_setup_path, "init_setup.crd"))
    shutil.copy(psf_path, os.path.join(new_setup_path, "init_setup.psf"))






