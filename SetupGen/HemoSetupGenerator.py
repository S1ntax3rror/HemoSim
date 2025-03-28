import os
import shutil
from gettext import lngettext
from os import mkdir
from distutils.dir_util import copy_tree


def generate_all_permutations_unbound(active_directory, base_setup_directory):
    line1 = "set co1"
    line2 = "set co2"
    line3 = "set co3"
    line4 = "set co4"

    for bond_1 in range(2):
        for bond_2 in range(2):
            for bond_3 in range(2):
                for bond_4 in range(2):
                    print("setting up folder for permuation:: ", bond_1, bond_2, bond_3, bond_4)
                    current_dir = active_directory + "/" + "permutation" + str(bond_1) + str(bond_2) + str(
                        bond_3) + str(bond_4)
                    mkdir(current_dir)
                    copy_tree(base_setup_directory, current_dir)
                    path_step4 = current_dir + "/" + "step_4_generate_valid_setup.inp"
                    with open(path_step4, "r") as f:
                        text = f.readlines()
                        text_copy = text.copy()
                        for i, line in enumerate(text):
                            if line1 in line:
                                if bond_1 == 0:
                                    text_copy[i] = line1 + " 1\n"
                                else:
                                    text_copy[i] = line1 + " 0\n"
                            if line2 in line:
                                if bond_2 == 0:
                                    text_copy[i] = line2 + " 1\n"
                                else:
                                    text_copy[i] = line2 + " 0\n"
                            if line3 in line:
                                if bond_3 == 0:
                                    text_copy[i] = line3 + " 1\n"
                                else:
                                    text_copy[i] = line3 + " 0\n"
                            if line4 in line:
                                if bond_4 == 0:
                                    text_copy[i] = line4 + " 1\n"
                                else:
                                    text_copy[i] = line4 + " 0\n"
                        f.close()
                    with open(path_step4, "w") as f:
                        f.writelines(text_copy)
                        f.close()


if __name__ == "__main__":

    overwrite = False
    live = True

    cwd = os.getcwd()

    if "cluster" in cwd:
        live = True

    hemo_base_path = cwd + "/" + "HemoBase"
    path_setup_destination = "all_perm_unbound_with_ghosts"

    if live:
        hemo_base_path = cwd + "/" + "V3_Kai_setup"

    if os.path.exists(cwd + "/" + path_setup_destination) and (not overwrite and not live):
        print("EITHER SET OVERWRITE TRUE OR CHANGE PATH")
        exit()
    elif os.path.exists(cwd + "/" + path_setup_destination):
        shutil.rmtree(cwd + "/" + path_setup_destination)

    mkdir(cwd + "/" + path_setup_destination)
    active_dir = cwd + "/" + path_setup_destination

    generate_all_permutations_unbound(active_dir, hemo_base_path)
