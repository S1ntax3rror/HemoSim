import os
import shutil
import sys
from os import mkdir
from distutils.dir_util import copy_tree

"""
var hemo_base_path: needs to be set
var path_setup_destination: can be chosen or left at default


"""



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


def gen_setups(hemo_base, dest_path):
    if os.path.exists(cwd + "/" + dest_path) and (not overwrite and not live):
        print("EITHER SET OVERWRITE TRUE OR CHANGE PATH")
        exit()
    elif os.path.exists(cwd + "/" + dest_path):
        shutil.rmtree(cwd + "/" + dest_path)

    mkdir(cwd + "/" + dest_path)
    active_dir = cwd + "/" + dest_path

    generate_all_permutations_unbound(active_dir, hemo_base)


def gen_multiple_setups(hemo_base, multiple_dests):
    # get all setup folders (hemo_bas_path should be set to parent dict of setups)
    # generate all destinations:
    setup_folders_uncleaned = os.listdir(hemo_base)
    print(setup_folders_uncleaned)
    setup_folders = []
    dest_parent_path = os.path.join(cwd, multiple_dests)
    dest_paths = []


    check1 = input(f"Confirm y/n that setup_folder dir is \n {hemo_base} \n ANSWER:  ")
    if check1 != "y":
        check2 = input(f"Enter full path to change it or enter n to close: ")
        if check2 == "n":
            exit()
        hemo_base = check2
        if not os.path.exists(hemo_base):
            print("Error path doesn't exist")
            exit()

    for setup_folder in setup_folders_uncleaned:
        if "toppar" in setup_folder:
            print("ERROR: Set the parent dictionary as hemo_base_path when using multiple setups!!")
            exit()
        print(setup_folder, os.path.isdir(setup_folder))
        if os.path.isdir(os.path.join(hemo_base, setup_folder)):
            setup_folders.append(os.path.join(hemo_base, setup_folder))
            dest_paths.append(os.path.join(dest_parent_path, setup_folder.replace("hemo_setup_", "")))


    for setup_folder, dest_path in zip(setup_folders, dest_paths):
        print("Generating all setups for ", dest_path[-10:], "using", setup_folder[-10:])
        mkdir(dest_path)

        # make sure folder doesn't exist yet to not overwrite any data
        if os.path.exists(dest_path) and (not overwrite and not live):
            print("EITHER SET OVERWRITE TRUE OR CHANGE PATH")
            exit()
        elif os.path.exists(dest_path) and (overwrite and not live):  # overwrite if forced and not live
            shutil.rmtree(dest_path)


        generate_all_permutations_unbound(dest_path, setup_folder)


if __name__ == "__main__":

    overwrite = False
    live = True
    multiple_setups = False

    cwd = os.getcwd()
    args = sys.argv

    if "cluster" in cwd:
        live = True
        overwrite = False

    hemo_base_path = cwd + "/" + "HemoBase"
    multiple_setup_destination = "Simulation_container/branches_of_main"
    path_setup_destination = "all_perm_unbound_with_ghosts"

    if live:
        if len(args) > 1:
            hemo_base_path = os.path.join(cwd, args[1])
            if not os.path.exists(hemo_base_path):
                hemo_base_path = cwd + "/" + "V4_multiple_setups"

            if len(args) > 2 and args[2] == "y":
                multiple_setups = True

        else:
            hemo_base_path = cwd + "/" + "V3_Kai_setup"


    if multiple_setups:
        print("starting multiple setups")
        gen_multiple_setups(hemo_base_path, multiple_setup_destination)
    else:
        print("starting single setup")
        gen_setups(hemo_base_path, path_setup_destination)
