import os
from os import mkdir

from HemoSetupGenerator import run_generation

""" REQUIRES HemoSetupGenerator """

def generateAllPermuationsForAllBranches(initial_condition_container_path, destination_container_name):
    cwd = os.getcwd()
    full_destination_container_path = os.path.join(cwd, destination_container_name)
    full_initial_condition_container_path = os.path.join(cwd, initial_condition_container_path)
    if os.path.exists(full_destination_container_path):
        print("ERROR: destination container path allready exists")
        print(full_destination_container_path)
        exit()

    if not os.path.exists(full_initial_condition_container_path):
        print("ERROR: initial condition container path points to nothing")
        print("GENERATE INITIAL SETUPS AND PASS THE PATH AS FIRST ARG")
        exit()

    mkdir(full_destination_container_path)

    potential_initial_conditions = os.listdir(initial_condition_container_path)
    initial_condition_folder_paths = []

    for potential_initial_condition in potential_initial_conditions:
        if os.path.isdir(os.path.join(initial_condition_container_path, potential_initial_condition)) and "hemo" in potential_initial_condition:
            initial_condition_folder_paths.append(os.path.join(initial_condition_container_path, potential_initial_condition))

    for initial_condition_path in initial_condition_folder_paths:
        dest_path_suffix = initial_condition_path.split("/")[-1]
        dest_path = os.path.join(full_destination_container_path, dest_path_suffix)

        run_generation(initial_condition_path, dest_path)


if __name__ == "__main__":
    init_cond_container_path = "init_cond_container"
    dest_container_path = "simulations_container"

    generateAllPermuationsForAllBranches(init_cond_container_path, dest_container_path)