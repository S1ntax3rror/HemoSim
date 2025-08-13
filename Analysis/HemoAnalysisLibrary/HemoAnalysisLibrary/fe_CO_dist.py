import os
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob


def plot_simulation(save_location="", prefix=""):
    cwd = os.getcwd()
    if "cluster" in cwd and prefix == "":
        prefix = "../permutation0000"
    else:
        if prefix == "":
            prefix = "test_data"

    dcd_files = sorted(glob.glob(os.path.join(prefix, "step5_*.dcd")))

    print(dcd_files)

    u = mda.Universe(
        os.path.join(prefix, "run_setup.psf"),
        dcd_files
    )

    # heme_segids = sorted(set(atom.segid for atom in u.select_atoms("resname HEME")))
    CO_segids = ["COA", "COB", "COC", "COD"]
    heme_segids = ['HETA', 'HETB', 'HETC', 'HETD']


    def get_CO_positions(CO_seg_atom_list):
        if isinstance(CO_seg_atom_list, list):
            return [atom[1].position for atom in CO_seg_atom_list]
        atom = u.select_atoms(f"segid {CO_seg_atom_list}")
        return [atom[0].position, atom[1].position]


    print(CO_segids)
    print(heme_segids)

    results = []
    for ts in u.trajectory:
        frame_data = {'frame': ts.frame}
        dist_array = np.zeros(4)
        cnt = 0

        for CO_segid, heme_segid in zip(CO_segids, heme_segids):
            heme = u.select_atoms(f"resname HEME and segid {heme_segid}")

            try:
                fe = heme.select_atoms(f"name FE")
                pos_fe = fe[0].position

                # CO atoms
                CO_positions = get_CO_positions(CO_segid)
                CO_position = np.mean(CO_positions, axis=0)

                dist_array[cnt] = np.linalg.norm(pos_fe-CO_position)
                if dist_array[cnt] < 24:
                    frame_data[CO_segid + "_dist"] = dist_array[cnt]


            except IndexError:
                print(f"Missing atoms in segid {heme_segid} at frame {ts.frame}")
                continue

            cnt += 1

        results.append(frame_data)

    df = pd.DataFrame(results)
    df.to_csv("heme_fe_geometry_all4.csv", index=False)

    for segid in CO_segids:
        col = f"{segid}_dist"  # change this if you want to plot something else such as distances, angles
        if col in df.columns:
            plt.plot(df['frame'], df[col], label=f'{segid}', linewidth=0.3)

    plt.axhline(0, color='gray', linestyle='--', linewidth=0.2)
    plt.title("Fe Out-of-Plane Displacement for All Heme Groups")
    plt.xlabel("Frame")  # can be changed to time
    plt.ylabel("Displacement (Å)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    if save_location == "":
        plt.savefig("all_fe_out_of_plane.png", dpi=300)
    else:
        plt.savefig(save_location, dpi=300)
    if "test_data" in prefix:
        plt.show()


if __name__ == "__main__":
    cwd = os.getcwd()
    if "cluster" in cwd:
        plot_simulation()
    else:
        ns_folders = os.listdir(cwd)
        for ns_folder in ns_folders:
            ns_folder_path = os.path.join(cwd, ns_folder)
            permutation_folders = os.listdir(ns_folder_path)

            for permutation_folder in permutation_folders:
                if "permutation" in permutation_folder:
                    permutation_folder_path = os.path.join(cwd, permutation_folder)
                    save_path = ns_folder_path + "_" + permutation_folder + ".png"
                    plot_simulation(save_location="", prefix=permutation_folder_path)
