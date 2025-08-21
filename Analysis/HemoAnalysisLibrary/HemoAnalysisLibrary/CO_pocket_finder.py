import importlib
import subprocess
import sys
import os

# TODO add the correct residues for each of the 4 subunits and adjust the code to only focus on the specified subunit
# TODO modify the saving process to work with the 4 subunits


def install_and_import(package):
    try:
        importlib.import_module(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    finally:
        globals()[package] = importlib.import_module(package)


# Example usage
install_and_import("numpy")
install_and_import("pandas")
install_and_import("MDAnalysis")
install_and_import("matplotlib")


import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

cwd = os.getcwd()
if "cluster" in cwd or "scratch" in cwd:
    prefix = ""
else:
    prefix = "test_data"

dcd_files = []
for i in range(1000):
    dcd_path = os.path.join(prefix, f"step5_{i}.dcd")
    if os.path.exists(dcd_path):
        dcd_files.append(dcd_path)
    else:
        break

print(dcd_files)

u = mda.Universe(
    os.path.join(prefix, "run_setup.psf"),
    dcd_files
)


def calc_pocket_pos(pocket_indices, frame_, mean_CO_pos):
    positions_pocket = frame_._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    dist_H2_pocket = np.linalg.norm(mean_CO_pos - mean_position_pocket)
    return dist_H2_pocket


def read_pocket_list(res_list, md_analysis_reference):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(md_analysis_reference.select_atoms(res).ix[0])
    return np.array(pocket_indices)


# heme_segids = sorted(set(atom.segid for atom in u.select_atoms("resname HEME")))
CO_segids = ["COA", "COB", "COC", "COD"]
heme_segids = ['HETA', 'HETB', 'HETC', 'HETD']

pocket_XE1 = ["resname LEU and resid 89 and name C", "resname HSD and resid 93 and name C",
              "resname LEU and resid 104 and name C", "resname PHE and resid 138 and name C",
              "resname ILE and resid 142 and name C", "resname TYR and resid 146 and name C"]

pocket_XE2 = ["resname LEU and resid 72 and name C", "resname ILE and resid 107 and name C",
              "resname SER and resid 108 and name C", "resname LEU and resid 135 and name C",
              "resname PHE and resid 138 and name C", "resname ARG and resid 139 and name C"]

pocket_XE3 = ["resname TRP and resid 7 and name C", "resname LEU and resid 76 and name C",
              "resname GLY and resid 80 and name C", "resname ALA and resid 134 and name C",
              "resname LEU and resid 137 and name C", "resname PHE and resid 138 and name C"]

pocket_XE4 = ["resname GLY and resid 25 and name C", "resname ILE and resid 28 and name C",
              "resname LEU and resid 29 and name C", "resname GLY and resid 65 and name C",
              "resname VAL and resid 68 and name C", "resname LEU and resid 72 and name C"]

pocket_XE5 = ["resname LEU and resid 29 and name C", "resname HEME and resid 1154 and name FE"]


# initialize arrays which will hold the distance from H2 to the corresponding pocket
num_timesteps = len(u.trajectory)

pocket_arr = np.zeros(num_timesteps)
pocket_XE1_dist_arr = np.zeros((4, num_timesteps))
pocket_XE2_dist_arr = np.zeros((4, num_timesteps))
pocket_XE3_dist_arr = np.zeros((4, num_timesteps))
pocket_XE4_dist_arr = np.zeros((4, num_timesteps))

# Get atom indices for pocket
XE1_pocket_indices = read_pocket_list(pocket_XE1, u)
XE2_pocket_indices = read_pocket_list(pocket_XE2, u)
XE3_pocket_indices = read_pocket_list(pocket_XE3, u)
XE4_pocket_indices = read_pocket_list(pocket_XE4, u)


def get_CO_positions(CO_seg_atom_list):
    if isinstance(CO_seg_atom_list, list):
        return [atom[1].position for atom in CO_seg_atom_list]
    atom = u.select_atoms(f"segid {CO_seg_atom_list}")
    return [atom[0].position, atom[1].position]

count = -1

results = []
pocket_storage = []

# for jj, frame in enumerate(dcd.trajectory[:num_timesteps]):
for ts in u.trajectory:
    count += 1
    frame = ts.frame
    frame_data = {'frame': frame}
    frame_pocket_data = {'frame': frame}
    print(f"reading for frame {ts}")
    reading_co = -1
    for CO_segid, heme_segid in zip(CO_segids, heme_segids):
        reading_co += 1
        heme = u.select_atoms(f"resname HEME and segid {heme_segid}")

        try:
            fe = heme.select_atoms(f"name FE")
            pos_fe = fe[0].position

            # CO atoms
            CO_positions = get_CO_positions(CO_segid)
            CO_position = np.mean(CO_positions, axis=0)

            pocket_XE1_dist_arr[reading_co][count] = calc_pocket_pos(XE1_pocket_indices, ts, CO_position)
            pocket_XE2_dist_arr[reading_co][count] = calc_pocket_pos(XE2_pocket_indices, ts, CO_position)
            pocket_XE3_dist_arr[reading_co][count] = calc_pocket_pos(XE3_pocket_indices, ts, CO_position)
            pocket_XE4_dist_arr[reading_co][count] = calc_pocket_pos(XE4_pocket_indices, ts, CO_position)


            # dist = np.linalg.norm(pos_fe-CO_position)
            # if dist < 6:
            #     frame_pocket_data[CO_segid + " Pocket"] = "A or B"
            # elif 6 <= dist <= 12.5:
            #     frame_pocket_data[CO_segid + " Pocket"] = "XE"
            # else:
            #     frame_pocket_data[CO_segid + " Pocket"] = "NA"
            #
            # if dist < 24:
            #     frame_data[CO_segid + "_dist"] = dist

        except IndexError:
            print(f"Missing atoms in segid {heme_segid} at frame {ts.frame}")
            continue

    # results.append(frame_data)
    # pocket_storage.append(frame_pocket_data)


for co in range(4):
    for distance in range(num_timesteps):
        frame_pocket_data[CO_segid + " Pocket"] = "A/B"


df = pd.DataFrame(results)
df.to_csv("distances.csv", index=False)

df2 = pd.DataFrame(pocket_storage)
df2.to_csv("pockets.csv", index=False)

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

plt.savefig("all_fe_out_of_plane.png", dpi=300)
if "test_data" in prefix:
    plt.show()

