import importlib
import subprocess
import sys
import os


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
if "cluster" in cwd or "scratch" in cwd or "Hemo" in cwd:
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

# heme_segids = sorted(set(atom.segid for atom in u.select_atoms("resname HEME")))
CO_segids = ["COA", "COB", "COC", "COD"]
heme_segids = ['HETA', 'HETB', 'HETC', 'HETD']

triangles = [
    ('NA', 'NB'),
    ('NB', 'NC'),
    ('NC', 'ND'),
    ('ND', 'NA'),
]



def get_CO_positions(CO_seg_atom_list):
    if isinstance(CO_seg_atom_list, list):
        return [atom[1].position for atom in CO_seg_atom_list]
    atom = u.select_atoms(f"segid {CO_seg_atom_list}")
    return [atom[0].position, atom[1].position]


results = []
pocket_storage = []
for ts in u.trajectory:
    frame_data = {'frame': ts.frame}
    frame_pocket_data = {'frame': ts.frame}
    print(f"reading for frame {ts}")
    for CO_segid, heme_segid in zip(CO_segids, heme_segids):
        heme = u.select_atoms(f"resname HEME and segid {heme_segid}")

        try:
            fe = heme.select_atoms(f"name FE")
            pos_fe = fe[0].position

            # CO atoms
            CO_positions = get_CO_positions(CO_segid)
            CO_position = np.mean(CO_positions, axis=0)

            dist = np.linalg.norm(pos_fe-CO_position)
            if dist < 6:
                frame_pocket_data[CO_segid + " Pocket"] = "A or B"
            elif 6 <= dist <= 12.5:
                frame_pocket_data[CO_segid + " Pocket"] = "XE"
            else:
                frame_pocket_data[CO_segid + " Pocket"] = "NA"

            if dist < 24:
                frame_data[CO_segid + "_dist"] = dist

        except IndexError:
            print(f"Missing atoms in segid {heme_segid} at frame {ts.frame}")
            continue

        n_atoms = {name: heme.select_atoms(f"name {name}")[0].position for name in ['NA', 'NB', 'NC', 'ND']}

        nitrogen_positions = np.array(list(n_atoms.values()))

        centroid = nitrogen_positions.mean(axis=0)
        centered = nitrogen_positions - centroid
        _, _, vh = np.linalg.svd(centered)  # to calculate unit normal vector, 2nd index
        normal = vh[2, :]  # normal to plane

        # Signed out-of-plane distance
        displacement = np.dot(CO_position - centroid, normal)
        frame_data[f"{heme_segid}_fe_out_of_plane"] = displacement

        co_vector = CO_position - centroid
        co_unit = co_vector / np.linalg.norm(co_vector)
        normal_unit = normal / np.linalg.norm(normal)

        # Angle in radians
        angle_rad = np.arccos(np.clip(np.dot(co_unit, normal_unit), -1.0, 1.0))
        # Convert to degrees
        angle_deg = np.degrees(angle_rad)

        frame_data[f"{heme_segid}_CO_plane_angle"] = angle_deg


        # n1 = n_atoms[n1_name]
        # n2 = n_atoms[n2_name]
        #
        # v1 = n1 - CO_position
        # v2 = n2 - CO_position
        # d1 = np.linalg.norm(v1)
        # d2 = np.linalg.norm(v2)
        #
        # cos_theta = np.dot(v1, v2) / (d1 * d2)  # angle in degree calculation
        # angle = np.arccos(np.clip(cos_theta, -1.0, 1.0)) * 180 / np.pi
        #
        # prefix = f"{heme_segid}_{n1_name.lower()}_fe_{n2_name.lower()}"
        # frame_data[f"{prefix}_angle"] = angle
        # frame_data[f"{prefix}_dist1"] = d1  # first one
        # frame_data[f"{prefix}_dist2"] = d2  # second one


    results.append(frame_data)
    pocket_storage.append(frame_pocket_data)

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

