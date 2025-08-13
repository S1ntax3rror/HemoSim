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

# heme_segids = sorted(set(atom.segid for atom in u.select_atoms("resname HEME")))
CO_segids = ["COA", "COB", "COC", "COD"]
heme_segids = ['HETA', 'HETB', 'HETC', 'HETD']


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

