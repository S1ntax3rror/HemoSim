import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import MDAnalysis as mda
from MDAnalysis.analysis import distances

cutoff = 3.2  #  # Angstrom cutoff based on VMD 

# Loop over each permutation directory
for perm_dir in sorted(glob.glob("permutation*")):
    if not os.path.isdir(perm_dir):
        continue
        
    pdb_or_psf = None
    traj_files = []
    for ext in ["pdb", "psf"]:
        files = glob.glob(os.path.join(perm_dir, f"*.{ext}"))
        if files:
            pdb_or_psf = files[0]
            break


    traj_files = sorted(glob.glob(os.path.join(perm_dir, "step5_*.dcd")))
    if not pdb_or_psf or not traj_files:
        print(f"  Missing files in {perm_dir}, skip.")
        continue

    u = mda.Universe(pdb_or_psf, traj_files)

    # Select acidic and basic residues
    acidic_residues = [res for res in u.select_atoms("resname ASP GLU").residues]
    basic_residues = [res for res in u.select_atoms("resname ARG LYS").residues]

    # distance calculations
    acidic_atoms = u.select_atoms("resname ASP GLU and name OD1 OD2 OE1 OE2")
    basic_atoms = u.select_atoms("resname ARG LYS and name NZ NH1 NH2 NE")

    bridge_ids = set()
    bridge_data = defaultdict(list)
    frame_numbers = []

    for ts in u.trajectory: # use [::20] to skip frames or save every 20 steps 
        frame_numbers.append(ts.frame)
        current_frame_bridges = {}

        # Distance matrix for this frame
        dists = distances.distance_array(basic_atoms.positions, acidic_atoms.positions)

        for a_res in acidic_residues:
            if a_res.resname == "ASP":
                acid_group = a_res.atoms.select_atoms("name OD1 OD2")
            elif a_res.resname == "GLU":
                acid_group = a_res.atoms.select_atoms("name OE1 OE2")
            else:
                continue
            com_acid = acid_group.center_of_mass()

            for b_res in basic_residues:
                if b_res.resname == "LYS":
                    base_group = b_res.atoms.select_atoms("name NZ")
                elif b_res.resname == "ARG":
                    base_group = b_res.atoms.select_atoms("name NH1 NH2 NE")
                else:
                    continue
                com_base = base_group.center_of_mass()

                # Check cutoff
                close = False
                for b_atom in base_group:
                    for a_atom in acid_group:
                        b_idx = basic_atoms.atoms.indices.tolist().index(b_atom.index)
                        a_idx = acidic_atoms.atoms.indices.tolist().index(a_atom.index)
                        if dists[b_idx, a_idx] <= cutoff:
                            close = True
                            break
                    if close:
                        break

                if close:
                    bridge = f"{b_res.resname}{b_res.resid}{b_res.segid}-{a_res.resname}{a_res.resid}{a_res.segid}"
                    bridge_ids.add(bridge)
                    current_frame_bridges[bridge] = round(np.linalg.norm(com_acid - com_base), 2)

        # Append data 
        for bridge in bridge_ids:
            bridge_data[bridge].append(current_frame_bridges.get(bridge, np.nan))

    # Pad missing data
    for bridge in bridge_data:
        if len(bridge_data[bridge]) < len(frame_numbers):
            bridge_data[bridge] += [np.nan] * (len(frame_numbers) - len(bridge_data[bridge]))

    # Save CSV for this permutation
    df = pd.DataFrame(bridge_data, index=frame_numbers)
    df.index.name = "Frame"
    out_csv = os.path.join(perm_dir, f"salt_bridge_COM_distances_{perm_dir}.csv")
    df.to_csv(out_csv)

    # plotting
    plt.figure(figsize=(12, 6))
    for bridge in df.columns:
        plt.plot(df.index, df[bridge], label=bridge)
    plt.title(f"Salt Bridge COM Distances Over Time ({perm_dir})")
    plt.xlabel("Frame #")
    plt.ylabel("Distance (Å)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="small")
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(os.path.join(perm_dir, f"salt_bridge_COM_distances_{perm_dir}.png"), dpi=300)
    plt.close()

