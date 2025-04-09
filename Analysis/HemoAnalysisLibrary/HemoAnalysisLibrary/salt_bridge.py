import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# Load file
u = mda.Universe("step5_999.pdb", "step5_999.dcd") #XX.pdb/psf, XX.dcd

# Get lists of acidic and basic ress
acidic_residues = [res for res in u.select_atoms("resname ASP GLU").residues]
basic_residues = [res for res in u.select_atoms("resname ARG LYS").residues]

# atom based list for the cut-off
acidic_atoms = u.select_atoms("resname ASP GLU and name OD1 OD2 OE1 OE2")
basic_atoms = u.select_atoms("resname ARG LYS and name NZ NH1 NH2 NE")

cutoff = 3.2  # Angstrom cutoff based on VMD 
bridge_ids = set()
bridge_data = defaultdict(list)
frame_numbers = []
com_all = []
for ts in u.trajectory: # use [::20] to skip frames or save every 20 steps 
    frame_numbers.append(ts.frame)
    current_frame_bridges = dict()

    # Distance matrix 
    dists = distances.distance_array(basic_atoms.positions, acidic_atoms.positions)

    for a_res in acidic_residues:
        if a_res.resname == "ASP":
            acid_group = a_res.atoms.select_atoms("name OD1 OD2")
        elif a_res.resname == "GLU":
            acid_group = a_res.atoms.select_atoms("name OE1 OE2")
        else:
            continue
        com_acid = acid_group.center_of_mass() # calculate center of mass

        for b_res in basic_residues:
            if b_res.resname == "LYS":
                base_group = b_res.atoms.select_atoms("name NZ")
            elif b_res.resname == "ARG":
                base_group = b_res.atoms.select_atoms("name NH1 NH2 NE")
            else:
                continue
            com_base = base_group.center_of_mass()

            # Check if any atom–atom pair between these groups is within cutoff
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
                com_all.append(round(np.linalg.norm(com_acid - com_base), 2))
    # Fill distances or NaN / can change
    for bridge in bridge_ids:
        bridge_data[bridge].append(current_frame_bridges.get(bridge, np.nan ))

# Padding
for bridge in bridge_data:
    if len(bridge_data[bridge]) < len(frame_numbers):
        bridge_data[bridge] += [np.nan] * (len(frame_numbers) - len(bridge_data[bridge]))

# Save to DataFrame
df = pd.DataFrame(bridge_data, index=frame_numbers)
df.index.name = "Frame"
df.to_csv("salt_bridge_COM_distances.csv")

# Plotting
plt.figure(figsize=(12, 6))
for bridge in df.columns:
    plt.plot(df.index, df[bridge], label=bridge)

plt.title("Salt Bridge COM Distances Over Time")
plt.xlabel("Frame #")
plt.ylabel("Distance (Å)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="small")
plt.tight_layout()
plt.grid(True)
plt.show()

