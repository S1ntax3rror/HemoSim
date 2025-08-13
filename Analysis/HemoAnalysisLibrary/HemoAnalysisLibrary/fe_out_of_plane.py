import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


coords = ["../permutation0000/step5_999.dcd", "../permutation0000/step5_998.dcd", "../permutation0000/step5_997.dcd"]
u = mda.Universe("../permutation0000/run_setup.psf", coords)
heme_segids = sorted(set(atom.segid for atom in u.select_atoms("resname HEME")))

# Triangular combinations
triangles = [
    ('NA', 'NB'),
    ('NB', 'NC'),
    ('NC', 'ND'),
    ('ND', 'NA'),
]

results = []
for ts in u.trajectory:
    frame_data = {'frame': ts.frame}

    for segid in heme_segids:
        heme = u.select_atoms(f"resname HEME and segid {segid}")
        if len(heme) == 0:
            continue  # skip if this segid doesn't have a heme

        try:
            fe = heme.select_atoms("name FE")[0]
            pos_fe = fe.position

            # nitrogen atoms
            n_atoms = {name: heme.select_atoms(f"name {name}")[0].position for name in ['NA', 'NB', 'NC', 'ND']}
            nitrogen_positions = np.array(list(n_atoms.values()))    
            
            centroid = nitrogen_positions.mean(axis=0)       
            centered = nitrogen_positions - centroid
            _, _, vh = np.linalg.svd(centered) #to calculate unit normal vector, 2nd index
            normal = vh[2, :]  # normal to plane

            # Signed out-of-plane distance
            displacement = np.dot(pos_fe - centroid, normal)
            frame_data[f"{segid}_fe_out_of_plane"] = displacement

            # angles and distances
            for n1_name, n2_name in triangles:
                n1 = n_atoms[n1_name]
                n2 = n_atoms[n2_name]

                v1 = n1 - pos_fe
                v2 = n2 - pos_fe
                d1 = np.linalg.norm(v1)
                d2 = np.linalg.norm(v2)

                cos_theta = np.dot(v1, v2) / (d1 * d2) #angle in degree calculation
                angle = np.arccos(np.clip(cos_theta, -1.0, 1.0)) * 180 / np.pi

                prefix = f"{segid}_{n1_name.lower()}_fe_{n2_name.lower()}"
                frame_data[f"{prefix}_angle"] = angle
                frame_data[f"{prefix}_dist1"] = d1 #first one
                frame_data[f"{prefix}_dist2"] = d2 #second one

        except IndexError:
            print(f"Missing atoms in segid {segid} at frame {ts.frame}")
            continue

    results.append(frame_data)


df = pd.DataFrame(results)
df.to_csv("heme_fe_geometry_all4.csv", index=False)

for segid in heme_segids:
    col = f"{segid}_fe_out_of_plane" # change this if you want to plot something else such as distances, angles
    if col in df.columns:
        plt.plot(df['frame'], df[col], label=f'{segid}')

plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.title("Fe Out-of-Plane Displacement for All Heme Groups")
plt.xlabel("Frame") #can be changed to time
plt.ylabel("Displacement (Å)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("all_fe_out_of_plane.png", dpi=300)
plt.show()

