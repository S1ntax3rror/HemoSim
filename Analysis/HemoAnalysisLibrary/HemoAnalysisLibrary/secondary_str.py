import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
import matplotlib.pyplot as plt
import numpy as np


# Load PDB and PSF files
u = mda.Universe('step5_999.pdb', 'step5_999.dcd')

counts_per_frame = []
# Analyze secondary structure using DSSP for each frame in the trajectory
for ts in u.trajectory[::50]:  # Skipping frames 
    dssp = ''.join(DSSP(u, ts).run().results.dssp[1]) 
    print(dssp)    
    counts = {'H': 0, 'E': 0, 'C': 0}
    for ss in dssp:
        if ss == 'H':
            counts['H'] += 1
        elif ss == 'E':
            counts['E'] += 1
        elif ss == 'C':
            counts['C'] += 1

    counts_per_frame.append(counts)

counts_array = np.array([[counts['H'], counts['E'], counts['C']] for counts in counts_per_frame])
frames = np.arange(len(counts_per_frame))  # X-axis: frames
print(frames)
width = 0.2  

# figure 
fig, ax = plt.subplots()

# Plot bars 
bar1 = ax.bar(frames, counts_array[:, 0], width, label='Helix (H)', color='r')
#bar2 = ax.bar(frames, counts_array[:, 1], width, label='Sheet (E)', color='b')
#bar3 = ax.bar(frames, counts_array[:, 2], width, label='Coil (C)', color='g') #+width

ax.set_xlabel('Frame Index')
ax.set_ylabel('Count of Residues')
ax.set_title('Secondary Structure Counts Over Time')
ax.legend()

plt.show()

