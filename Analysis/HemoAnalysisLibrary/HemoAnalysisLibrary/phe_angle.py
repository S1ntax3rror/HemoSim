import os
import glob
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral

# Loop over all permutation directories
for perm_dir in sorted(glob.glob("permutation*")):
    if not os.path.isdir(perm_dir):
        continue  
    dcd_files = sorted(glob.glob(os.path.join(perm_dir, "step5_*.dcd")))  
    u = mda.Universe(
        os.path.join(perm_dir, "run_setup.psf"),
        dcd_files
    )

    # phe98 are in chain a and c (not phe97)
    for chain, segid in [("A", "PROA"), ("C", "PROC")]:
        sel = u.select_atoms(f"segid {segid} and resname PHE and resid 98 and name N CA CB CG")
        dih = Dihedral([sel]).run()
        out_file = os.path.join(perm_dir, f"chi1_Phe98_chain{chain}_{perm_dir}.dat") # saved in the same directory as the permutation with the filename
        np.savetxt(out_file, dih.angles)

