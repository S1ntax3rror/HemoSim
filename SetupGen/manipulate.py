import os
import sys

from MyoSim.SimulationSetupAndRun.getPocketCenterOfPDB import read_pocket_list


def get_hem5_identifier_and_charges():
    with open(cwd + "/toppar/toppar_all36_prot_heme.str") as f:
        topparlines = f.readlines()
        f.close()

    HEM5_lines = []
    heme_started = False
    for line in topparlines:  # read all atoms in HEM5
        if heme_started and "GROUP" not in line:
            if "BOND" in line:
                break

            split_line = line.split(" ")
            while '' in split_line:
                split_line.remove('')

            split_line = split_line[1:4]

            if "\n" in split_line[-1]:
                split_line[-1] = split_line[-1].strip("\n")

            HEM5_lines.append(split_line)
        elif "RESI HEM5" in line:
            heme_started = True

    return HEM5_lines


# Check input
if len(sys.argv) < 2:
    raise SyntaxError("Input to small: ", sys.argv)

# Get case label
case = sys.argv[1]
cwd = os.getcwd()

# For case 'heme':
# Replace residue labels of the particular segment with the defined residue
if case.strip().lower() == 'heme':

    # Check input
    if len(sys.argv) < 4:
        raise SyntaxError("Input to small for HEME manipulation: ", sys.argv)

    segment = sys.argv[2]
    residue = sys.argv[3]

    # Check psf file
    psffile = os.path.join(cwd, "init_setup.psf")
    if not os.path.exists(psffile):
        print(psffile)
        raise SyntaxError(f"Missing psf file: {psffile}")

    # Check crd file
    crdfile = os.path.join(cwd, "init_setup.crd")
    if not os.path.exists(psffile):
        raise SyntaxError(f"Missing psf file: {crdfile}")

    Hem5data = get_hem5_identifier_and_charges()

    # Read psf file
    with open(psffile, 'r') as fpsf:
        psflines = fpsf.readlines()
        fpsf.close()
    newlines = psflines.copy()
    cnt = 0
    # Get lines to replace
    psf_atoms = []
    for il, line in enumerate(psflines):
        if segment in line and "HEME" in line:
            psf_atoms.append([line, il])

    # Replace lines
    for psf_atom_i, hem5_atom in zip(psf_atoms, Hem5data):
        i = psf_atom_i[1]
        psf_atom = psf_atom_i[0]
        replacer = psf_atom[:29] + f"{residue:4s}" + psf_atom[33:54] + f"{float(hem5_atom[-1]):10f}" + "       " + psf_atom[71:]
        newlines[i] = replacer

    # Collect sodium
    sod_list = []
    num_dumdum = 1
    for i, line in enumerate(psflines):
        if "SOD      SOD      SOD" in line:
            sod_list.append([line, i])
        if "DUM      DUM      DUM" in line:
            num_dumdum += 1

    # remove one charge by making one sodium a ghost
    dumdum_i = sod_list[-1]
    i = dumdum_i[1]
    dumdum = dumdum_i[0]
    replacer = dumdum[:11] + "DUM" + dumdum[14:20] + f"{str(num_dumdum):2s}" + dumdum[22:29] + f"{'DUM':3s}" + "      DUM      DUM    " + f"{float(0.00000):11f}" + "      " + dumdum[71:]
    newlines[i] = replacer

    dum_id = replacer[len("     "):len("     67954")]
    # print(num_dumdum)
    # print(dumdum)
    # print(replacer)

    # Write manipulated psf file
    with open(psffile, 'w') as fpsf:
        fpsf.writelines(newlines)

    # Read crd file
    with open(crdfile, 'r') as fcrd:
        crdlines = fcrd.readlines()
    newlines = crdlines.copy()

    # Assign segment the defined residue
    for il, line in enumerate(crdlines):
        if segment in line and not "NTHETA" in line:
            newlines[il] = line[:22] + f"{residue:4s}" + line[26:]
        if dum_id in line[:len("     67926    ")]:
            replacer = line[:22] + "DUM       DUM" + line[22 + len("DUM       DUM"):102] + "DUM" + line[105:]
            newlines[il] = replacer

    # Write manipulated crd file
    with open(crdfile, 'w') as fcrd:
        fcrd.writelines(newlines)

