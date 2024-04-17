import sys
import numpy as np
from pymatgen.core.periodic_table import Element

l_units = 1.889725992
recip_e_units = 1.0 / 0.036749309
element_dict = {element.Z: element.symbol for element in Element}

args = sys.argv
if len(args) < 2:
    print("usage: cube2chg.py <.cube file>")
    sys.exit(1)

inputfilename = args[1]

print("Reading .cube file:", inputfilename)

with open(inputfilename, 'r') as f:
    cube_data = f.readlines()

total_atoms = int(cube_data[2].split()[0])
origin = np.array([float(datum) for datum in cube_data[2].split()[1:]])

nx = int(cube_data[3].split()[0])
basis_x = nx * np.array([float(datum) for datum in cube_data[3].split()[1:]])
ny = int(cube_data[4].split()[0])
basis_y = ny * np.array([float(datum) for datum in cube_data[4].split()[1:]])
nz = int(cube_data[5].split()[0])
basis_z = nz * np.array([float(datum) for datum in cube_data[5].split()[1:]])
basis = np.vstack((basis_x, basis_y, basis_z))

vol = np.linalg.det(basis)
vol = abs(vol) / (l_units ** 3)

atom_types = []
Zs = []
atom_counts = {}
coordinates = np.zeros((3, total_atoms))
density = []

line_index = 6
for i in range(total_atoms):
    line = cube_data[line_index + i].split()

    Z = int(line[0])
    if Z in Zs:
        atom_counts[Z] += 1
    else:
        Zs.append(Z)
        atom_counts[Z] = 1

    coordinates[:, i] = np.array(list(map(float, line[2:])))
coordinates = np.linalg.inv(basis.T) @ coordinates
elements = [element_dict[Z] for Z in Zs]

outputfilename = args[2] if len(args) >= 3 else inputfilename + ".chg"

with open(outputfilename, 'w') as f:
    f.write("VolumetricData\n")
    f.write(f"{1.0:19.14f}\n")
    f.write(f"{basis[0,0]:12.6f} {basis[0,1]:12.6f} {basis[0,2]:12.6f}\n")
    f.write(f"{basis[1,0]:12.6f} {basis[1,1]:12.6f} {basis[1,2]:12.6f}\n")
    f.write(f"{basis[2,0]:12.6f} {basis[2,1]:12.6f} {basis[2,2]:12.6f}\n")

    for element in elements:
        f.write(f"{element:<5} ")
    f.write("\n")

    for Z in Zs:
        f.write(f"{atom_counts[Z]:>6} ")
    f.write("\n")

    f.write("Direct\n")

    j = 0
    for Z in Zs:
        for i in range(atom_counts[Z]):
            f.write(f"{coordinates[0][j + i]:10.6f} {coordinates[1][j + i]:10.6f} {coordinates[2][j + i]:10.6f}\n")
        j += atom_counts[Z]

print("POSCAR file conversion complete.")
