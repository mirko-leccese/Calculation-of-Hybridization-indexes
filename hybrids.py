'''This Python script reads the geometrical structure of a molecule, as given in the xyz file, and computes the hybridization indexes and pyramid
angle of C atoms, if any'''

from collections import Counter
from collections import defaultdict
import sys
import numpy as np

# Defining function read_xyz() that reads a xyz file (filename)
def read_xyz(filename):
    '''This function return a list of atoms and coordinates
    
    If the number of coordinates does not agree with the number
    of atoms within the xyz file, it raises a ValueError'''

    atoms = []
    coordinates = []

    xyz_file = open(filename)
    
    # Reading number of atoms 
    n_atoms = int(xyz_file.readline())

    # Reading title
    title = xyz_file.readline()

    # Reading atom labels and xyz coordinates 
    for line in xyz_file:
        atom, x, y, z = line.split()
        atoms.append(atom)
        coordinates.append([float(x),float(y), float(z)])

    xyz_file.close()

    # Check if number of atoms and number of coordinates match 
    if n_atoms != len(coordinates):
        raise ValueError("File says %d atoms but read %d points." % (n_atoms, len(coordinates)))

    return atoms, coordinates

# Reading filename from stdin
file = input("xyz? ")

atoms, coordinates = read_xyz(file)

# Molecule info:
n_atoms = len(atoms)

atom_counter = Counter(atoms)
n_carbons = atom_counter.get("C", 0)
if n_carbons == 0:
    sys.exit("ERROR! No C atoms in the xyz file!")

print("")
print("############################################")
print("")
print("filename:                %s" % file)
print("number of atoms:         %d" % n_atoms)
print("number of distinct C:    %d" % n_carbons)
print("")
print("############################################")

coord_np = np.array(coordinates)

#dist = np.linalg.norm(coord_np[0]-coord_np[1])

dist_ind_dict = defaultdict(list)
dist_val_dict = defaultdict(list)

for i in range(n_atoms):
    dist_list = [np.linalg.norm(coord_np[i]-coord_np[j]) for j in range(n_atoms)]
    dist_np = np.array(dist_list)
    indx = list(np.argsort(dist_np)[1:4])

    dist_ind_dict[i]=indx
    dist_val_dict[i]=dist_np[indx]


for key, val in dist_ind_dict.items():
    print("Index of nn of", key, " are: ", val)
