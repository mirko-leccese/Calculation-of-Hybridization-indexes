'''This Python script reads the geometrical structure of a molecule, as given in the xyz file, and computes the hybridization indexes and pyramid
angle of C atoms, if any'''

from collections import Counter
from collections import defaultdict
from ssl import HAS_TLSv1_1
import sys
import numpy as np
import math

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

# Initiliaze empty dictionary for nn indexes and nn distances:
dist_ind_dict = defaultdict(list)
dist_val_dict = defaultdict(list)

def nn(i, coord):
    ''' This function determines the first three nearest-neighbors of an atom of index i 
    in the set of atoms specified by the coordinates coord (a Numpy array) '''

    # Finding number of atoms
    natom, ncoord = coord.shape

    # Computing Euclidean distances 
    dist_list = np.array([np.linalg.norm(coord[i]-coord[j]) for j in range(natom) ])

    nn_indx = list(np.argsort(dist_list)[1:4])

    return nn_indx, dist_list[nn_indx]

print(nn(1, coord_np))

nn_1, dist_nn1_ = nn(1,coord_np)
# Computing angle 
vector_1 = coord_np[1] - coord_np[nn_1[0]]
vector_2 = coord_np[1] - coord_np[nn_1[1]]
vector_3 = coord_np[1] - coord_np[nn_1[2]]

unit_vector_1  = vector_1 / np.linalg.norm(vector_1)
unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
unit_vector_3 = vector_3 / np.linalg.norm(vector_3)


dot_product_12 = np.dot(unit_vector_1, unit_vector_2)
angle_12 = math.degrees(np.arccos(dot_product_12))

dot_product_13 = np.dot(unit_vector_1, unit_vector_3)
angle_13 = math.degrees(np.arccos(dot_product_13))

dot_product_23 = np.dot(unit_vector_2,unit_vector_3)
angle_23 = math.degrees(np.arccos(dot_product_23))


print(angle_12, angle_13, angle_23)

h1 = - dot_product_23 / (dot_product_12*dot_product_13)
h2 =  - dot_product_13 / (dot_product_23 * dot_product_12)
h3 = - dot_product_12 / (dot_product_23*dot_product_13)
print(h1, h2, h3)
