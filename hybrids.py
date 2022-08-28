'''This Python script reads the geometrical structure of a molecule, as given in the xyz file, and computes the hybridization indexes and pyramid
angle of C atoms, if any'''

from collections import Counter
from collections import defaultdict
from ssl import HAS_TLSv1_1
import sys
import numpy as np
import math
import matplotlib.pyplot as plt

# Defining function read_xyz() that reads a xyz file (filename)
def read_xyz(filename):
    """This function reads atom labels and coordinates from a xyz file.
    
    Args:
        filename (string): the name of the xyz file
    
    Returns:
        atoms (list): a list of atom labels 
        coordinates (list): a list of [x,y,z] coordinates 
    
    Raises:
        ValueError: the function raises an error if it reads a number of atoms different from what 
            specified by the header of the xyz file.
    """

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
print('#'*80)
print("")
print("filename:                %s" % file)
print("number of atoms:         %d" % n_atoms)
print("number of distinct C:    %d" % n_carbons)
print("")
print('#'*80)

coord_np = np.array(coordinates)

def nn(i, coord):
    """ This function determines the first three nearest-neighbors of an atom.

    Args:
        i (integer): index of the target atom 
        coord (NumPy array): array specifying xyz coordinates of atoms in the molecule

    Returns:
        nn_indx (list): list of nearest-neighbors' indexes
        nn_dist (list): list of corresponding distances from target atom (i.e. bond distances)
    """

    # Finding number of atoms
    natom, ncoord = coord.shape

    # Computing Euclidean distances 
    dist_list = np.array([np.linalg.norm(coord[i]-coord[j]) for j in range(natom) ])
    nn_indx = list(np.argsort(dist_list)[1:4])
    nn_dist = dist_list[nn_indx]

    return nn_indx, nn_dist

nn_1, dist_nn1_ = nn(1,coord_np)

def angle(i, coord_np):
    '''Function that computes bond angles for an atom of index i'''
    
    # Getting first three-nearest neighbors calling function nn():
    nn_i, dist_nn_i = nn(i, coord_np)
    
    # Building vectors connecting atom i to its nearest neighbors and normalizing them:
    vector_nn = np.array([coord_np[i] - coord_np[nn_i[j]] for j in range(3)])
    unit_vector_nn = np.array([vector_nn[i]/np.linalg.norm(vector_nn[i]) for i in range(3)])
    
    # Computing dot products between units vetors:
    dot_prods = []
    prod1 = np.dot(unit_vector_nn[0], unit_vector_nn[1])
    prod2 = np.dot(unit_vector_nn[1], unit_vector_nn[2])
    prod3 = np.dot(unit_vector_nn[0], unit_vector_nn[2])
    dot_prods.extend([prod1, prod2, prod3])

    angles = [math.degrees(np.arccos(x)) for x in dot_prods]

    return list(dot_prods), angles

def hybrid(i, coord_np):
    '''This function computes hybrid orbitals for atom of index i'''
    p_hybrids = []
    s_hybrids = []

    # Getting angles:
    cosines, bond_angles = angle(i, coord_np)

    # Computing p and s weights
    p1 = - cosines[2]/(cosines[0]*cosines[1])
    s1 = 1/(1+p1)

    p2 = -cosines[1]/(cosines[0]*cosines[2])
    s2 = 1/(1+p2)

    p3 = -cosines[0]/(cosines[1]*cosines[2])
    s3 = 1/(1+p3)

    p_hybrids.extend([p1,p2,p3])

    # Computing p4
    sh = 0.0
    for j in p_hybrids:
        sh += (1+j)**(-1)
    p4 = 1/(1-sh)-1
    s4 = 1 - (s1+s2+s3)

    p_hybrids.extend([p4])
    s_hybrids.extend([s1, s2, s3, s4])

    return s_hybrids, p_hybrids

# Computing s-weights for any C in the lattice

print('')
print('')
atom_header = "Atom"
label_header = "Label"
hybrid_header = "s-weights"

print('{:<8}{:<8}{:<8}'.format(atom_header, label_header, hybrid_header))

perc_s_pi_like = []
for i in range(n_atoms):
    if atoms[i] == 'C':
        s_weights, p_weights = hybrid(i, coord_np)
        perc_s_pi_like.append(s_weights[3]*100)
        print("{:<8}{:<8}{}".format(atoms[i], i, s_weights))
    else:
        none_string = "Not a C! Not computed"
        print("{:<8}{:<8}{}".format(atoms[i], i, none_string))
