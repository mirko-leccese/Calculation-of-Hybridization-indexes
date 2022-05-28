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

#print(nn(1, coord_np))

nn_1, dist_nn1_ = nn(1,coord_np)

def angle(i, coord_np):
    '''Function that computes bond angles for an atom of index i'''
    
    # Getting first three-nearest neighbors calling function nn():
    nn_i, dist_nn_i = nn(i, coord_np)
    
    # Building vectors connecting atom i to its nearest neighbors and normalizing them:
    vector_nn = np.array([coord_np[i] - coord_np[nn_i[j]] for j in range(3)])
    unit_vector_nn = np.array([vector_nn[i]/np.linalg.norm(vector_nn[i]) for i in range(3)])
    
    # Computing dot products between units vetors:
    dot_prods = set(
            [np.dot(unit_vector_nn[i], unit_vector_nn[j]) for i in range(3) for j in range(3) if i != j]
            )
    # Computing angles (in  degrees):
    angles = [math.degrees(np.arccos(x)) for x in dot_prods]

    return list(dot_prods), angles

def hybrid(i, coord_np):
    '''This function computes hybrid orbitals for atom of index i'''
    p_hybrids = []

    # Getting angles:
    cosines, bond_angles = angle(i, coord_np)
    p1 = - cosines[2]/(cosines[0]*cosines[1])
    p2 = -cosines[1]/(cosines[0]*cosines[2])
    p3 = -cosines[0]/(cosines[1]*cosines[2])
    p_hybrids.extend([p1,p2,p3])

    # Computing h4
    sh = 0.0
    for j in p_hybrids:
        sh += (1+j)**(-1)
    
    print(sh)
    p4 = 1/(1-sh)-1

    p_hybrids.extend([p4])

    return p_hybrids

phyb_atom1 = hybrid(1,coord_np)

print(phyb_atom1)

#print(atoms)

#for i in len(atoms):
#    if atoms[i] == 'C':

