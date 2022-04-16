'''This Python script reads the geometrical structure of a molecule, as given in the xyz file, and computes the hybridization indexes and pyramid
angle of C atoms, if any'''

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
    title = xyz_file.readline()
    for line in xyz_file:
        atom, x, y, z = line.split()
        atoms.append(atom)
        coordinates.append([float(x),float(y), float(z)])
    xyz_file.close()

    if n_atoms != len(coordinates):
        raise ValueError("File says %d atoms but read %d points." % (n_atoms, len(coordinates)))

    return atoms, coordinates

file = input("xyz? ")

atoms, coordinates = read_xyz(file)



