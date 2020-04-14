#!/usr/bin/env python

"""
Calculate RMSD between two XYZ files
by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Andersen Bratholm <larsbratholm@gmail.com>
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE
"""

import numpy as np
import re



def kabsch_rmsd(P, Q, outputCoords = False):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = rotate(P, Q)
    
    if outputCoords:
        return rmsd(P, Q), P
    else:
        return rmsd(P, Q)


def rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.
    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.
    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix
    Returns:
    U -- Rotation matrix
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def write_coordinates_xyz(atomInfo, V, argsOut):
    """
    Print coordinates V
    """
    N, D = V.shape

    print str(N)
    print

    with open(argsOut, 'w') as outputFile:
        for i in xrange(N):
            line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}\n".format(atomInfo[i], V[i, 0], V[i, 1], V[i, 2])
            outputFile.write(line)

def write_coordinates_pdb(atomInfo, V, argsOut):
    """
    Print coordinates V
    """
    N = len(atomInfo)

    
    with open(argsOut, 'a') as outputFile:
        for i in xrange(N):
            if atomInfo[i][0].startswith("TER") or atomInfo[i][0].startswith("END"):
                line = atomInfo[i][0] + "\n"
                outputFile.write(line)
            else:
                line = "{0:6s}{1:6s}{2:4s}{3:1s}{4:3s} {5:1s}{6:4s}{7:1s}   \
{8:8.3f}{9:8.3f}{10:8.3f}{11:22s}{12:2s}\n".format(atomInfo[i][0],  atomInfo[i][1], 
                        atomInfo[i][2],  atomInfo[i][3],  atomInfo[i][4], atomInfo[i][5],
                        atomInfo[i][6],  atomInfo[i][7],  V[i, 0], V[i, 1], V[i, 2],
                        atomInfo[i][11], atomInfo[i][12] )
        #line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atomInfo[i], V[i, 0], V[i, 1], V[i, 2])
                outputFile.write(line)

def write_coordinates(atoms, V, fmt, argsOut):
    """
    Print coordinates V
    """
    if fmt == "xyz":
        return write_coordinates_xyz(atoms, V, argsOut)
    elif fmt == "pdb":
        return write_coordinates_pdb(atoms, V, argsOut)
    exit("Could not recognize file format: {:s}".format(fmt))


def get_coordinates(filename, fmt, ignore_hydrogens):
    """
    Get coordinates from filename.
    """
    if fmt == "xyz":
        return get_coordinates_xyz(filename, ignore_hydrogens)
    elif fmt == "pdb":
        return get_coordinates_pdb(filename, ignore_hydrogens)
    exit("Could not recognize file format: {:s}".format(fmt))


def get_coordinates_pdb(filename, ignore_hydrogens):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.
    """

    V =[]
    skipThis = False

    with open(filename) as f:
        lines = f.readlines()
        myPDB = [ [] for _ in xrange(len(lines)) ]
        for idx, line in enumerate(lines):
            if line.startswith("TER") or line.startswith("END"):
                skipThis = True
                myPDB[idx].append(line.strip())
            elif line.startswith("ATOM"):
                myPDB[idx].append(line[0:6])
                myPDB[idx].append(line[6:12])
                myPDB[idx].append(line[12:16])
                myPDB[idx].append(line[16:17])
                myPDB[idx].append(line[17:20])
                myPDB[idx].append(line[21:22])
                myPDB[idx].append(line[22:26])
                myPDB[idx].append(line[26:27])
                myPDB[idx].append(line[30:38])
                myPDB[idx].append(line[38:46])
                myPDB[idx].append(line[46:54])
                myPDB[idx].append(line[54:76])
                myPDB[idx].append(line[76:78])

                if skipThis == True:
                    skipThis = False
                else:
                    V.append(np.asarray([myPDB[idx][8],myPDB[idx][9],myPDB[idx][10]],dtype=float))

                #atom = myPDB[idx][-1].strip()
                #if ignore_hydrogens and atom == "H":
                #    continue
                #elif atom in ["H", "C", "N", "O", "S", "P"]:
                #    atoms.append(atom)

    V = np.asarray(V)
    return myPDB, V


def get_coordinates_xyz(filename, ignore_hydrogens):
    """
    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.
    This function has been written to parse XYZ files, but can easily be
    written to parse others.
    """
    f = open(filename, 'r')
    V = []
    atoms = []
    n_atoms = 0
    lines_read = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.next())
    except ValueError:
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.next()

    # Use the number of atoms to not read beyond the end of a file
    for line in f:

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        # ignore hydrogens
        if ignore_hydrogens and atom.lower() == "h":
            continue

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

        lines_read += 1

    f.close()
    V = np.array(V)
    return atoms, V


if __name__ == "__main__":
    import argparse

    description = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ or PDB format.
The order of the atoms *must* be the same for both structures.
"""

    epilog = """
The script will return three RMSD values:
1) Normal: The RMSD calculated the straight-forward way.
2) Kabsch: The RMSD after the two coordinate sets are translated and rotated onto each other.
"""

    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('structure_a', metavar='structure_a.xyz', type=str)
    parser.add_argument('structure_b', metavar='structure_b.xyz', type=str)
    parser.add_argument('-o', '--output', action='store', help='print out structure A, centered and rotated unto structure B\'s coordinates in XYZ format')
    parser.add_argument('-n', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    parser.add_argument('-f', '--format', action='store', help='Format of input files. Supports xyz or pdb.', default="xyz")

    args = parser.parse_args()

    atomInfoP, P = get_coordinates(args.structure_a, args.format, args.no_hydrogen)
    atomInfoQ, Q = get_coordinates(args.structure_b, args.format, args.no_hydrogen)

    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if args.output:
        V = rotate(P, Q)
        V += Qc
        write_coordinates(atomInfoP, V, args.format, args.output)

    #print "Normal RMSD:", normal_rmsd
    print kabsch_rmsd(P, Q)
    
