#!/usr/bin/env python3

import numpy as np

from dmdpy.utility import constants

__all__=[
    'Atom'
]

class Atom:

    __slots__ = ['element', 'coords', 'id', 'residue', 'chain', 'number', 'bonds']

    def __init__(self, line: str = None, element: str = None, coords: np.array = None, id=None, number=None):

        if line is None:
            self.element = element
            self.coords = coords
            self.id = id
            self.number = number

        elif line is not None:
            self.element = line[76:78].strip().lower()
            self.coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            self.id = line[12:16].strip().upper()
            self.number = int(line[6:11])

        if self.element.lower() == 'eh':
            self.element = 'h'

        self.residue = None
        self.chain = None
        self.bonds = []


    def write_inConstr(self):
        if self.residue.name in constants.AMINO_ACID_RESIDUES:
            return f"{ord(self.chain.name)-ord('A')+1}.{self.residue.inConstr_number}.{self.id.upper()}"

        else:
            return f"{ord(self.chain.name) - ord('A') + self.residue.inConstr_number}.1.{self.id.upper()}"

    def pdb_line(self):
        return '{:<6}{:>5} {:<4} {} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00          {:>2}\n'.format(
            'ATOM' if self.residue.name in constants.AMINO_ACID_RESIDUES else "HETATM",
            self.number, self.id if len(self.id) > 3 else f" {self.id}", self.residue.name, self.residue.chain.name, self.residue.number,
            self.coords[0], self.coords[1], self.coords[2], self.element.capitalize())

    def __str__(self):
        return f"{self.id} {self.coords} {self.element}"

    def add_bond(self, atom):
        self.bonds.append(atom)
        atom.bonds.append(self)
