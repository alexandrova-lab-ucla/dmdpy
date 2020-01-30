#!/usr/bin/env python3

from dmdpy.protein.atom import Atom
from dmdpy.utility import constants

#TODO have a heavy atom checker class to call a function that fixes
# have it pipe out to a file first and then have setupjob.py call a protein function to call chimera swapaa function

__all__=[
    'Residue'
]

class Residue:

    __slots__ = ["number", "name", "chain", "atoms", "inConstr_number"]

    def __init__(self, line: str=None, name: str=None, number: int=None):
        if line is None:
            self.name = name
            self.number = number

        elif line is not None:
            self.name = line[17:20]
            self.number = int(line[22:26])

        self.atoms = []
        self.chain = None
        self.inConstr_number = self.number

    def add_atom(self, atom: Atom):
        atom.residue = self
        atom.chain = self.chain
        self.atoms.append(atom)

    def set_chain(self, chain):
        pass

    def get_atom(self, name):
        for atom in self.atoms:
            if atom.id == name:
                return atom

        raise ValueError(f"Could not find atom: {name}")

    def write_inConstr(self):
        if self.name in constants.AMINO_ACID_RESIDUES:
            return f"{ord(self.chain.name) - ord('A') + 1}.{self.inConstr_number}.*"

        else:
            return f"{ord(self.chain.name)- ord('A') + self.inConstr_number}.1.*"

    def __str__(self):
        return f"{self.name} {self.number}"
