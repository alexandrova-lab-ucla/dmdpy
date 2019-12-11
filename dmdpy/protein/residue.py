#!/usr/bin/env python3

from dmdpy.protein.atom import Atom

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

    def write_inConstr(self):
        return f"{ord(self.chain.name) - ord('A') + 1}.{self.inConstr_number}.*"

    def __str__(self):
        return f"{self.name} {self.number}"
