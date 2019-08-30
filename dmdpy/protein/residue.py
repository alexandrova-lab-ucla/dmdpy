#!/usr/bin/env python3

from dmdpy.protein.atom import Atom


class Residue:

    __slots__ = ["number", "name", "chain", "atoms"]

    def __init__(self, line: str=None, name: str=None, number: int=None):
        if line is None:
            self.name = name
            self.number = number

        elif line is not None:
            self.name = line[17:20]
            self.number = int(line[22:26])

        self.atoms = []
        self.chain = None

    def add_atom(self, atom: Atom):
        atom.residue = self
        atom.chain = self.chain
        self.atoms.append(atom)

    def __str__(self):
        return f"{self.name} {self.number}"
