#!/usr/bin/env python3

import logging

from dmdpy.protein.atom import Atom

class Residue:

    __slots__ = ["number", "name", "chain", "atoms", "_logger"]

    def __init__(self, line: str=None, name: str=None, number: int=None):
        self._logger = logging.getLogger(__name__)

        if line is None:
            self.name = name
            self.number = number

        elif line is not None:
            self.name = line[17:20]
            self.number = int(line[22:26])

        self.atoms = []
        self.chain = None
        self._logger.debug(f"Made residue: {str(self)}")

    def add_atom(self, atom: Atom):
        self._logger.debug(f"Adding atom {str(atom)} to residue: {str(self)}")
        atom.set_residue(self)
        atom.set_chain(self.chain)
        self.atoms.append(atom)

    def set_chain(self, chain):
        self._logger.debug(f"Setting chain for residue: {str(self)} to {str(chain)}")
        self.chain = chain

    def __str__(self):
        return f"{self.name} {self.number}"
