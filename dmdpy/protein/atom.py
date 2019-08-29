#!/usr/bin/env python3

import logging
import numpy as np

import dmdpy.utilities.constants as constants


class Atom:

    __slots__ = ['element', 'coords', 'id', 'residue', '_logger', 'chain', 'molecule', 'number']

    def __init__(self, line: str = None, element: str = None, coords: np.array = None, id=None):

        self._logger = logging.getLogger(__name__)

        if line is None:
            self._logger.debug("Adding elements individually")
            self.element = element
            self.coords = coords
            self.id = id

        elif line is not None:
            self._logger.debug("Adding elements line")
            self.element = line[76:78].strip().lower()
            self.coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            self.id = line[12:16]

        if self.element.lower() == 'eh':
            self.element = 'h'

        self._logger.debug("Initializing rest of vars")
        self.residue = None
        self.chain = None
        self.molecule = None
        self.number = None

        self._logger.debug(f"Made atom: {str(self)}")

    def set_residue(self, res):
        self._logger.debug(f"Setting residue for atom: {str(self)} to {str(res)}")
        self.residue = res

    def set_chain(self, chain):
        self._logger.debug(f"Setting chain for atom: {str(self)} to {str(chain)}")
        self.chain = chain

    def set_number(self, num: int):
        self._logger.debug(f"Setting number for atom: {str(self)} to {str(num)}")
        self.number = num

    def pdb_line(self):
        return '{:<6}{:>5} {:>4} {} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00          {:>2}\n'.format(
        'ATOM' if self.residue.name in constants.AMINO_ACID_RESIDUES else "HETATM",
        self.number, self.id, self.residue.name, self.chain.name, self.residue.number,
        self.coords[0], self.coords[1], self.coords[2], self.element.capitalize())

    def __str__(self):
        return f"{self.id} {self.coords} {self.element}"


