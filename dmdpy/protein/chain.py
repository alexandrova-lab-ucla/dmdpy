#!/usr/bin/env python3

import logging

from dmdpy.protein.residue import Residue

class Chain:

    __slots__ = ['name', 'residues', '_logger']

    def __init__(self, name: str=''):
        self._logger = logging.getLogger(__name__)

        self.name = name
        self.residues = []

        self._logger.debug(f"Made chain {str(self)}")

    def add_residue(self, res: Residue):
        res.set_chain(self)
        self._logger.debug(f"Adding residue {str(res)} to chain {str(self)}")
        self.residues.append(res)

    def __str__(self):
        return self.name
