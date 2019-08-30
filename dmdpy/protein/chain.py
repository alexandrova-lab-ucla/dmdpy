#!/usr/bin/env python3

from dmdpy.protein.residue import Residue


class Chain:

    __slots__ = ['name', 'residues']

    def __init__(self, name: str=''):
        self.name = name
        self.residues = []

    def add_residue(self, res: Residue):
        res.chain = self
        self.residues.append(res)

    def __str__(self):
        return self.name
