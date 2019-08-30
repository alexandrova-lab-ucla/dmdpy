#!/usr/bin/env python3


import logging

import dmdpy.utilities.constants as constants
from dmdpy.protein.chain import Chain
from dmdpy.protein.residue import Residue


class Protein:

    __slots__ = ['chains', '_logger', 'non_residues', 'metals', 'name', 'sub_chain']

    def __init__(self, name: str, chains: [Chain]):

        #TODO: add a flag that is enabled if reformat is called for the protein to speed up algorithm of finding atoms/residues

        self._logger = logging.getLogger(__name__)

        self._logger.debug(f"Initializing vars for protein {name}")
        self.name = name
        self.chains = chains
        self.non_residues = []
        self.metals = []
        self.sub_chain = Chain()

        self._logger.debug(f"Created protein {str(self)}")

    def reformat_protein(self):
        res_num = 1
        chain_let = 'A'
        for (chain_num, chain) in enumerate(self.chains):
            atom_num = 1
            for (residue_num, residue) in enumerate(chain.residues):

                # Case where we have a non-residue/metal, most likely substrate
                if residue.name not in constants.AMINO_ACID_RESIDUES:
                    for (atom_num, atom) in enumerate(residue.atoms):
                        # Remove metal from this, and make it its own residue essentially
                        if atom.element in constants.METALS:
                            self._logger.debug(f"Found a metal: {atom} in residue {residue}")
                            self.metals.append(atom)
                            del residue[atom_num]

                    self._logger.debug(f"Removing residue: {residue}")
                    self.non_residues.append(residue)
                    del chain.residues[residue_num]

                else:
                    #update residue, and atomic numbering for normal atoms
                    self._logger.debug(f"Renumbering residue: {residue} to {res_num} and its atoms starting at {atom_num}")
                    residue.number = res_num
                    for atom in residue.atoms:
                        atom.number = atom_num
                        atom_num += 1

                    res_num += 1

            # Check if we now just have an empty chain
            if not chain.residues:
                self._logger.debug(f"Removing chain: {chain}")
                del self.chains[chain_num]

            else:
                chain.name = chain_let
                chain_let = chr(ord(chain_let) + 1)

        #Now we add in the non-residues and the metals into a new chain!
        res_num = 1
        if self.non_residues or self.metals:
            self._logger.debug("Creating a substrate/metal chain")
            atom_num = 1
            sub_chain = Chain(chr(ord(self.chains[-1].name) + 1))

            # Want to sort the metals!
            self.metals.sort(key=lambda metal: metal.element)

            cur_metal = ""
            metal_num = 1
            for metal in self.metals:
                metal.number = atom_num
                if cur_metal != metal.element:
                    metal_num = 1
                    cur_metal = metal.element

                if len(metal.element) == 1:
                    name = metal.element + f"{metal_num:02d}"

                elif len(metal.element) == 2:
                    name = metal.element + f"{metal_num:01d}"

                else:
                    self._logger.error(f"Encountered a metal with an unusual element ID: {metal}")
                    raise ValueError

                if len(name) != 3:
                    self._logger.error(f"The name for this is too long: {metal} with name: {name}")
                    raise ValueError

                metal_residue = Residue(name=name, number=res_num)
                self._logger.debug(f"Adding residue {metal_residue} to substrate chain")
                metal_residue.add_atom(metal)
                sub_chain.add_residue(metal_residue)

                res_num += 1
                metal_num += 1
                atom_num += 1

            for residue in self.non_residues:
                self._logger.debug(f"Adding residue {residue} to substrate chain")
                residue.number = res_num
                for atom in residue.atoms:
                    atom.number = atom_num
                    atom_num += 1

                sub_chain.add_residue(residue)
                res_num += 1

            self._logger.debug("Adding substrate chain to master chain")
            self.sub_chain = sub_chain

    def get_atom(self, identifier):
        for chain in self.chains:
            if chain.name == identifier[0]:
                for residue in chain.residues:
                    if residue.number == identifier[1]:
                        for atom in residue.atoms:
                            if atom.id == identifier[2]:
                                return atom

        self._logger.error(f"Could not find requested atom {identifier}")
        raise ValueError

    def get_residue(self, identifier):
        for chain in self.chains:
            if chain.name == identifier[0]:
                for residue in chain.residues:
                    if residue.number == identifier[1]:
                        return residue

        self._logger.error("Could not find requested residue")
        raise ValueError

    def write_pdb(self):
        self._logger.debug(f"Writing out pdb: {self}")
        with open(self.name, 'w') as pdb:
            for chain in self.chains:
                for residue in chain.residues:
                    for atom in residue.atoms:
                        pdb.write(atom.pdb_line())
            pdb.write('TER\n')
            for residue in self.sub_chain.residues:
                for atom in residue.atoms:
                    pdb.write(atom.pdb_line())
                pdb.write('TER\n')

            pdb.write('ENDMDL')
