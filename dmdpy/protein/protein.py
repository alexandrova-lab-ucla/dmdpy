#!/usr/bin/env python3


import logging
import pkg_resources
import csv

from dmdpy.utility import constants
from dmdpy.protein.chain import Chain
from dmdpy.protein.residue import Residue

__all__=[
    'Protein'
]

#TODO Add a function to remove Hydrogens that were deprotonated or protonated, ie standard protonation states

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
        # This is the BIG BIG BIG function that fixes EVERYTHING of a pdb for DMD
        res_renum = 1
        chain_let = 'A'

        chain_num = 0
        while chain_num < len(self.chains):
            atom_renum = 1
            residue_num = 0
            while residue_num < len(self.chains[chain_num].residues):
                # Case where we have a non-residue/metal, most likely substrate
                if self.chains[chain_num].residues[residue_num].name not in constants.AMINO_ACID_RESIDUES:
                    atom_num = 0
                    while atom_num < len(self.chains[chain_num].residues[residue_num].atoms):
                        # Remove metal from this, and make it its own residue essentially
                        if self.chains[chain_num].residues[residue_num].atoms[atom_num].element in constants.METALS:
                            self._logger.debug(f"Found a metal: {self.chains[chain_num].residues[residue_num].atoms[atom_num]} in residue {self.chains[chain_num].residues[residue_num]}")
                            self.metals.append(self.chains[chain_num].residues[residue_num].atoms.pop(atom_num))
                            atom_num -= 1

                        atom_num += 1

                    self._logger.debug(f"Removing residue: {self.chains[chain_num].residues[residue_num]}")
                    if self.chains[chain_num].residues[residue_num].atoms:
                        self.non_residues.append(self.chains[chain_num].residues.pop(residue_num))

                    else:
                        del self.chains[chain_num].residues[residue_num]

                else:
                    # update residue, and atomic numbering for normal atoms
                    self._logger.debug(f"Renumbering residue: {self.chains[chain_num].residues[residue_num]} to {res_renum} and its atoms starting at {atom_renum}")
                    self.chains[chain_num].residues[residue_num].number = res_renum
                    self.chains[chain_num].residues[residue_num].inConstr_number = residue_num
                    for atom in self.chains[chain_num].residues[residue_num].atoms:
                        atom.number = atom_renum
                        atom_renum += 1

                    res_renum += 1
                    residue_num += 1

            if not self.chains[chain_num].residues:
                self._logger.debug(f"Removing chain: {self.chains[chain_num]}")
                del self.chains[chain_num]

            else:
                self.chains[chain_num].name = chain_let
                chain_let = chr(ord(chain_let) + 1)
                chain_num += 1

        #Now we add in the non-residues and the metals into a new chain!
        res_num = 1
        if self.non_residues or self.metals:
            self._logger.debug("Creating a substrate/metal chain")
            atom_num = 1

            # in the event that we have only substrate
            if not len(self.chains):
                self.sub_chain.name = 'A'

            else:
                self.sub_chain.name = chr(ord(self.chains[-1].name) + 1)

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
                    name = metal.element.upper() + f"{metal_num:02d}"

                elif len(metal.element) == 2:
                    name = metal.element.upper() + f"{metal_num:01d}"

                else:
                    self._logger.error(f"Encountered a metal with an unusual element ID: {metal}")
                    raise ValueError

                if len(name) != 3:
                    self._logger.error(f"The name for this is too long: {metal} with name: {name}")
                    raise ValueError

                # DMD does not know how to handle any other metal but zinc
                if metal.element.lower() != "zn":
                    metal.element = 'Zn'

                metal_residue = Residue(name=name, number=res_num)
                self._logger.debug(f"Adding residue {metal_residue} to substrate chain")
                self.sub_chain.add_residue(metal_residue)
                metal_residue.add_atom(metal)

                res_num += 1
                metal_num += 1
                atom_num += 1

            # Issue with assigning chain to atom
            for residue in self.non_residues:
                start = 100
                for atom in residue.atoms:
                    atom.id = f"{atom.element.upper()}{start}"
                    start += 1

                self._logger.debug(f"Adding residue {residue} to substrate chain")
                self.sub_chain.add_residue(residue)
                residue.number = res_num
                residue.inConstr_number = res_num
                for atom in residue.atoms:
                    atom.number = atom_num
                    atom_num += 1

                res_num += 1

            self._logger.debug("Adding substrate chain to master chain")
            self.chains.append(self.sub_chain)

        self.relabel()

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

    def get_chain(self, identifier):
        for chain in self.chains:
            if chain.name == identifier:
                return chain

        self._logger.error("Could not find the requested chain")
        raise ValueError

    def write_pdb(self):
        self._logger.debug(f"Writing out pdb: {self}")
        try:
            with open(self.name, 'w') as pdb:
                for chain in self.chains[:-1]:
                    for residue in chain.residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                    pdb.write('TER\n')

                if not self.sub_chain.residues:
                    for residue in self.chains[-1].residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                    pdb.write('TER\n')

                else:
                    for residue in self.sub_chain.residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                        pdb.write('TER\n')

                pdb.write('ENDMDL\n')

        except IOError:
            self._logger.exception(f"Error writing out to file {self.name}")
            raise

    def relabel(self, format: str="DMD"):

        atom_label_dict = {}
        with open(pkg_resources.resource_filename('dmdpy.resources', 'atom_label.csv')) as csvfile:
            csvreader = csv.reader(csvfile)
            schemenames = next(csvreader)[1:]
            try:
                newid = schemenames.index(format)

            except ValueError:
                raise ValueError("Format key not found in atom_label.csv")

            for row in csvreader:
                if row[0] in atom_label_dict:
                    atom_label_dict[row[0]].append(row[1:])

                else:
                    atom_label_dict[row[0]] = [row[1:]]

        def rename_residue(residue, terminus: str = None):

            # Checks for any amino acids/molecules not in the csv file first!
            if residue.name not in atom_label_dict.keys():
                self._logger.warning(f"Residue: {residue.name}{residue.number} not in atom_label.csv!")
                return

            #Find the column that has the current naming scheme present
            for schemeid in range(len(schemenames)):
                scheme = []
                for namelist in atom_label_dict[residue.name]:
                    scheme.append(namelist[schemeid])

                if terminus is not None:
                    for namelist in atom_label_dict[terminus]:
                        scheme.append(namelist[schemeid])
                # Check to see if this is the naming scheme
                for atom in residue.atoms:
                    if atom.id not in scheme:
                        break

                else:
                    break # This is the correct naming scheme!

            else:
                raise ValueError(f"Could not find the naming scheme for {residue.name}{residue.number}")

            #Loop over all of the atoms
            for atom in residue.atoms:
                old_atomid = scheme.index(atom.id)
                if terminus is not None:
                    atom.id = (atom_label_dict[residue.name] + atom_label_dict[terminus])[old_atomid][newid]

                else:
                    atom.id = atom_label_dict[residue.name][old_atomid][newid]

            # TODO add Jacks glycine hydrogen fixed so that naming convention is always the same with Chimera

        for chain in self.chains:
            index = 0
            rename_residue(chain.residues[0], "NTERM")
            rename_residue(chain.residues[-1], "CTERM")
            cterm = chain.residues[-1]

            for residue in chain.residues[1:-1]:
                rename_residue(residue)


