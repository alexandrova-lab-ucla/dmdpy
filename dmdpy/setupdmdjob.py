#!/usr/bin/env python3

import os
import logging
from subprocess import Popen, PIPE


import dmdpy.utilities.utilities as utilities

logger = logging.getLogger(__name__)


class setupDMDjob:

    def __init__(self, parameters: dict=None, dir: str="./"):
        logger.debug("Entered setupDMDjob ")

        # Instance variables
        logger.debug("Initializing Variables")
        self._raw_parameters = parameters
        self._run_directory = dir
        self._initial_directory = os.getcwd()
        self._protein = utilities.load_pdb("test.pdb")

    def make_topparam(self):
        with open('topparam', 'w') as topparam_file:
            for residue in self._protein.sub_chain.residues:
                with open(f"{residue.name}.pdb", 'w') as mol2_file:
                    for atom in residue.atoms:
                        mol2_file.write(atom.pdb_line())

                    mol2_file.write('TER\nENDMDL')
                # Now we execute the babel command here
                with Popen(['babel', f"{residue.name}.pdb", f"{residue.name}.mol2"], stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
                    while shell.poll() is None:
                        logger.debug(shell.stdout.readline().strip())
                        logger.debug(shell.stderr.readline().strip())

                # Add the name to the topparm_file!
                topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2")



