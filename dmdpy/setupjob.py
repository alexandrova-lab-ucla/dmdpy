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
        self._dmd_config = utilities.load_dmd_config()

        self._protein = utilities.load_pdb("test.pdb")
        self.make_topparam()
        self.make_inConstr()

    def make_inConstr(self):
        # Need to also have the user constraints placed into here at some point, but that can wait for the time being!
        with open('inConstr', 'a') as inConstr_file:
            with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                       stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                       env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stderr.readline().strip())

        logger.debug("Finished making the inConstr file!")

    def make_topparam(self):
        with open('topparam', 'w') as topparam_file:
            for residue in self._protein.sub_chain.residues:
                successful = False
                with open(f"{residue.name}.pdb", 'w') as mol2_file:
                    for atom in residue.atoms:
                        mol2_file.write(atom.pdb_line())

                    mol2_file.write('TER\nENDMDL')
                # Now we execute the babel command here
                with Popen(f"babel {residue.name}.pdb {residue.name}.mol2", stdin=PIPE, stdout=PIPE, stderr=PIPE,
                           universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
                    while shell.poll() is None:
                        logger.debug(shell.stdout.readline().strip())
                        output = shell.stderr.readline().strip()
                        logger.debug(output)
                        if "1 molecule converted" in output:
                            successful = True
                
                if not successful:
                    logger.error("Could not create {residue.name} mol2 file!")
                    raise Exception

                logger.info(f"Successfuly made: {residue.name} mol2")
                # Add the name to the topparm_file!
                topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2")



