#!/usr/bin/env python3

import os
import logging
import json
from subprocess import Popen, PIPE


import dmdpy.utilities.utilities as utilities

logger = logging.getLogger(__name__)


class setupDMDjob:

    def __init__(self, parameters: dict=None, dir: str="./"):
        logger.debug("Entered setupDMDjob ")

        # Instance variables
        logger.debug("Initializing Variables")
        with open('dmdinput.json') as inp:
            self._raw_parameters = json.load(inp)

        self._run_directory = dir
        self._initial_directory = os.getcwd()
        self._dmd_config = utilities.load_dmd_config()

        self._protein = utilities.load_pdb("test.pdb")
        self.make_topparam()
        self.make_inConstr()
        self.make_state_file()
        self.make_start_file()

    def make_topparam(self):
        with open('topparam', 'w') as topparam_file:
            for residue in self._protein.sub_chain.residues:
                utilities.make_mol2(residue)
                topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2")

    def make_inConstr(self):
        # Need to also have the user constraints placed into here at some point, but that can wait for the time being!
        with open('inConstr', 'a') as inConstr_file:
            with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                       stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                       env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stderr.readline().strip())

        # TODO: add the user specified constraints now!
        logger.debug("Finished making the inConstr file!")

    def make_start_file(self):
        with open("dmd_start", 'w') as dmdstart:
            dmdstart.write(f"THERMOSTAT     {self._raw_parameters['Thermostat']}\n")
            dmdstart.write(f"T_NEW          {self._raw_parameters['Initial Temperature']}\n")
            dmdstart.write(f"T_LIMIT        {self._raw_parameters['Final Temperature']}\n")
            dmdstart.write(f"HEAT_X_C       {self._raw_parameters['HEAT_X_C']}\n")
            dmdstart.write(f"RESTART_FILE   {self._raw_parameters['Restart File']}\n")
            dmdstart.write(f"RESTART_DT     {self._raw_parameters['Restart dt']}\n")
            dmdstart.write(f"ECHO_FILE      {self._raw_parameters['Echo File']}\n")
            dmdstart.write(f"ECHO_DT        {self._raw_parameters['Echo dt']}\n")
            dmdstart.write(f"MOVIE_FILE     {self._raw_parameters['Movie File']}\n")
            dmdstart.write(f"START_TIME     {self._raw_parameters['Start time']}\n")
            dmdstart.write(f"MOVIE_DT       {self._raw_parameters['Movie dt']}\n")
            dmdstart.write(f"MAX_TIME       {self._raw_parameters['Max time']}\n")

        logger.debug("made the start file!")

    def make_state_file(self):
        with Popen(f"complex-1.linux -P {self._dmd_config['PATHS']['parameters']} -I {self._protein.name} -T topparam -D 200 -p param -s state -S 123 -C inConstr -c outConstr",
                   stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())
                logger.debug(shell.stderr.readline().strip())

        logger.debug("Made the state file!")
