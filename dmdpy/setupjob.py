#!/usr/bin/env python3

import os
import logging
import json
from subprocess import Popen, PIPE


import dmdpy.utilities.utilities as utilities
import dmdpy.protein.protein as protein

logger = logging.getLogger(__name__)


class setupDMDjob:

    # Have it do a small DMD step after everything is donew

    def __init__(self, parameters: dict=None, dir: str="./", pro: protein.Protein=None):
        logger.debug("Entered setupDMDjob")

        # Instance variables
        logger.debug("Initializing Variables")
        self._run_directory = dir
        logger.debug(f"Set run directory to: {self._run_directory}")
        self._initial_directory = os.getcwd()
        logger.debug(f"Set initial directory to: {self._initial_directory}")
        self._dmd_config = utilities.load_dmd_config()

        try:
            logger.debug(f"Changing to run directory: {self._run_directory}" )
            os.chdir(self._run_directory)

        except OSError:
            logger.exception(f"Failed moving to run directory: {self._run_directory}")
            raise

        # Read in the parameters file, or use the passed parameters
        if parameters is None:
            logger.debug("Searching for a dmdinput.json file")
            if not os.path.isfile("dmdinput.json"):
                logger.error("Could not find dmdinput.json")
                logger.error("I will copy over a default dmdinput.json for you to edit")
                # TODO have it copy over a default dmdinput.json
                raise ValueError("dmdinput.json")

            try:
                with open('dmdinput.json') as inp:
                    self._raw_parameters = json.load(inp)

            except IOError:
                logger.exception("Error reading in dmdinput.json file!")
                raise

            except ValueError:
                logger.exception("dmdinput.json not formatted correctly")
                raise

        else:
            logger.debug("Using passed parameters")
            self._raw_parameters = parameters

        # TODO: check to see if parameters is good!

        # TODO change test.pdb to a searched pdb!
        if pro is None:
            self._protein = utilities.load_pdb("test.pdb")

        else:
            self._protein = pro

        # TODO have it reupdate the/add params to the atoms via inConstr file stuff prior to reformatting pdb

        # These should be pointers to these objects so that if they change, it is updated in this list automatically
        self._static = []
        for static_atom in self._raw_parameters["Frozen atoms"]:
            self._static.append(self._protein.get_atom(static_atom))

        self._protonation_states = []
        for prot_atom in self._raw_parameters["Custom protonation states"]:
            self._protonation_states.append([self._protein.get_residue(prot_atom[:1]), prot_atom[2]])

        self._displacement = []
        for atom_pair in self._raw_parameters["Restrict Displacement"]:
            self._displacement.append([self._protein.get_atom(atom_pair[0]), self._protein.get_atom(atom_pair[1]), atom_pair[2]])

        logger.debug("Changing protein name to initial.pdb and writing out")
        self._protein.reformat_protein()
        self._protein.name = 'initial.pdb'
        self._protein.write_pdb()

        self.make_topparam()
        self.make_inConstr()
        self.make_state_file()
        self.make_start_file()

        # TODO Have it run a quick DMD step to ensure it works properly
        # Can also use this to reformat the pdb!

    def make_topparam(self):
        try:
            logger.debug("Making topparam file")
            with open('topparam', 'w') as topparam_file:
                for residue in self._protein.sub_chain.residues:
                    try:
                        utilities.make_mol2(residue)

                    except OSError:
                        logger.error("Error in making mol2 file")
                        raise

                    topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2")

        except IOError:
            logger.exception("Error with writing topparam file!")
            raise

    def make_inConstr(self):
        # Need to also have the user constraints placed into here at some point, but that can wait for the time being!
        # TODO: wrap this is sexy try/except
        with open('inConstr', 'a') as inConstr_file:
            with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                       stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                       env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stderr.readline().strip())

        # TODO: add the user specified constraints now!
        logger.debug("Finished making the inConstr file!")

    def make_start_file(self):
        try:
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

        except IOError:
            logger.exception("Error writing out dmd_start file")
            raise

        logger.debug("made the start file!")

    def make_state_file(self):
        logger.debug("Calling complex-1.linux")
        try:
            with Popen(f"complex-1.linux -P {self._dmd_config['PATHS']['parameters']} -I {self._protein.name} -T topparam -D 200 -p param -s state -S 123 -C inConstr -c outConstr",
                       stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stdout.readline().strip())
                    logger.debug(shell.stderr.readline().strip())

        except OSError:
            logger.exception("Error calling complex-1.linux!")
            raise

        logger.debug("Made the state file!")
