#!/usr/bin/env python3

import logging
import os
import json

import dmdpy.protein.protein as protein
import dmdpy.utility.utilities as utilities
from dmdpy.setupjob import setupDMDjob

logger=logging.getLogger(__name__)

__all__ = [
    'calculation'
]

class calculation:

    def __init__(self, cores, run_dir: str='./', time=-1, pro: protein.Protein=None, parameters: dict=None, commands: dict=None):
        logger.info("Beginning DMD calculation")

        logger.debug("Initializing variables")
        self._submit_directory = os.getcwd()
        self._scratch_directory = run_dir
        self._config = utilities.load_dmd_config()
        self._cores = cores
        self._time_to_run = time
        self._timer_went_off = False
        self._dmd_config = utilities.load_dmd_config()
        os.environ["PATH"] += os.pathsep + self._dmd_config["PATHS"]["DMD_DIR"]

        # Want to make sure that we make the scratch directory!!
        try:
            logger.debug(f"Checking if scratch directory: {self._scratch_directory} is created")
            if not os.path.isdir(self._scratch_directory) and os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
                logger.info(f"Creating scratch directory: {self._scratch_directory}")
                os.mkdir(self._scratch_directory)

        except OSError:
            logger.warning("Could not make scratch directory : {self._scratch_directory}")
            logger.warning("I will run job in the current directory.")
            self._scratch_directory = './'

        logger.debug("Setting up DMD Environment")

        if parameters is None:
            logger.debug("Checking for a dmdinput.json file")
            if os.path.isfile("dmdinput.json"):
                self._parameter_file = os.path.join(self._submit_directory, "dmdinput.json")

            else:
                logger.error("No parameters specified for the job!")
                raise FileNotFoundError("dmdinput.json")

            #Now we read in the parameters here to a dictionary
            try:
                logger.debug("Loading in parameters")
                with open(self._parameter_file, 'r') as inputfile:
                    self._raw_parameters = json.load(inputfile)

            except IOError:
                logger.exception("Could not open the parameter file correctly!")
                raise

        else:
            logger.debug("Using parameters passed")
            self._raw_parameters = parameters

        # TODO check for valid parameters passed

        if not os.path.isfile("initial.pdb"):
            logger.info("initial.pdb not found, will try setting up from scratch")
            sj = setupDMDjob(parameters=self._raw_parameters)

        #We have an initial.pdb, now we have to check the jobs and update the parameters and then create the files and continue

        if commands is None:
            self.run_dmd(parameters)

        else:
            if commands["continue"]:
                if os.path.isfile(self._raw_parameters["Restart File"]):
                    # We continue from we last left off and then continue forward with the commands
                    pass

                #else:


    def run_dmd(self, parameters):
        # This will create the state file and the other files!
        utilities.make_state_file(parameters, "initial.pdb")
        utilities.make_start_file(parameters)

        #Now we execute the command to run the dmd


        

