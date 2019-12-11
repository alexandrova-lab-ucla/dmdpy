#!/usr/bin/env python3

import logging
import os
import shutil
import json
import signal
import subprocess
from subprocess import Popen, PIPE

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
        self._commands = commands
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

        # TODO check for any exceptions raised from setupDMDjob
        if pro is None:
            if not os.path.isfile("initial.pdb"):
                logger.info("initial.pdb not found, will try setting up from scratch")
                sj = setupDMDjob(parameters=self._raw_parameters)

        else:
            logger.debug("Will setup the protein for DMD")
            sj = setupDMDjob(parameters=self._raw_parameters, pro=pro)


        # TODO move to the scratch directory

        # TODO Arm the timer

        # We loop over the steps here and will pop elements off the beginning of the dictionary
        while len(self._commands.values()) != 0:
            if self._timer_went_off:
                break

            # Grab the next step dictionary to do
            steps = self._commands[self._commands.values(0)]
            updated_parameters = self._raw_parameters

            for changes in steps.keys():
                if changes == "continue":
                    continue

                updated_parameters[changes] = steps[changes]

            if not steps["continue"]:
                # We need to reset up the job, first convert the movie to a pdb and then remake everything
                try:
                    with Popen(
                            f"complex_M2P.linux {self._dmd_config['PATHS']['parameters']} initial.pdb topparam {self._raw_parameters['Movie File']} initial.pdb inConstr",
                            stdout=PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True, shell=True,
                            env=os.environ) as shell:
                        while shell.poll() is None:
                            logger.debug(shell.stdout.readline().strip())
                except OSError:
                    logger.exception("Error calling complex_M2P.linux")
                    raise

                # Now we resetup intermediate.pdb to initial.pdb
                self.run_dmd(updated_parameters, False)

            else:
                # We can just run the job with no issues other than those raised from the above
                self.run_dmd(updated_parameters, True)


            # Assuming we finished correctly, we pop off the last issue
            self._commands.pop(self._commands.values(0))

        # Now we save the remaining commands and transfer everything back and forth between the necessary locations!

    def run_dmd(self, parameters, skip: bool):
        # This will create the state file and the other files!
        if not skip:
            utilities.make_state_file(parameters, "initial.pdb")
            utilities.make_start_file(parameters)

        #Now we execute the command to run the dmd
        try:
            with Popen(f"pdmd.linux -i dmd_start -s state -p param -c outConstr -m {self._cores}",
                    stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stdout.readline().strip())

        except OSError:
            logger.exception("Error calling pdmd.linux")
            raise

    def calculation_alarm_handler(self, signum, frame):
        """
        Called if time is almost up in the dmd job! Will copy all files to a backup directory and write out the
        remaining commands to perform in the remaining_commands.json file.
        """
        logger.warning("Creating a backup directory!")
        self._timer_went_off = True

        logger.info("Placing the remaining commands into remaining_commands.json")
        with open("remaining_commands.json") as rc:
            json.dump(self._commands, rc)

        logger.debug("Checking if scratch directory is different from submit directory")
        if os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
            logger.warning("Creating dmd_backup in the submit directory")
            backup_dir = os.path.join(os.path.abspath(self._submit_directory), 'dmd_backup')
            if os.path.isdir(backup_dir):
                logger.warning("Removing backup directory already present in the submit directory")
                shutil.rmtree(backup_dir)

            logger.warning(f"Copying files from {os.path.abspath(self._scratch_directory)} to {backup_dir}")
            shutil.copytree(self._scratch_directory, backup_dir)

        logger.info("Turning off alarm")
        signal.alarm(0)