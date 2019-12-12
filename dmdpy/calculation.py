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

    def __init__(self, cores: int = 1, run_dir: str='./', time=-1, pro: protein.Protein=None, parameters: dict=None):
        logger.info("Beginning DMD calculation")

        logger.debug("Initializing variables")
        self._submit_directory = os.getcwd()
        self._scratch_directory = run_dir
        self._config = utilities.load_dmd_config()
        self._cores = cores
        self._time_to_run = time
        self._timer_went_off = False
        self._dmd_config = utilities.load_dmd_config()
        self._start_time = 0
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

        if os.path.isfile(self._raw_parameters["Echo File"]):
            with open(self._raw_parameters["Echo File"]) as echofile:
                lines = []
                for line in echofile:
                    lines.append(line)

                last_line = lines[-1].split()
                self._start_time = int(float(last_line[0]))
                logger.debug(f"Last recorded time: {self._start_time}")

        if self._raw_parameters["Remaining Commands"]:
            self._commands = self._raw_parameters["Remaining Commands"].copy()

            if self._raw_parameters["Commands"]:
                all_commands = self._raw_parameters["Commands"].copy()

                remove = len(self._commands)
                for i in range(remove):
                    all_commands.pop(list(self._commands.keys())[-1])

                time_elapsed = 0
                for step in all_commands:
                    if "Time" in all_commands[step].keys():
                        time_elapsed += all_commands[step]["Time"]

                    else:
                        time_elapsed += self._raw_parameters["Time"]

                logger.debug(f"Total time completed: {time_elapsed}")

                diff = self._start_time - time_elapsed
                logger.debug(f"We are off by: {diff}")

                if "Time" in self._commands[list(self._commands.keys())[0]].keys():
                    new_time = self._commands[list(self._commands.keys())[0]]["Time"] - diff

                else:
                    new_time = self._raw_parameters["Time"] - diff

                if new_time < 0:
                    logger.error("Somehow we moved onto a later step then what is reported in remaining calculations.")
                    raise ValueError("Invalid time")

                logger.debug(f"Setting new time for the first step to: {new_time}")
                self._commands[list(self._commands.keys())[0]]["Time"] = new_time

            else:
                logger.warning("Unknown how many steps prior to this one!")
                logger.warning("Will just start continue from where we left off then")

        elif self._raw_parameters["Commands"]:
            # There are no remaining commands, so continue like normal more or less
            self._commands = self._raw_parameters["Commands"].copy()

        else:
            self._commands = {"1": {}}
            logger.info("Commands passed, using those")

        # TODO move to the scratch directory

        # We can arm the timer
        if self._time_to_run != -1:
            logger.info("Starting the timer")
            signal.signal(signal.SIGALRM, self.calculation_alarm_handler)
            signal.alarm((self._time_to_run * 60 - 30) * 60)

        # We loop over the steps here and will pop elements off the beginning of the dictionary
        while len(self._commands.values()) != 0:
            if self._timer_went_off:
                logger.info("Timer went off, not continuing onto next command")
                break

            # Grab the next step dictionary to do
            steps = self._commands[list(self._commands.keys())[0]]
            logger.info(f"On step: {steps}")
            updated_parameters = self._raw_parameters.copy()

            for changes in steps:
                logger.debug(f"Updating {changes}: changing {updated_parameters[changes]} to {steps[changes]}")
                updated_parameters[changes] = steps[changes]

            if updated_parameters["titr"]["titr on"]:
                logger.critical("Not implemented yet, dingus")
                pass

            elif "Custom protonation states" in steps.keys():
                logger.critical("Not implemented yet, dingus")
                pass

            elif "Frozen atoms" in steps.keys() or "Restrict Displacement" in steps.keys():
                #TODO implement
                #make new inConstr file, run complex-linux.1, and then use the restart file
                logger.critical("Not implemented yet, dingus")

                #self.run_dmd(updated_parameters, self._start_time, True)

            else:
                # We can just run the job with no issues other than those raised from the above
                self.run_dmd(updated_parameters, self._start_time, True)


            # Assuming we finished correctly, we pop off the last issue
            self._commands.pop(list(self._commands.keys())[0])
            # Update the new start time!
            self._start_time += updated_parameters["Time"]

        #TODO
        # Now we save the remaining commands and transfer everything back and forth between the necessary locations!

        #TODO change back to initial directory if necessary

    def run_dmd(self, parameters, start_time: int, use_restart: bool):
        # Remake the start file with any changed parameters
        utilities.make_start_file(parameters, start_time)

        if use_restart:
            state_file = self._raw_parameters["Restart File"] if os.path.isfile(self._raw_parameters["Restart File"]) else "state"

        else:
            state_file = "state"

        #Now we execute the command to run the dmd
        try:
            logger.info(f"Issuing command: pdmd.linux -i dmd_start -s {state_file} -p param -c outConstr -m {self._cores} -fa")
            logger.info("*****************************************************************************************************")
            with Popen(f"pdmd.linux -i dmd_start -s {state_file} -p param -c outConstr -m {self._cores} -fa",
                    stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                while shell.poll() is None:
                    logger.info(shell.stdout.readline().strip())
            logger.info("*****************************************************************************************************")
        except OSError:
            logger.exception("Error calling pdmd.linux")
            raise

    def last_frame(self):
        #TODO implement
        """Returns a protein of the last from from the movie file"""
        pass

    def calculation_alarm_handler(self, signum, frame):
        """
        Called if time is almost up in the dmd job! Will copy all files to a backup directory and write out the
        remaining commands to perform in the remaining_commands.json file.
        """
        logger.warning("Creating a backup directory!")
        self._timer_went_off = True

        logger.info("Placing the remaining commands into remaining_commands.json")
        with open("remaining_commands.json", 'w') as rc:
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