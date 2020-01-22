#!/usr/bin/env python3

import os
import logging
import json
import pkg_resources
import shutil
import subprocess
from subprocess import Popen, PIPE


import dmdpy.utility.utilities as utilities
import dmdpy.protein.protein as protein

logger = logging.getLogger(__name__)

__all__=[
    'setupDMDjob'
]


class setupDMDjob:

    def __init__(self, parameters: dict=None, dir: str="./", pro: protein.Protein=None):
        logger.debug("Entered setupDMDjob")

        # Instance variables
        logger.debug("Initializing Variables")
        self._run_directory = dir
        logger.debug(f"Set run directory to: {self._run_directory}")
        self._initial_directory = os.getcwd()
        logger.debug(f"Set initial directory to: {self._initial_directory}")
        self._dmd_config = utilities.load_dmd_config()
        os.environ["PATH"] += os.pathsep + self._dmd_config["PATHS"]["DMD_DIR"]
        self._protein = None

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
                try:
                    shutil.copy(pkg_resources.resource_filename('dmdpy.resources', 'dmdinput.json'), './')

                except OSError:
                    logger.exception("Could not copy over default dmdinput.json!")
                    raise

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

        if pro is None:
            logger.debug("Checking for a pdb")
            files = [f for f in os.listdir('./') if os.path.isfile(os.path.join('./', f))]
            logger.debug(f"The files in {self._run_directory} are :")
            logger.debug(files)
            for f in files:
                if os.path.splitext(f)[-1].lower() == ".pdb":
                    logger.debug(f"Found a pdb file to use: {f}")
            # We store the residue here and then we backtrack and figure out the correct atom to protonate later after
            # proper relabeling/renaming

        self._displacement = []
        for atom_pair in self._raw_parameters["Restrict Displacement"]:
            self._displacement.append([self._protein.get_atom(atom_pair[0]), self._protein.get_atom(atom_pair[1]), atom_pair[2]])

        #TODO get the chains and the atoms in the frozen atoms here...
        self._static = {"chains": [], "residues": [], "atoms": []}
        if self._raw_parameters["Frozen atoms"]:
            for chain in self._raw_parameters["Frozen atoms"]["Chains"]:
                self._static["chains"].append(self._protein.get_chain(chain))

            #TODO Finish this

        logger.debug("Changing protein name to initial.pdb and writing out")
        self._protein.reformat_protein()
        self._protein.name = 'initial.pdb'
        self._protein.write_pdb()

        self.make_topparam()
        self.make_inConstr()
        utilities.make_state_file(self._raw_parameters, self._protein.name)
        self.short_dmd()
        utilities.make_start_file(self._raw_parameters)

        logger.info("#####################")
        logger.info("##                 ##")
        logger.info("##    SUCCESS!!    ##")
        logger.info("##                 ##")
        logger.info("#####################")

    def short_dmd(self):
        try:
            with open("dmd_start_short", 'w') as dmdstart:
                dmdstart.write(f"THERMOSTAT     {self._raw_parameters['Thermostat']}\n")
                dmdstart.write(f"T_NEW          0.001\n")
                dmdstart.write(f"T_LIMIT        0.001\n")
                dmdstart.write(f"HEAT_X_C       1\n")
                dmdstart.write(f"RESTART_FILE   {self._raw_parameters['Restart File']}\n")
                dmdstart.write(f"RESTART_DT     1\n")
                dmdstart.write(f"ECHO_FILE      {self._raw_parameters['Echo File']}\n")
                dmdstart.write(f"ECHO_DT        1\n")
                dmdstart.write(f"MOVIE_FILE     {self._raw_parameters['Movie File']}\n")
                dmdstart.write(f"START_TIME     0\n")
                dmdstart.write(f"MOVIE_DT       1\n")
                dmdstart.write(f"MAX_TIME       1\n")

        except IOError:
            logger.exception("Error writing out dmd_start file")
            raise

        # Here we do a short run
        try:
            with Popen(f"pdmd.linux -i dmd_start_short -s state -p param -c outConstr -m 1",
                    stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stdout.readline().strip())

        except OSError:
            logger.exception("Error calling pdmd.linux")
            raise

        if not os.path.isfile("movie"):
            logger.error("movie file was not made, dmd seems to be anrgy at your pdb")
            raise ValueError("initial.pdb")

        logger.debug("Finished the short DMD step successfully")

        utilities.make_movie("initial.pdb", "movie", "check.pdb")

        if not os.path.isfile("check.pdb"):
            logger.error("check.pdb not found, complex_M2P.linux did not run properly")
            raise ValueError("check.pdb")

        else:
            with open("check.pdb") as pdb:
                if len(pdb.readlines()) == 0:
                    logger.error("check.pdb is empty, complex_M2P.linux did not run properly")
                    raise ValueError("check.pdb")

        logger.debug("Was able to create a pdb from the short DMD run")
        logger.debug("Good to go, removing old file now")
        os.remove("check.pdb")
        os.remove(self._raw_parameters['Movie File'])
        os.remove(self._raw_parameters['Echo File'])
        os.remove(self._raw_parameters['Restart File'])
        os.remove("dmd_start_short")

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

                    topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2\n")

        except IOError:
            logger.exception("Error with writing topparam file!")
            raise

    def make_inConstr(self):
        try:
            with open('inConstr', 'w') as inConstr_file:
                try:
                    with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                               stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                               env=os.environ) as shell:
                        while shell.poll() is None:
                            if len(shell.stderr.readline().strip()) > 0:
                                logger.debug(shell.stderr.readline().strip())
                except OSError:
                    logger.exception("Error calling genESC.linux")
                    raise

                if self._raw_parameters["Freeze Non-Residues"]:
                    logger.debug("Freeze Non-residues turned on, freezing residues")
                    for residue in self._protein.sub_chain.residues:
                        logger.debug(f"Freezing residue: {residue}")
                        inConstr_file.write(f"Static {residue.write_inConstr()}\n")

                #TODO, I don't believe I every specify how to create this static list...
                for static_chain in self._static["chains"]:
                    logger.debug(f"Freezing chain: {static_chain}")
                    inConstr_file.write(f"Static {static_chain.write_inConstr()}\n")

                for static_residue in self._static["residues"]:
                    logger.debug(f"Freezing residue: {static_residue}")
                    inConstr_file.write(f"Static {static_residue.write_inConstr()}\n")

                for static_atom in self._static["atoms"]:
                    logger.debug(f"Freezing atom: {static_atom}")
                    inConstr_file.write(f"Static {static_atom.write_inConstr()}\n")

                # TODO Fix This for custom protonation states

                # for pro_atom in self._protonation_states:
                #     if pro_atom[1].lower() == "protonate":
                #         logger.debug(f"Protonating atom: {pro_atom[0]}")
                #         inConstr_file.write(f"Protonate {pro_atom[0].write_inConstr()}\n")
                #
                #     elif pro_atom[1].lower() == "deprotonate":
                #         logger.debug(f"Deprotonating atom: {pro_atom[0]}")
                #         inConstr_file.write(f"Deprotonate {pro_atom[0].write_inConstr()}\n")
                #

                # TODO include restrict Metal ligands (find all atoms that are some distance away from the metals)!
                for disp_atom in self._displacement:
                    logger.debug(
                        f"Restricting motion of atom: {disp_atom[0]} and atom {disp_atom[1]} by {disp_atom[2]}")
                    inConstr_file.write(
                        f"AtomPairRel {disp_atom[0].write_inConstr()} {disp_atom[1].write_inConstr()} -{disp_atom[2]} +{disp_atom[2]}\n")

        except IOError:
            logger.exception("Error opening inConstr file")
            raise

        logger.debug("Finished making the inConstr file!")
