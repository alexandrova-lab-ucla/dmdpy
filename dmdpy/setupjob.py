#!/usr/bin/env python3

import os
import logging
import json
import pkg_resources
import shutil
import subprocess
from subprocess import Popen, PIPE


import dmdpy.utilities.utilities as utilities
import dmdpy.protein.protein as protein

logger = logging.getLogger(__name__)

os.environ["PATH"] += os.pathsep + "/home/mhennefarth/Programs/dmd/bin"

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
                    try:
                        self._protein = utilities.load_pdb(f)

                    except:
                        logger.exception("You idiot!")
                        raise

            if self._protein is None:
                logger.error("No pdb file in the run directory")
                raise ValueError

        else:
            self._protein = pro

        # These should be pointers to these objects so that if they change, it is updated in this list automatically
        # TODO: change this to the way it is in dmdinput.json
        self._static = []
        for static_atom in self._raw_parameters["Frozen atoms"]:
            self._static.append(self._protein.get_atom(static_atom))

        self._protonation_states = []
        for prot_atom in self._raw_parameters["Custom protonation states"]:
            self._protonation_states.append([self._protein.get_atom(prot_atom[:2]), prot_atom[3]])

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
        self.short_dmd()
        self.make_start_file()

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
                dmdstart.write(f"MAX_TIME       5\n")

        except IOError:
            logger.exception("Error writing out dmd_start file")
            raise

        # Here we do a short run
        # TODO wrap in some try/except clauses
        with Popen(
                f"pdmd.linux -i dmd_start_short -s state -p param -c outConstr -m 1",
                stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

        if not os.path.isfile("movie"):
            logger.error("movie file was not made, dmd seems to be anrgy at your pdb")
            raise ValueError("initial.pdb")

        logger.debug("Finished the short DMD step successfully")

        with Popen(
                f"complex_M2P.linux {self._dmd_config['PATHS']['parameters']} initial.pdb topparam movie check.pdb inConstr",
                stdout=PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

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
        # TODO: wrap this is sexy try/except
        with open('inConstr', 'a') as inConstr_file:
            with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                       stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                       env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stderr.readline().strip())

            for static_atom in self._static:
                logger.debug(f"Freezing atom: {static_atom}")
                inConstr_file.write(f"Static {static_atom.write_inConstr()}\n")

            for pro_atom in self._protonation_states:
                if pro_atom[1].lower() == "protonate":
                    logger.debug(f"Protonating atom: {pro_atom[0]}")
                    inConstr_file.write(f"Protonate {pro_atom[0].write_inConstr()}\n")

                elif pro_atom[1].lower() == "deprotonate":
                    logger.debug(f"Deprotonating atom: {pro_atom[0]}")
                    inConstr_file.write(f"Deprotonate {pro_atom[0].write_inConstr()}\n")

            for disp_atom in self._displacement:
                logger.debug(f"Restricting motion of atom: {disp_atom[0]} and atom {disp_atom[1]} by {disp_atom[2]}")
                inConstr_file.write(f"AtomPairRel {disp_atom[0].write_inConstr()} {disp_atom[1].write_inConstr()} -{disp_atom[2]} +{disp_atom[2]}\n")

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
        logger.debug("Calling complex.linux")
        try:
            # TODO: There is an issue here with complex.linux not actually running
            # It is coming from the mol2 of the substrate...interesting...
            # There is a Segmentation fault (core dumped) error that occurs, asking Jack if he knows what the issue is
            # Jack thinks it is the segfault mike wrote in his HACK ALERT section
            # I have emailed the Dohkyan group regarding it...its only for certain pdbs...
            with Popen(f"complex-1.linux -P {self._dmd_config['PATHS']['parameters']} -I {self._protein.name} -T topparam -D 200 -p param -s state -C inConstr -c outConstr",
                       stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stdout.readline().strip())

            attempts = 3
            if not os.path.isfile("state"):
                logger.warning("Could not make state file first time around, could be a segfault error")
                logger.warning("Going to reorder the bond list in the mol2 files!")
                attempts = 1

            while not os.path.isfile("state") and attempts < 3:
                logger.warning(f"complex fix attempt {attempts}")
                mol2_files = []
                with open("topparam") as topparm:
                    for line in topparm:
                        mol2_files.append(line.split()[2])

                for mol2 in mol2_files:
                    with open(mol2, 'r') as mf:
                        bonds = []
                        save = []
                        bond_section = False
                        for line in mf:
                            if not bond_section and "BOND" in line:
                                save.append(line)
                                bond_section = True
                                continue

                            if bond_section:
                                bonds.append(line.split())

                            else:
                                save.append(line)

                        #now we reorder based upon the attempt!
                        bonds.sort(key=lambda x: int(x[attempts]))

                    with open(mol2, 'w+') as mf:
                        for line in save:
                            mf.write(line)

                        for bond in bonds:
                            mf.write(f"{bond[0]}\t{bond[1]}\t{bond[2]}\t{bond[3]}\n")

                    with Popen(f"complex-1.linux -P {self._dmd_config['PATHS']['parameters']} -I {self._protein.name} -T topparam -D 200 -p param -s state -C inConstr -c outConstr",
                            stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                        while shell.poll() is None:
                            logger.debug(shell.stdout.readline().strip())

                attempts += 1

            if not os.path.isfile("state"):
                logger.critical("could not create state file, something is very wrong!")
                raise ValueError("state file")

        except OSError:
            logger.exception("Error calling complex-1.linux!")
            raise

        logger.debug("Made the state file!")
