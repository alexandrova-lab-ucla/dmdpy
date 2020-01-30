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
from dmdpy.utility.exceptions import ParameterError
import dmdpy.utility.constants as constant

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

        try:
            utilities.valid_parameters(self._raw_parameters)

        except ValueError:
            logger.exception("Missing a parameter definition!")
            raise ValueError("definition")

        except ParameterError:
            logger.exception("Invalid parameter specification")
            raise

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
                        continue
                    break
            # We store the residue here and then we backtrack and figure out the correct atom to protonate later after
            # proper relabeling/renaming

        else:
            self._protein = pro

        if self._protein is None:
            raise ValueError("No Protein!")

        self._displacement = []
        if "Restrict Displacement" in self._raw_parameters.keys():
            for atom_pair in self._raw_parameters["Restrict Displacement"]:
                self._displacement.append([self._protein.get_atom(atom_pair[0]), self._protein.get_atom(atom_pair[1]), atom_pair[2]])

        self._static = []
        if "Frozen atoms" in self._raw_parameters.keys():
            for chain in self._raw_parameters["Frozen atoms"]["Chains"]:
                try:
                    protein_chain =self._protein.get_chain(chain)
                    for residue in protein_chain.residues:
                        self._static.extend(residue.atoms)

                except ValueError:
                    logger.exception("Could not find the chain!")
                    raise

            for residue in self._raw_parameters["Frozen atoms"]["Residues"]:
                try:
                    #TODO check to see if the residue HAS a metal in it...if so it will not work!!!!!!
                    #self._static["residues"].append(self._protein.get_residue(residue))
                    protein_residue = self._protein.get_residue(residue)
                    self._static.extend(protein_residue.atoms)

                except ValueError:
                    logger.exception("Could not find the residue!")
                    raise

            for atom in self._raw_parameters["Frozen atoms"]["Atoms"]:
                try:
                    self._static.append(self._protein.get_atom(atom))

                except ValueError:
                    logger.exception("Could not find the atom!")
                    raise

        # These holds all of the residues with weird protonation or deprotonation states
        self._protonate = []
        if "Custom protonation states" in self._raw_parameters.keys():
            for item in self._raw_parameters["Custom protonation states"]:
                id = [item[0], item[1]]
                self._protonate.append([self._protein.get_residue(id), item[2:]])

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
                # This is apparently not needed
                # try:
                #     with Popen(f"genESC.linux {self._dmd_config['PATHS']['parameters']} {self._protein.name} topparam",
                #                stdout=inConstr_file, stderr=PIPE, universal_newlines=True, shell=True, bufsize=1,
                #                env=os.environ) as shell:
                #         while shell.poll() is None:
                #             if len(shell.stderr.readline().strip()) > 0:
                #                 logger.debug(shell.stderr.readline().strip())
                # except OSError:
                #     logger.exception("Error calling genESC.linux")
                #     raise

                if self._raw_parameters["Freeze Non-Residues"]:
                    logger.debug("Freeze Non-residues turned on, freezing residues")
                    for residue in self._protein.sub_chain.residues:
                        logger.debug(f"Freezing residue: {residue}")
                        inConstr_file.write(f"Static {residue.write_inConstr()}\n")

                for static_atom in self._static:
                    logger.debug(f"Freezing atom: {static_atom}")
                    inConstr_file.write(f"Static {static_atom.write_inConstr()}\n")

                for state in self._protonate:
                    logger.debug(f"Adding protonation state: {state[0]} and {state[1]}")
                    atom_id = ""
                    #TODO try and except for weird atoms or residues if it cannot find it
                    if len(state[1]) > 1:
                        #Then we had a number specify
                        logger.debug("Specified which atom specifically to use!")
                        if state[1][0] == "protonate":
                            atom_id = constant.PROTONATED[state[0].name][state[1][1]]

                        elif state[1][0] == "deprotonate":
                            atom_id = constant.DEPROTONATED[state[0].name][state[1][1]]

                    else:
                        if state[1][0] == "protonate":
                            atom_id = constant.PROTONATED[state[0].name][0]

                        elif state[1][0] == "deprotonate":
                            atom_id = constant.DEPROTONATED[state[0].name][0]

                    if atom_id == "":
                        raise ValueError("Did not specify to protonate or deprotonate correctly")

                    try:
                        pro_atom = state[0].get_atom(atom_id[0])

                    except ValueError:
                        logger.exception("Could not find the correct atom to protonate or deprotonate in the residue")
                        logger.warning(f"{state}")
                        raise

                    if state[1][0] == "protonate":
                        inConstr_file.write(f"Protonate {pro_atom.write_inConstr()}\n")

                    else:
                        inConstr_file.write(f"Deprotonate {pro_atom.write_inConstr()}\n")

                if self._raw_parameters["Restrict Metal Ligands"]:
                    logger.debug("Restricting distance between atoms and metals!")
                    for metal in self._protein.metals:
                        logger.debug(f"Looking at metal: {metal}")
                        atoms_near_metal = self._protein.atoms_near_metal(metal, 3.1)
                        for atoms in atoms_near_metal:
                            logger.debug(f"Freezing atom: {atoms} since too close to a metal")
                            inConstr_file.write(f"Static {atoms.write_inConstr()}\n")

                            for bonded_atoms in atoms.bonds:
                                if bonded_atoms.element != "h" and bonded_atoms.element not in constant.METALS:
                                    logger.debug(f"Restricting motion of atom {bonded_atoms} and atom {metal} by {0.05}")
                                    inConstr_file.write(
                                        f"AtomPairRel {bonded_atoms.write_inConstr()} {metal.write_inConstr()} -{0.05} +{0.05}\n")

                for disp_atom in self._displacement:
                    logger.debug(
                        f"Restricting motion of atom: {disp_atom[0]} and atom {disp_atom[1]} by {disp_atom[2]}")
                    inConstr_file.write(
                        f"AtomPairRel {disp_atom[0].write_inConstr()} {disp_atom[1].write_inConstr()} -{disp_atom[2]} +{disp_atom[2]}\n")

        except IOError:
            logger.exception("Error opening inConstr file")
            raise

        logger.debug("Finished making the inConstr file!")

    def updated_parameters(self):
        #TODO make sure that this is correct
        new_parameters = self._raw_parameters.copy()

        # Update the custom protonation states
        for new_state, state in zip(new_parameters["Custom protonation states"], self._protonate):
            new_state[0] = state[0].chain.name
            new_state[1] = state[0].residue.number
            new_state[2] = state[1][0]
            if len(new_state) == 4:
                new_state[3] = state[1][1]

        # Update the frozen atoms
        for new_state, state in zip(new_parameters["Frozen atoms"]["Chains"], self._static["chains"]):
            new_state = state.name

        for new_state, state in zip(new_parameters["Frozen atoms"]["Residues"], self._static["residues"]):
            new_state[0] = state.chain.name
            new_state[1] = state.number

        for new_state, state in zip(new_parameters["Frozen atoms"]["Atoms"], self._static["atoms"]):
            new_state[0] = state.chain.name
            new_state[1] = state.residue.number
            new_state[2] = state.id

        # Update the displacement atoms
        for new_state, state in zip(new_parameters["Restrict Displacement"], self._displacement):
            new_state[0][0] = state[0].chain.name
            new_state[0][1] = state[0].residue.number
            new_state[0][2] = state[0].id
            new_state[1][0] = state[1].chain.name
            new_state[1][1] = state[1].residue.number
            new_state[1][2] = state[1].id
            new_state[2] = state[2]

        return new_parameters
