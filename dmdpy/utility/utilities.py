#!/usr/bin/env python3

import json
import pkg_resources
import logging
import shutil
import os
from logging.config import dictConfig
import subprocess
from subprocess import Popen, PIPE


from dmdpy.protein import atom, chain, residue, protein
import dmdpy.utility.constants as constants
from dmdpy.utility.exceptions import ParameterError

__all__=[
    'load_logger_config',
    'load_dmd_config',
    'load_pdb',
    'make_mol2',
    'create_config',
    'make_start_file',
    'make_state_file',
    'make_movie',
    'setup_dmd_environ',
    'valid_parameters'
]

logger = logging.getLogger(__name__)


def load_logger_config():
    """Loads in the Logger Config File"""

    home = os.path.expanduser('~')
    path_to_config = os.path.join(home, ".dmdpy/logger_config.json")

    if not os.path.isdir(os.path.join(home, '.dmdpy')):
        create_config()
        raise ValueError(".dmdpy missing")

    if not os.path.isfile(path_to_config):
        print("Logger file not in .dmdpy directory")
        print("Copying over default logger now")
        shutil.copy(pkg_resources.resource_filename('dmdpy.resources', 'logger_config.json'), os.path.join(home, '.dmdpy'))

    try:
        with open(path_to_config) as logger_config:
            config = json.load(logger_config)

        dictConfig(config)
        logger.debug("Loaded in logger parameters")
        logger.debug("Logger has started!")

    except IOError:
        logger.critical("Could not open logger_config.json")
        raise

    except ValueError:
        logger.critical("logger_config.json is not formatted properly")
        raise


def load_dmd_config():

    home = os.path.expanduser('~')
    path_to_config = os.path.join(home, ".dmdpy/dmd_config.json")

    if not os.path.isdir(os.path.join(home, '.dmdpy')):
        create_config()
        raise ValueError(".dmdpy missing")

    if not os.path.isfile(path_to_config):
        logger.critical("DMD_config file not in .dmdpy directory")
        logger.critical("Copying over default file now")
        logger.critical("Ensure that the file is correct before continuing!")
        shutil.copy(pkg_resources.resource_filename('dmdpy.resources', 'dmd_config.json'), os.path.join(home, '.dmdpy'))
        raise ValueError("dmd_config.json file missing")

    try:
        logger.debug("Loading in the dmd_config parameters")
        with open(pkg_resources.resource_filename('dmdpy.resources', 'dmd_config.json')) as logger_config:
            config = json.load(logger_config)

        return config

    except IOError:
        logger.critical("Could not open dmd_config.json")
        raise

    except ValueError:
        logger.critical("dmd_config.json is not formatted properly")
        raise

dmd_config = load_dmd_config()

def setup_dmd_environ():
    """ setups of an os.environ so that DMD can be run """
    logger.debug("Setting up DMD Environment")
    os.environ["PATH"] += os.pathsep + dmd_config["PATHS"]["DMD_DIR"]
    return os.environ

def valid_parameters(parameters: dict):
    #TODO Finish this
    """Checks to see if the dmd parameters passed are valid at all"""
    logger.debug("Checking if parameters are valid")
    if "Thermostat" not in parameters.keys():
        raise ParameterError("No Thermostat")

    elif parameters["Thermostat"] not in constants.THERMOSTATS:
        raise ValueError("Invalid thermostat")

    if "Initial Temperature" not in parameters.keys():
        raise ParameterError("No Initial Temperature")

    else:
        try:
            assert(parameters["Initial Temperature"].isdigit() and parameters["Initial Temperature"] > 0)

        except ValueError:
            raise ParameterError("Incorrect initial temperature provided")

    if "Final Temperature" not in parameters.keys():
        raise ParameterError("No Final Temperature")

    else:
        try:
            assert (parameters["Final Temperature"].isdigit() and parameters["Final Temperature"] > 0)

        except ValueError:
            raise ParameterError("Incorrect Final temperature provided")

    if "HEAT_X_C" not in parameters.keys():
        raise ParameterError("No HEAT_X_C")

    else:
        try:
            assert (parameters["HEAT_X_C"].isdigit() and parameters["HEAT_X_C"] > 0)

        except ValueError:
            raise ParameterError("Incorrect HEAT_X_C provided")

    if "Echo File" not in parameters.keys():
        raise ParameterError("You did not specify what to label the echo file")

    if "Movie File" not in parameters.keys():
        raise ParameterError("You did not specify what to call the movie file")

    if 'Restart File' not in parameters.keys():
        raise ParameterError("You did not specify what to call the restart file")

    if "dt" not in parameters.keys():
        raise ParameterError("You did not specify the time in which save data")

    else:
        try:
            assert(parameters["dt"].isdigit() and parameters["dt"] > 0)

        except ValueError:
            raise ParameterError(f"Invalid value provided for dt: {parameters['dt']}")

    if "Time" not in parameters.keys():
        raise ParameterError("You did not specify the default length of time to be for the DMD simulation")

    else:
        try:
            assert (parameters["Time"].isdigit() and parameters["Time"] > 0)

        except ValueError:
            raise ParameterError(f"Invalid value provided for Time: {parameters['Time']}")

    if "titr" in parameters.keys():
        try:
            assert(type(parameters["titr"]["titr on"]) == bool)

        except ValueError:
            raise ParameterError("You did not correctly type the titr on parameters, must be a bool")

    if "Freeze Non-Residues" not in parameters.keys():
        raise ValueError("Missing Freeze Non-Residues parameters")

    else:
        try:
            assert(type(parameters["Freeze Non-Residues"]) == bool)

        except ValueError:
            raise ParameterError("Freeze Non-Residues MUST be a bool")

    if "Restrict Metal Ligands" not in parameters.keys():
        raise ValueError("Missing Restrict Metal Ligands parameters")

    else:
        try:
            assert (type(parameters["Restrict Metal Ligands"]) == bool)

        except ValueError:
            raise ParameterError("Restrict Metal Ligands MUST be a bool")

    if "Custom protonation states" not in parameters.keys():
        raise ValueError("Missing Custom protonation state parameters")

    else:
        for state in parameters["Custom protonation states"]:
            try:
                assert(len(state) == 3)
                assert(type(state[0]) == str)
                assert(type(state[1]) == int and state[1] > 0)
                assert(type(state[2]) == str)

            except ValueError:
                raise ParameterError(f"Invalid specification of protonation states: {state}")

    if "Frozen atoms" in parameters.keys():
        if "Chains" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Chains in frozen atom definition")

        else:
            for chain in parameters["Frozen atoms"]["Chains"]:
                try:
                    assert(type(chain) == str)
                    assert(chain.isalpha())

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen chain: {chain}")

        if "Residues" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Residue in frozen atom definition")

        else:
            for residue in parameters["Frozen atoms"]["Residues"]:
                try:
                    assert(type(residue[0]) == str and residue[0].isalpha())
                    assert(type(residue[1]) == int and residue[1] > 0)

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen residue: {residue}")

        if "Atoms" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Atom in frozen atom definition")

        else:
            for atom in parameters["Frozen atoms"]["Atoms"]:
                try:
                    assert(type(atom[0]) is str and atom[0].isalpha())
                    assert(type(atom[1]) is int and atom[1] > 0)
                    assert(type(atom[2]) is str and atom[2].isalpha())

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen atom: {atom}")

    if "Restrict Displacement" in parameters.keys():
        for id in parameters["Restrict Displacement"]:
            try:
                assert(len( id) == 3)
                assert(type(id[2]) is float and id[2] > 0)

                assert(len(id[0]) == 3)
                assert(type(id[0][0]) is str and id[0][0].isalpha())
                assert(type(id[0][1]) is int and id[0][1] > 0)
                assert(type(id[0][2]) is str and id[0][0].isalpha())

                assert(len(id[1]) == 3)
                assert (type(id[1][0]) is str and id[1][0].isalpha())
                assert (type(id[1][1]) is int and id[1][1] > 0)
                assert (type(id[1][2]) is str and id[1][0].isalpha())

                assert(id[0] != id[1])

            except ValueError:
                raise ParameterError(f"Invalid definition of Restrict Displacement: {id}")

    if "Commands" not in parameters.keys():
        raise ValueError("Missing commands definition")

    else:
        for command in parameters["Commands"]:
            try:
                assert(type(command) is dict)

            except ValueError:
                raise ParameterError(f"Command provided is not a dictionary: {command}")


def load_pdb(file: str):
    logger.debug(f"Finding file: {file}")

    chains = []

    resNum = 0
    chainLet = ""
    lineNumber = 1
    try:
        with open(file, 'r') as pdb:
            for line in pdb:
                if "ATOM" in line or "HETATM" in line:
                    tmpAtom = atom.Atom(line)
                    if chainLet != line[21:22]:
                        chainLet = line[21:22]
                        chains.append(chain.Chain(chainLet))

                    if resNum != int(line[22:26]):
                        resNum = int(line[22:26])
                        chains[-1].add_residue(residue.Residue(line))

                    chains[-1].residues[-1].add_atom(tmpAtom)
                lineNumber += 1
    except IOError:
        logger.exception(f"Error opening {file}")
        raise

    return protein.Protein(file, chains)


def make_mol2(res: residue, reformat: bool=True):
    # TODO: try and except wrap
    # TODO: add logger stuff to this
    successful = False
    with open(f"{res.name}.pdb", 'w') as mol2_file:
        for atom in res.atoms:
            mol2_file.write(atom.pdb_line())

        mol2_file.write('TER\nENDMDL')
    # Now we execute the babel command here
    with Popen(f"babel {res.name}.pdb {res.name}.mol2", stdin=PIPE, stdout=PIPE, stderr=PIPE,
               universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
        while shell.poll() is None:
            logger.debug(shell.stdout.readline().strip())
            output = shell.stderr.readline().strip()
            logger.debug(output)
            if "1 molecule converted" in output:
                successful = True

    if not successful:
        logger.error("Could not create {residue.name} mol2 file!")
        raise OSError("mol2_file")

    if reformat:
        # Now we fix the mol2 file
        no_atoms = []
        hatoms = []
        polarhatoms = []
        atom_section = False
        bond_section = False

        mol_file_lines = []

        with open(f"{res.name}.mol2") as mol_file:
            for line in mol_file:
                if "ATOM" in line:
                    atom_section = True
                    bond_section = False

                elif "BOND" in line:
                    atom_section = False
                    bond_section = True

                elif atom_section:
                    if line[47:49] == 'N.' or line[47:49] == 'O.':
                        no_atoms.append(int(line[:7]))

                    elif line[47:49] == 'H ':
                        hatoms.append(int(line[:7]))

                elif bond_section:
                    if int(line[6:12]) in no_atoms and int(line[12:18]) in hatoms:
                        polarhatoms.append(int(line[12:18]))

                    elif int(line[12:18]) in no_atoms and int(line[6:12]) in hatoms:
                        polarhatoms.append(int(line[6:12]))

                mol_file_lines.append(line)

        atom_section = False
        atomline = 0
        with open(f"{res.name}.mol2", 'w+') as mol_file:
            for line in mol_file_lines:
                if "ATOM" in line:
                    atom_section = True

                elif "BOND" in line:
                    atom_section = False

                if atom_section:
                    if atomline in hatoms and atomline not in polarhatoms:
                        line = line[:47] + 'Eh' + line[49:]

                    atomline += 1

                mol_file.write(line)

    os.remove(f"{res.name}.pdb")
    logger.info(f"Successfuly made: {res.name} mol2")

def create_config():
    # TODO verify the accuracy of the config file!
    logger.debug(f"Creating .dmdpy in the home directory: {os.path.expanduser('~')}")
    home = os.path.expanduser('~')

    try:
        os.mkdir(os.path.join(home, '.dmdpy'))

    except FileExistsError:
        logger.debug(".dmdpy directory already exists, continuing")

    logger.debug("Placing default dmdpy_config.json in the .dmdpy directory")
    shutil.copy(pkg_resources.resource_filename('dmdpy.resources', 'dmd_config.json'), os.path.join(home, '.dmdpy'))
    shutil.copy(pkg_resources.resource_filename('dmdpy.resources', 'logger_config.json'), os.path.join(home, '.dmdpy'))

    logger.info(f"Please ensure that {os.path.join(home, '.dmdpy')} has the correct values")

def make_start_file(parameters: dict, start_time: int =0):
    logger.debug("Making the Start File")
    try:
        with open("dmd_start", 'w') as dmdstart:
            dmdstart.write(f"THERMOSTAT     {parameters['Thermostat']}\n")
            dmdstart.write(f"T_NEW          {parameters['Initial Temperature']}\n")
            dmdstart.write(f"T_LIMIT        {parameters['Final Temperature']}\n")
            dmdstart.write(f"HEAT_X_C       {parameters['HEAT_X_C']}\n")
            dmdstart.write(f"RESTART_FILE   {parameters['Restart File']}\n")
            dmdstart.write(f"RESTART_DT     {parameters['dt']}\n")
            dmdstart.write(f"ECHO_FILE      {parameters['Echo File']}\n")
            dmdstart.write(f"ECHO_DT        {parameters['dt']}\n")
            dmdstart.write(f"MOVIE_FILE     {parameters['Movie File']}\n")
            dmdstart.write(f"START_TIME     {start_time}\n")
            dmdstart.write(f"MOVIE_DT       {parameters['dt']}\n")
            dmdstart.write(f"MAX_TIME       {parameters['Time'] + start_time}\n")

    except IOError:
        logger.exception("Error writing out dmd_start file")
        raise

    logger.debug("made the start file!")

def make_state_file(parameters: dict, pdbName):
    logger.debug("Calling complex.linux")
    try:
        # TODO: There is an issue here with complex.linux not actually running
        # It is coming from the mol2 of the substrate...interesting...
        # There is a Segmentation fault (core dumped) error that occurs, asking Jack if he knows what the issue is
        # Jack thinks it is the segfault mike wrote in his HACK ALERT section
        # I have emailed the Dohkyan group regarding it...its only for certain pdbs...
        with Popen(
                f"complex-1.linux -P {dmd_config['PATHS']['parameters']} -I {pdbName} -T topparam -D 200 -p param -s state -C inConstr -c outConstr",
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

                    # now we reorder based upon the attempt!
                    bonds.sort(key=lambda x: int(x[attempts]))

                with open(mol2, 'w+') as mf:
                    for line in save:
                        mf.write(line)

                    for bond in bonds:
                        mf.write(f"{bond[0]}\t{bond[1]}\t{bond[2]}\t{bond[3]}\n")

                with Popen(
                        f"complex-1.linux -P {dmd_config['PATHS']['parameters']} -I {pdbName} -T topparam -D 200 -p param -s state -C inConstr -c outConstr",
                        stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                        env=os.environ) as shell:
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

def make_movie(initial_pdb, movie_file, output_pdb):
    """

    :param initial_pdb: name of the initial pdb for the dmd run
    :param movie_file: name of the movie file created from dmd
    :param output_pdb: name of the output pdb that is generated from the movie file
    :return:
    """
    try:
        logger.debug("Creating movie file")
        with Popen(
                f"complex_M2P.linux {dmd_config['PATHS']['parameters']} {initial_pdb} topparam {movie_file} {output_pdb} inConstr",
                stdout=PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True, shell=True,
                env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())
    except OSError:
        logger.exception("Error calling complex_M2P.linux")
        raise