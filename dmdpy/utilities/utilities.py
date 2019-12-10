#!/usr/bin/env python3

import json
import pkg_resources
import logging
import os
from logging.config import dictConfig
from subprocess import Popen, PIPE


from dmdpy.protein import atom, chain, residue, protein


logger = logging.getLogger(__name__)


def load_logger_config():
    """Loads in the Logger Config File"""
    try:
        with open(pkg_resources.resource_filename('dmdpy.resources', 'logger_config.json')) as logger_config:
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

    logger.info(f"Successfuly made: {res.name} mol2")

