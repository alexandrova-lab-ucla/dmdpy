

import sys
import Bio.PDB as bp
import xpdb # For the sloppy parser


def make_atom_naming_dict(atom_name_file, out_format, atom_type):
    out_names_dict = {}

    # Read in the names repository csvs
    name_file = open(atom_name_file, 'r+')
    name_file_lines = name_file.readlines()

    # Determine which column to read for output names
    labels = name_file_lines.pop(0)
    labels = labels.strip()
    labels = labels.split(',')
    out_num = labels.index(out_format)

    # For saving the N- and C- terminal atom names
    nterm_names = []
    cterm_names = []

    # Define all possible atom types by corresponding output name
    for line in name_file_lines:
        line = line.strip()
        line = line.split(',')
        res = line[0] + ':'
        if len(res) == 3:
            res = ' ' + res
        for i in range(len(line) - 1):
            out_names_dict[res + line[i + 1]] = line[out_num]

        if line[0] == 'NTERM': # And collect N- and C- terminal names
            nterm_names += line[1:]
        elif line[0] == 'CTERM':
            cterm_names += line[1:]

    # For titratable hydrogen this is identical, except it also saves the number of hydrogen to the dictionary under the standard res name
    if atom_type == 'titr':
        last_res = ''

        for line in name_file_lines:
            line = line.strip()
            line = line.split(',')
            res = line[0] + ':'
            for i in range(len(line) - 1):
                out_names_dict[res + line[i + 1]] = line[out_num]

            if line[0] != last_res:
                last_res = line[0]
                out_names_dict[last_res] = 1
            else:
                out_names_dict[last_res] += 1

            if line[0] == 'NTERM': # And collect N- and C- terminal names
                nterm_names += line[1:]
            elif line[0] == 'CTERM':
                cterm_names += line[1:]
        
    # Save the lists of N-terminal and C-terminal atom names
    out_names_dict['nterm_names'] = nterm_names
    out_names_dict['cterm_names'] = cterm_names

    return out_names_dict

def make_res_naming_dict(titr_res_name_file, out_format):
    res_name_dict = {} # Hold protonated / deprotonated name pairings for the given out_format
    res_name_standardize_dict = {} # Converts any residue name to the standard

    # Read in the names repository csvs
    titr_file = open(titr_res_name_file, 'r+')
    titr_file_lines = titr_file.readlines()

    # Determine which column to read for output names
    labels = titr_file_lines.pop(0)
    labels = labels.strip()
    labels = labels.split(',')
    out_num = labels.index(out_format)

    # In case there is a distinction between different protonation states which have the same number of hydrogen
    last_res_name = ''
    prot_state_num = 0

    # Define all possible atom types by corresponding output name
    for line in titr_file_lines:
        line = line.strip()
        line = line.split(',')

        if last_res_name != line[0]: # Construct the res name if there are multiple prot states
            last_res_name = line[0]
            prot_state_num = 0
        else:
            prot_state_num += 1
        if prot_state_num == 0:
            res_name_prot_state_num = ''
        else:
            res_name_prot_state_num = str(prot_state_num)

        res_name_dict[line[0] + res_name_prot_state_num + ':' + 'deprot'] = line[out_num] # Construct the res name pairing dictionary
        res_name_dict[line[0] + res_name_prot_state_num + ':' + 'prot'] = line[out_num + 1]

        for i in range(len(line) - 1): # Now construct standardizing dictionary
            res_name_standardize_dict[line[i + 1]] = line[0]

    return res_name_dict, res_name_standardize_dict

def find_res_bounds(chain, res_name_standardize_dict):
    bound_high_num = 'start'
    bound_low_num = 'start'

    for res in chain:
        if res.get_resname() in res_name_standardize_dict: # Check if it is a recognized residue to exclude ligands and solvent
            for atom in res:
                res_num = int(atom.get_full_id()[3][1])
                if bound_high_num == 'start':
                    bound_high_num = res_num
                    bound_low_num = res_num

                if bound_high_num < res_num:
                    bound_high_num = res_num
                elif bound_low_num > res_num:
                    bound_low_num = res_num

                break

    return bound_high_num, bound_low_num

def add_spacing(name): # Adds the appropriate amount of whitespace to a given atom name to create atom fullname
    if len(name) == 3:
        name = ' ' + name
    elif len(name) == 2:
        name = ' ' + name + ' '
    elif len(name) == 1:
        name = ' ' + name + '  '

    return name

def gen_struct(pdb_path):
    #parser = bp.PDBParser(PERMISSIVE=1)
    #structure = parser.get_structure('ID', pdb_path)

    # OR sloppy parser to deal with lots of water molecules
    sloppyparser = bp.PDBParser(PERMISSIVE=True, structure_builder=xpdb.SloppyStructureBuilder())
    structure = sloppyparser.get_structure('ID', pdb_path)

    return structure

def rename_pdb_atoms(structure, atom_name_dict, res_name_standardize_dict):
    for model in structure:
        for chain in model: # Strange algorithm here because BioPDB doesn't store the N- and C- termini, and doesn't store residue numbers in an easy to access spot
            cterm_num, nterm_num = find_res_bounds(chain, res_name_standardize_dict) # Identify N- and C- termini (always the first and last residues, but not always residue 1 and highest residue)
            nterm = False
            cterm = False

            for res in chain:
                res_name = res.get_resname()
                if res_name in res_name_standardize_dict:
                    res_name = res_name_standardize_dict[res.get_resname()] # Convert the res name to the standard to search the atom_name_dict properly
                res.resname = res_name

                for atom in res:
                    if nterm_num == int(atom.get_full_id()[3][1]):
                        nterm = True
                    elif cterm_num == int(atom.get_full_id()[3][1]):
                        cterm = True
                    else:
                        nterm = False
                        cterm = False
                    break

                # Now rename atoms, with process differing depending on whether residue is C- or N- terminus, or if it is neither
                if nterm == False and cterm == False:
                    for atom in res:
                        atom.name = atom_name_dict[res_name + ':' + atom.get_name()]
                        atom.fullname = add_spacing(atom_name_dict[res_name + ':' + atom.get_name()])
                elif nterm == True:
                    for atom in res:
                        if atom.get_name() in atom_name_dict['nterm_names']:
                            atom.name = atom_name_dict['NTERM:' + atom.get_name()]
                            atom.fullname = add_spacing(atom_name_dict['NTERM:' + atom.get_name()])
                        else:
                            atom.name = atom_name_dict[res_name + ':' + atom.get_name()]
                            atom.fullname = add_spacing(atom_name_dict[res_name + ':' + atom.get_name()])
                elif nterm == True:
                    for atom in res:
                        if atom.get_name() in atom_name_dict['cterm_names']:
                            atom.name = atom_name_dict['CTERM:' + atom.get_name()]
                            atom.fullname = add_spacing(atom_name_dict['CTERM:' + atom.get_name()])
                        else:
                            atom.name = atom_name_dict[res_name + ':' + atom.get_name()]
                            atom.fullname = add_spacing(atom_name_dict[res_name + ':' + atom.get_name()])

    return structure

def rename_pdb_res(structure, titr_atom_name_dict, titr_res_name_dict):
    for model in structure:
        for chain in model:
            for res in chain:
                res_name = res.get_resname()
                titr_atom_num = 0 # Count up the number of titratable hydrogen

                if res_name in titr_atom_name_dict:
                    if res_name != 'HIS':
                        for atom in res:
                            atom_ID = res_name + ':' + atom.get_name()
                            if atom_ID in titr_atom_name_dict:
                                titr_atom_num += 1

                        if titr_atom_name_dict[res_name] == titr_atom_num: # Next, determine and assign the protonation state based on this number
                            res.resname = titr_res_name_dict[res_name + ':' + 'prot']
                        else:
                            res.resname = titr_res_name_dict[res_name + ':' + 'deprot']
                    else: # Take special consideration for histidine (given the distinction between its E and D forms)
                        single_prot_form = ''
                        for atom in res:
                            D_form_ID = 'HISD' + ':' + atom.get_name()
                            E_form_ID = 'HISE' + ':' + atom.get_name()

                            if D_form_ID in titr_atom_name_dict:
                                single_prot_form = 'D'
                                titr_atom_num += 1
                            
                            elif E_form_ID in titr_atom_name_dict:
                                single_prot_form = 'E'
                                titr_atom_num += 1

                        if titr_atom_num == 2: # In this case it is doubly protonated
                            res.resname = titr_res_name_dict['HIS' + ':' + 'prot']
                        elif single_prot_form == 'D':
                            res.resname = titr_res_name_dict['HIS' + ':' + 'deprot']
                        elif single_prot_form == 'E':
                            res.resname = titr_res_name_dict['HIS1' + ':' + 'deprot']

    return structure


# User argument handling
in_pdb_path = sys.argv[1]
out_pdb_path = sys.argv[2]
name_file_path = sys.argv[3]
out_format = sys.argv[4]

atom_name_file = name_file_path + 'atom_names.csv'
titr_atom_name_file = name_file_path + 'titr_atom_names.csv'
titr_res_name_file = name_file_path + 'titr_res_names.csv'

# Load the atom and res naming conventions as dictionaries
atom_name_dict = make_atom_naming_dict(atom_name_file, out_format, 'all')
titr_atom_name_dict = make_atom_naming_dict(titr_atom_name_file, out_format, 'titr')
titr_res_name_dict, res_name_standardize_dict = make_res_naming_dict(titr_res_name_file, out_format)

# Make modifications to the pdb
structure = gen_struct(in_pdb_path)
structure = rename_pdb_atoms(structure, atom_name_dict, res_name_standardize_dict)
structure = rename_pdb_res(structure, titr_atom_name_dict, titr_res_name_dict)

# And print it
#io = bp.PDBIO()
io = xpdb.SloppyPDBIO() # Sloppy to deal with water molecules
io.set_structure(structure)
io.save(out_pdb_path)
