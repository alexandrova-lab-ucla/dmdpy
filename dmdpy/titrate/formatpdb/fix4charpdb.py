import sys

def read_pdb_lines(pdb_path):
    pdb_file = open(pdb_path, 'r+')
    pdb_lines = pdb_file.readlines()

    return pdb_lines

def remove_4char_resname_whitespace(pdb_in_lines):
    pdb_out_lines = []

    for line in pdb_in_lines:
        if len(line) > 20:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                if line[20] != ' ':
                    pdb_out_lines += [line[0:21] + line[22:]]
                else:
                    pdb_out_lines += [line]
            else:
                pdb_out_lines += [line]
        else:
            pdb_out_lines += [line]

    return pdb_out_lines

def write_pdb_lines(pdb_lines, pdb_path):
    output = ''
    for line in pdb_lines:
        output += line

    pdb_file = open(pdb_path, 'w')
    pdb_file.write(output)
    pdb_file.close()


pdb_path = sys.argv[1]

pdb_lines = read_pdb_lines(pdb_path) # Read the pdb lines
fixed_pdb = remove_4char_resname_whitespace(pdb_lines) # Fix the pdb
write_pdb_lines(fixed_pdb, pdb_path) # Write the pdb from lines
