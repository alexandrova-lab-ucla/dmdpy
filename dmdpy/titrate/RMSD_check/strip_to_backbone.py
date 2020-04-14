#!/usr/bin/env python

import sys

def strip_to_backbone(pdb_file, out_file):
    infile = open(pdb_file, 'r+')
    lines = infile.readlines()
    output = ''
    pdb_outfile = ''
    
    # List of backbone atom labels
    backbone_atom_labels = ['N', 'CA', 'C', 'O']

    for line in lines:
        line_entries = line.split()
        if len(line_entries) > 3:
            if line_entries[0] == 'ATOM':
                if line_entries[2] in backbone_atom_labels:
                    output += line

    output += 'TER\n'
    #output += 'ENDMDL\n'

    outfile = open(out_file, 'w+')
    outfile.write(output)


pdb_file = sys.argv[1]
out_file = sys.argv[2]

strip_to_backbone(pdb_file, out_file)
