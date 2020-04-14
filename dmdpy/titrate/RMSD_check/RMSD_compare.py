#!/usr/bin/env python

import sys
import os

def strip_to_path(inp_arg):
    curr_dir = ''

    inp_arg_split = inp_arg.split('/')
    inp_arg_split.pop(-1)

    for directory in inp_arg_split:
        curr_dir += directory + '/'

    return curr_dir

def calc_RMSD(script_dir, in_pdb_path, ref_pdb_path):
    RMSD_val = 0

    os.system(script_dir + 'newRMSD.py ' + ref_pdb_path + ' ' + in_pdb_path + ' -f pdb >> CURR_RMSD_VAL')

    with open('CURR_RMSD_VAL') as f:
        lines = f.readlines()
        RMSD_val = lines[0][0:10]

    return RMSD_val

def print_RMSD(RMSD_VAL):
    with open('PREV_RMSD_VAL', 'w') as f:
        f.write(str(RMSD_VAL))

def compare_RMSD(NEW_RMSD_VAL, PREV_RMSD_VAL_path):
    with open(PREV_RMSD_VAL_path) as f:
        lines = f.readlines()
        PREV_RMSD_VAL = lines[0][0:10]

    RMSD_diff = abs(float(NEW_RMSD_VAL) - float(PREV_RMSD_VAL))

    if NEW_RMSD_VAL == PREV_RMSD_VAL:
        open('RMSD_CHECK_FAILED', 'a').close()
        is_failure = True
    elif RMSD_diff >= 1.5: # Assumes that a change in RMSD this big in just one step is also problem...
        open('RMSD_CHECK_FAILED', 'a').close()
        is_failure = True
    else:
        is_failure = False

    return is_failure

script_dir = strip_to_path(sys.argv[0])
in_pdb_path = sys.argv[1]
ref_pdb_path = sys.argv[2]

NEW_RMSD_VAL = calc_RMSD(script_dir, in_pdb_path, ref_pdb_path)

if len(sys.argv) == 3:
    print_RMSD(NEW_RMSD_VAL)
else:
    PREV_RMSD_VAL_path = sys.argv[3]
    is_failure = compare_RMSD(NEW_RMSD_VAL, PREV_RMSD_VAL_path)
    if is_failure == False:
        print_RMSD(NEW_RMSD_VAL)

os.system('rm CURR_RMSD_VAL')
