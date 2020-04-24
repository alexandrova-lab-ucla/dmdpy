#!/usr/bin/env python

import os
import sys

# Unfortunately, this is the most condensed I can make my titratable feature as it exists now. Will clean this up to work with dmdpy better later.

def run_titr_feature(updated_parameters, inp_pdb, pH, buried_cutoff, partner_dist):
    script_dir = os.path.dirname(os.path.realpath(__file__))

    #TODO use python import modules to fix this

    #TODO, dmdpy already has protein reformatter built in
    os.system(script_dir + '/formatpdb/formatpdb.sh -i ' + inp_pdb + ' -o ' + inp_pdb + ' -f Standard')

    #TODO pip install propka and call it that way...
    os.system('python3 ' + script_dir + '/propka31 ' + inp_pdb + ' > propka.stdout') # Run propka itself on the pdb file
    
    #TODO python modularize this..
    os.system(script_dir + '/main.py ' + inp_pdb + ' ./ prot.pdb ' + str(pH) + ' ' + str(buried_cutoff) + ' ' + str(partner_dist) + ' inConstr' + ' > titrConstr')

    # Old procedure bash calls:
    #predict-pka new.pdb out_all.pdb avepka curr_step
    #titration-inConstr new.pdb pH buried_cutoff partner_dist
    #prot2inConstr.py orignew.pdb prot.pdb new.pdb >> protsConstr

    new_prot_com = []
    with open('titrConstr', 'r') as titrCom_file:
        for line in titrCom_file:
            if len(line) > 3:
                temp_prot_com = []
                split_line_space = line.split()
                split_line_period = split_line_space[1].split('.')
                temp_prot_com = [chr(ord('A') - 1 + int(split_line_period[0])),int(split_line_period[1]),split_line_space[0].lower()]#, split_line_period[2]]
                new_prot_com += [temp_prot_com]

    updated_parameters['Custom protonation states'] = new_prot_com

    return updated_parameters

# For testing
'''
updated_parameters = 'stuff'
inp_pdb = 'initial.pdb'
pH = 3.0
buried_cutoff = 0.75
partner_dist = 3.5
run_titr_feature(updated_parameters, inp_pdb, pH, buried_cutoff, partner_dist)
'''
