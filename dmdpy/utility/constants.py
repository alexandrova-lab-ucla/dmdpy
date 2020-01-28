#!/usr/bin/env python3

__all__=[
    'AMINO_ACID_RESIDUES',
    'METALS',
    'THERMOSTATS',
    'PROTONATED',
    'DEPROTONATED',
    'HEAVY_ATOMS'
]

AMINO_ACID_RESIDUES = ('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL')

THERMOSTATS = ('ANDERSON')

# Given a Residue, it can get the heavy atom and the proton that can be either protonated or deprotonated
PROTONATED = {
    "ASP": [("OD1", "2HND"), ("OD2", "1HND")],
    "GLU": [("OE1", "2HNE"), ("OE2", "1HNE")],
    "HIS": [("NE2", "1HNE")]
}

#Need to add arganine and lysine
DEPROTONATED = {
    "SER": [("OG", "HO")],
    "CYS": [("SG", "HG1")],
    "THR": [("OG1", "HO")],
    "ASN": [("ND2", "1HND"), ("ND2", "2HND")],
    "GLN": [("NE2", "1HNE"), ("NE2", "2HNE")],
    "TYR": [("OH", "HO")],
    "TRP": [("NE1", "HE1")],
    "HIS": [("ND1", "HD1")],
    "ARG": [("NH1", "2HH1"), ("NH1", "1HH1"), ("NH2", "2HH2"), ("NH2", "1HH2"), ("NE", "HE")],
    "LYS": [("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3")]
}

HEAVY_ATOMS = ['n', 'o', 's', 'se']

#All metals through bismuth
METALS = ('li','be','na','mg','al','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn','ga','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os','ir','pt','au','hg','tl','pb','bi')