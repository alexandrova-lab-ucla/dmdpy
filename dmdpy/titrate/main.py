#!/usr/bin/env python

# This script initializes the class structure for titratable residues
# and performs all necessary operations on it making the appropriate
# calls to other scripts

# input: pdb file
# output: pdb file with updated protonation states

import tab_data as DATA
import subprocess
import math
import random
import sys
from itertools import combinations
import re

class titr_res:
	def __init__(self, amino_acid, res_num, chain, atoms, ter):
		self.name = amino_acid + res_num + chain
		self.amino_acid = amino_acid
		self.res_num = res_num
		self.chain = chain
		if ter == 'No':
			ter_res = amino_acid
		else:
			ter_res = ter
		self.ter_name = ter_res
		self.full_name = ter_res + res_num + chain # Like self.name but with 'C' or 'N' for C and N-terminal residues
		self.atoms = atoms # List of lists with atom type followed by cartesian coordinates for every atom in the residue
		self.prots = [] # List of lists with titratable protons, stored in same format as self.atoms
		self.heteroatoms = [] # List of lists with heteroatoms bound to the titratable protons, stored in the same format as self.atoms
		self.prot_state = [] # List with two entries, with the first a string + or - for protonated or deprotonated, and the second a number corresponding to its form
		self.partners = [] # List of titr_res instances, containing nearby amino acids that this amino acid can interact with
		self.pKa = 'none'
		self.new_prot_state = ['-', 0] # New protonation state decided by script, same formatting as self.prot_state
		self.change = ['None'] # List with two entries, the first saying whether there is a change to protonation state and what it is, the second as a list specifying the added or removed proton(s) if there is a change
		self.change_heteroatom = [] # List: If there is an added proton, name of the heteroatom it is bound to (1st entry) and a series of reference atoms for hydrogen placement (following entries)
		self.covalent_link = 'False'
		self.change_prob = 0.0
		self.change_roll = 0.0
		
	def define_prot_state(self):
		titr_atoms = []
		titr_het_atoms = []
		
		for atom in self.atoms: # Find protons and heteroatoms that define the titration state of the residue
			atom_name = self.ter_name + ':' + atom[0]
			if atom_name in DATA.titr_form_prots:
				self.prots += [atom]
				titr_atoms += [atom[0]]
			elif atom_name in DATA.titr_form_heteroatoms:
				self.heteroatoms += [atom]
				titr_het_atoms += [atom[0]]
				
		titr_atoms_ID = self.ter_name # Sort present titratable protons and heteroatoms to define protonation state
		titr_atoms.sort()
		titr_het_atoms.sort()
		for i in range(len(titr_atoms)):
			titr_atoms_ID += ':' + titr_atoms[i]
		for i in range(len(titr_het_atoms)):
			titr_atoms_ID += ':' + titr_het_atoms[i]

		#print(self.full_name)
		prot_state = DATA.prots2titr_form[titr_atoms_ID]
		prot_state = prot_state.split(':')
		self.prot_state = [prot_state[1], prot_state[2]]
		
	def update_ter(self, ter_res):
		self.ter = ter_res
		
	def assign_pKa(self, pKa_data, method): # Looks up pKa from dictionary pKa_data and saves this to class instance
		if method == 'propka31':
			if self.full_name in pKa_data:
				self.pKa = pKa_data[self.full_name]
			else: # If for some reason PropKa3.1 didn't calculate the pKa for this residue, set the pKa so that the protonation state can't change and notify the user
				#print('pKa not found for ' + self.full_name + ', skipping for now but the formatting for this residue should be checked')
				if self.prot_state[0] == '-':
					self.pKa = 1.0
				else:
					self.pKa = 20.0
					
			if self.covalent_link == 'True': # Do not allow residues involved in covalent interactions with other residues change protonation state
				if self.prot_state[0] == '-':
					self.pKa = 1.0
				else:
					self.pKa = 20.0
		
	def calc_pKa(self): # UNUSED: call external script to calculate pKa for single residue and save it to class instance
		# Currently this doesn't do anything
		self.pKa = 'None'

	def create_new_prot(self, new_atom_name):
		ref_coords = [0.0, 0.0, 0.0] # Need to generate cartesian coordinates for the new proton
		heteroatom_coords = [0.0, 0.0, 0.0]
		reference_atoms = self.change_heteroatom[1:]

		for atom in self.atoms: # Generate vector pointing from average of reference atoms to the heteroatom and project that from the heteroatom at a distance of 1 angstrom for the proton position
			if atom[0] == self.change_heteroatom[0]:
				heteroatom_coords[0] += atom[1]
				heteroatom_coords[1] += atom[2]
				heteroatom_coords[2] += atom[3]
			elif atom[0] in reference_atoms:
				ref_coords[0] += atom[1]/len(reference_atoms)
				ref_coords[1] += atom[2]/len(reference_atoms)
				ref_coords[2] += atom[3]/len(reference_atoms)
					
		prot_vec = [0.0, 0.0, 0.0]
		vec_length = 0.0
		for i in range(3):
			prot_vec[i] += heteroatom_coords[i] - ref_coords[i]
			vec_length += prot_vec[i]**2
		vec_length = math.sqrt(vec_length)
		
		new_prot_pos = [0.0, 0.0, 0.0]
		for i in range(3): # Place proton with vector projected from the heteroatom normalized to 1 angstrom
			new_prot_pos[i] += (prot_vec[i] / vec_length) + heteroatom_coords[i]
				
		self.atoms += [[new_atom_name, new_prot_pos[0], new_prot_pos[1], new_prot_pos[2]]]
		
	def update_prots(self):
		if self.change[0] == 'Remove':
			for i in range(len(self.atoms)):
				if self.atoms[i][0] == self.change[1][0]:
					self.atoms.pop(i)
					break
		elif self.change[0] == 'Add':
			for new_atom_name in self.change[1]:
				self.create_new_prot(new_atom_name)
		
# Titratable residue defining functions	
		
def assess_pdb_line_format(pdb_lines): # Gives parsing scripts the correct positions for things like atom coordinates as pdb formatting can vary
	pdb_line_atom_num = 1
	pdb_line_atom_type = 2
	pdb_line_res_type = 3
	pdb_line_res_num = 0
	pdb_line_x_coord = 0
	pdb_line_element = 0
	pdb_line_chain = 'none'
	for line in pdb_in_lines:
		if len(line) > 20:
			if line[0:4] == 'ATOM' or line[0:6] == 'HETATM': # Assess based on the first line about an atom in a titratable residue in the pdb
				split_line = line.split()
				if split_line[pdb_line_res_type] in DATA.titr_amino_acids_3let:
					if split_line[4] == 'A' or split_line[4] == 'B': # if the 5th column holds the chain rather than the residue number
						pdb_line_chain = 4
						pdb_line_res_num = 5
						pdb_line_x_coord = 6
						pdb_line_element = 11
					else:
						pdb_line_res_num = 4
						pdb_line_x_coord = 5
						pdb_line_element = 10
					current_res = split_line[pdb_line_res_type] + split_line[pdb_line_res_num]
					break
	
	return pdb_line_atom_num, pdb_line_atom_type, pdb_line_res_type, pdb_line_res_num, pdb_line_x_coord, pdb_line_element, pdb_line_chain		
		
def process_in_pdb(pdb_in_lines):
	all_titr_res = [] # List of titratable residues
	current_res = 0 # Place holders for iterating
	current_res_allatoms = []
	current_res_name = ''
	hit_ter = 'yes'
	ter_res = [] # List of terminal residue names
	ter_type = [] # And whether they are C or N terminal

	# Assess pdb line format
	pdb_line_atom_num, pdb_line_atom_type, pdb_line_res_type, pdb_line_res_num, pdb_line_x_coord, pdb_line_element, pdb_line_chain = assess_pdb_line_format(pdb_in_lines)

	# Find all titratable amino acids in pdb and initialize class instances
	for i in range(len(pdb_in_lines)):
		line = pdb_in_lines[i]
		if i != 0: # Define the previous line to identify N-terminal residues
			prev_line = pdb_in_lines[i - 1]
		
		if len(line) > 20:
			if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
				split_line = line.split()
				if split_line[pdb_line_res_type] in DATA.titr_amino_acids_3let:
					if split_line[pdb_line_res_num] != current_res:
						if current_res != 0:
							current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, 'No')
							current_titr_res.define_prot_state()
							all_titr_res += [current_titr_res]

						current_res = split_line[pdb_line_res_num] # Reset variables for the next titratable residue
						current_res_name = split_line[pdb_line_res_type]
						if pdb_line_chain != 'none':
							current_chain = split_line[pdb_line_chain]
						else:
							current_chain = ''
						current_res_allatoms = []
					current_res_allatoms += [[split_line[pdb_line_atom_type], float(split_line[pdb_line_x_coord]), float(split_line[pdb_line_x_coord + 1]), float(split_line[pdb_line_x_coord + 2])]]

				if hit_ter == 'yes':
					if split_line[pdb_line_res_type] in DATA.amino_acids_3let:
						if pdb_line_chain != 'none':
							ter_res += [split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain]]
						else:
							ter_res += [split_line[pdb_line_res_type] + split_line[pdb_line_res_num]]
						ter_type += ['N+']
						hit_ter = 'no'
					
		if len(line) >= 4:
			if line[0:3] == 'TER':
				split_line = prev_line.split()
				if split_line[pdb_line_res_type] in DATA.amino_acids_3let:
					if pdb_line_chain != 'none':
						ter_res += [split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain]]
					else:
						ter_res += [split_line[pdb_line_res_type] + split_line[pdb_line_res_num]]
					ter_type += ['C-']
					hit_ter = 'yes'
				
	# Enter the last titratable residue in
	current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, 'No')
	current_titr_res.define_prot_state()
	all_titr_res += [current_titr_res]
				
	current_res = 0
	
	# Loop back through the pdb to save the C and N terminal residues
	for line in pdb_in_lines:
		if len(line) > 20:
			if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
				split_line = line.split()
				if pdb_line_chain != 'none':
					current_res_full_name = split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain]
				else:
					current_res_full_name = split_line[pdb_line_res_type] + split_line[pdb_line_res_num]
					
				if current_res_full_name in ter_res:
					if split_line[pdb_line_res_num] != current_res:
						if current_res != 0:
							current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, ter_type[0])
							current_titr_res.define_prot_state()
							all_titr_res += [current_titr_res] # Store terminal residues that are titratable redundantly
							ter_res.pop(0) # Clear the 0th entry because they are found sequentially
							ter_type.pop(0)

						current_res = split_line[pdb_line_res_num] # Reset variables for the next titratable residue
						current_res_name = split_line[pdb_line_res_type]
						if pdb_line_chain != 'none':
							current_chain = split_line[pdb_line_chain]
						else:
							current_chain = ''
						current_res_allatoms = []
					current_res_allatoms += [[split_line[pdb_line_atom_type], float(split_line[pdb_line_x_coord]), float(split_line[pdb_line_x_coord + 1]), float(split_line[pdb_line_x_coord + 2])]]
		
	# Enter the last terminal residue in
	current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, ter_type[0])
	current_titr_res.define_prot_state()
	all_titr_res += [current_titr_res]
		
	return all_titr_res
	
# Amino acid interaction network definition
	
def distance(atom1, atom2):
	dist12 = (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2
	dist12 = math.sqrt(dist12)
	
	return dist12
	
def define_connections(all_titr_res, int_cutoff):
	all_heteroatoms = []
	
	for res in all_titr_res: # Construct list of heteroatoms
		for heteroatom in res.heteroatoms:
			all_heteroatoms += [[res, heteroatom]]
			
	# Sort list of heteroatoms by x position
	all_heteroatoms.sort(key=lambda x: x[1][1])
	
	for i in range(len(all_heteroatoms) - 1):
		j = 1

		while (all_heteroatoms[i + j][1][1] - all_heteroatoms[i][1][1]) <= int_cutoff: # If the x axis distance between residues is less than the cutoff, there could be an interaction
			if all_heteroatoms[i][0].full_name != all_heteroatoms[i + j][0].full_name: # Only look for interactions between different residues
				if all_heteroatoms[i + j][0] not in all_heteroatoms[i][0].partners: # No need to check for more interactions if there already is one between two residues
					if distance(all_heteroatoms[i + j][1], all_heteroatoms[i][1]) <= int_cutoff: # If the distance between the two atoms is within the cutoff, register as an interaction
						all_heteroatoms[i][0].partners += [all_heteroatoms[i + j][0]]
						all_heteroatoms[i + j][0].partners += [all_heteroatoms[i][0]]
						if all_heteroatoms[i][0].amino_acid + ':' + all_heteroatoms[i][1][0] in DATA.cov_link_heteroatoms and all_heteroatoms[i + j][0].amino_acid + ':' + all_heteroatoms[i + j][1][0] in DATA.cov_link_heteroatoms: # Check if the connection is a covalent link, like a disulfide bridge
							all_heteroatoms[i][0].covalent_link = 'True' # If this is the case, this amino acid is involved in a covalent link and its probability of protonation state change is greatly reduced
							all_heteroatoms[i + j][0].covalent_link = 'True'
			if (i + j) == (len(all_heteroatoms) - 1):
				break
			else:
				j += 1

	
def define_aa_networks(all_titr_res):
	all_networks = []
	current_network = []
	already_in_network = []
	already_in_current_network = []
	
	while all_titr_res != []: # Empty stack of all residues to trace all networks
		if all_titr_res[0] in already_in_network: # Don't consider amino acids already part of a network
			all_titr_res.pop(0)
		else:
			current_network += [all_titr_res[0]]
			current_partners = all_titr_res[0].partners
			already_in_current_network += [all_titr_res[0]]
			all_titr_res.pop(0)
		
		while current_partners != []: # If there are partners to the current residue, it is part of a large network
			if current_partners[0] in already_in_current_network: # Takes care of loops, when multiple residues interact with eachother it add duplicates to the current_partner list
				current_partners.pop(0)
			else:
				current_network += [current_partners[0]]
				current_partners += current_partners[0].partners # Trace out any chain of partners
				already_in_network += [current_partners[0]]
				already_in_current_network += [current_partners[0]]
				current_partners.pop(0) # This partner has now been checked
		
		if current_network != []:
			all_networks += [current_network]
		current_network = []
		already_in_current_network = []
		
	return all_networks
	
# pKa calculation scripts (because python refuses to acknowledge the existence of functions in imported files even though modules work fine)

def calc_pKa_total_pdb(pdb_in_path, all_titr_res, method, chains):
	if method == 'propka31': # Run propka3.1 on the input pdb and generate the dictionary of pKa by ionizable residue
		#subprocess.call("propka31 " + pdb_in_path)
		
		propka_out_path = pdb_in_path[:-3] + 'pka' # Just change the extension for propka output
		with open(propka_out_path, 'r+') as propka_out: # Open propka output
			propka_in_lines = propka_out.readlines()
			propka_out.close()

		calc_pKa_data = {} # Dictionary that stores pKa by residue name consistent with the established titr_res class
		
		start_read = False
		end_read = False
		for line in propka_in_lines:
			split_line = line.split()
			if start_read == False: # Finds the location in the file to start reading pKa values
				if len(split_line) > 2:
					if split_line[0] == 'Group' and split_line[1] == 'pKa':
						start_read = True
			elif start_read == True and end_read == False:
				if len(split_line) > 4:
					if split_line[0] == 'Free' and split_line[1] == 'energy': # And find location to stop reading
						end_read = True
					else:
						if chains == 'True':
							calc_pKa_data[split_line[0] + split_line[1] + split_line[2]] = float(split_line[3])
						else:
							calc_pKa_data[split_line[0] + split_line[1]] = float(split_line[3])

	return calc_pKa_data
	
# Solvation calculation scripts (because python refuses to acknowledge the existence of functions in imported files even though modules work fine)

def find_solv_shell(PATH, base_name, method, chains):
	if method == 'placevent':
		# Unfortunately placevent requires an AMBER .dx external potential file as input, so we must generate that first

		# Placevent initialization settings
		water_conc = 55.5 # in Molarity
		NaCl_conc = 0.005 # in Molarity

		dx_path = PATH + base_name + '.O.1' + '.dx' # Path to future dx file
		solv_pdb_path = PATH + base_name + '_solv.pdb' # Path and name of output solvated pdb

		subprocess.call("pdb_generate_dx.sh " + PATH + " " + base_name + " " + str(water_conc) + " " + str(NaCl_conc)) # Bash script which generates the .dx file, based on tutorial found at dansindhikara.com/Tutorials
		subprocess.call("placevent.py " + dx_path + " " + str(water_conc) + "> " + solv_pdb_path) # Now we run placevent itself
	elif method == 'propka31_prerun': # The simplest method (probably too much so), assumes that the solvent accessibility can be approximated as % buried data calculated by propka
		# The prerun suffix means that this method assumes that propka31 has already been run as the pKa prediction method
		propka_out_path = PATH + base_name + '.pka'
		with open(propka_out_path, 'r+') as propka_out: # Open propka output
			propka_in_lines = propka_out.readlines()
			propka_out.close()

		run_solv_data = {} # Dictionary that stores pKa by residue name consistent with the established titr_res class
		
		start_read = False
		end_read = False
		for line in propka_in_lines:
			split_line = line.split()
			if start_read == False: # Finds the location in the file to start reading pKa values
				if len(split_line) > 9:
					if split_line[0] == 'RESIDUE' and split_line[1] == 'pKa':
						start_read = True
			elif start_read == True and end_read == False:
				if len(split_line) > 15:
					if split_line[5] == '%':
						if chains == 'True':
							run_solv_data[split_line[0] + split_line[1] + split_line[2]] = (float(split_line[4]) / 100) # Store the buried % as a decimal
						else:
							run_solv_data[split_line[0] + split_line[1]] = (float(split_line[4]) / 100) # Store the buried % as a decimal
				elif len(split_line) > 8:
					if split_line[0] == 'Coupled' and split_line[1] == 'residues': # And find location to stop reading
						end_read = True

	return run_solv_data
	
# Amino acid network solvent accessibility definition

def find_network_solvent_access(old_networks, run_solv_data, solv_cutoff, solv_prob, method):
	new_networks = []
	
	for network in old_networks:
		solvated = 'False'
		atom = 0
		
		while solvated == 'False' and atom < len(network): # If any residue is within the solvation cutoff, save network as solvent accessible
			if network[atom].full_name in run_solv_data:
				if run_solv_data[network[atom].full_name] < solv_cutoff:
					solvated = 'True'	
			atom += 1
			
		if solvated == 'True':
			new_networks += [['Solv', network]]
		else:
			new_networks += [['NoSolv', network]]
			
	return new_networks
	
# New protonation state generation scripts

def MC_decide_solv(residue, resevoir_pH):
	prob_add = (10.0**(residue.pKa - resevoir_pH)) / (1.0 + 10.0**(residue.pKa - resevoir_pH))
	MC_prot_state_roll = float(random.randint(0, 1000000) / 1000000.0)

	if MC_prot_state_roll <= prob_add:
		if residue.prot_state[0] == '-': # If this is a change to protonation state, take proper action
			possible_prot_states = DATA.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Add' + ':' + str(residue.prot_state[1])]
			#print(possible_prot_states)
			MC_prot_form_roll = random.randint(1, len(possible_prot_states))
			residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
			residue.change = ['Add', possible_prot_states[MC_prot_form_roll - 1][1]]
			residue.change_heteroatom = DATA.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(residue.prot_state[1])] # Selects the first hydrogen if two are added (which only affects N-terminus)
			#print(residue.change_heteroatom)
	else:
		if residue.prot_state[0] == '+': # If this is a change to protonation state, take proper action
			possible_prot_states = DATA.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Remove' + ':' + str(residue.prot_state[1])]
			MC_prot_form_roll = random.randint(1, len(possible_prot_states))
			residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
			residue.change = ['Remove', possible_prot_states[MC_prot_form_roll - 1][1]]
			residue.change_heteroatom = DATA.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(possible_prot_states[MC_prot_form_roll - 1][0][1])] # Selects the first hydrogen if two are added (which only affects N-terminus)

			
	residue.change_prob = prob_add
	residue.change_roll = MC_prot_state_roll

def MC_decide_nosolv(network):
	available_prots = 0
	
	for residue in network: # Tally the number of mobile protons in the system
		if residue.prot_state[0] == '+':
			available_prots += 1
		
	prot_combos = combinations(network, available_prots) # Enumerate all possible total protonation states of the system
	probabilities = []
	partition_func = 0.0
	for combo in prot_combos: # Construct the partition function for the system
		temp_prob = 0.0
		for residue in combo:
			temp_prob += residue.pKa
		temp_prob = 10.0**(temp_prob)
		probabilities += [[temp_prob, combo]]
		partition_func += temp_prob
	
	MC_prot_state_roll = float(random.randint(0, 1000000) / 1000000.0) # Roll the dice and decide the state
	current_prob_total = 0.0
	unchanged_protonated_residues = []
	for state in probabilities:
		previous_prob_total = current_prob_total
		current_prob_total += state[0] / partition_func
		if MC_prot_state_roll < current_prob_total and MC_prot_state_roll > previous_prob_total: # If the MC roll falls within the current state probability range, choose this state
			for residue in state[1]:
				if residue.prot_state[0] == '-': # If this involves a protonation state change to this residue, update it
					possible_prot_states = DATA.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Add' + ':' + str(residue.prot_state[1])]
					MC_prot_form_roll = random.randint(1, len(possible_prot_states))
					residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
					residue.change = ['Add', possible_prot_states[MC_prot_form_roll - 1][1]]
					residue.change_heteroatom = DATA.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(residue.prot_state[1])] # Selects the first hydrogen if two are added (which only affects N-terminus)
				if residue.prot_state[0] == '+': # If this is an unchanged residue, record that
					unchanged_protonated_residues += [residue]
				residue.change_prob = current_prob_total
				residue.change_roll = MC_prot_state_roll
				
	for residue in network: # Go back around and update the residues losing protons
		if residue.prot_state[0] == '+' and residue not in unchanged_protonated_residues: # Remove protons from newly deprotonated residues
			possible_prot_states = DATA.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Remove' + ':' + str(residue.prot_state[1])]
			MC_prot_form_roll = random.randint(1, len(possible_prot_states))
			residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
			residue.change = ['Remove', possible_prot_states[MC_prot_form_roll - 1][1]]
			residue.change_heteroatom = DATA.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(possible_prot_states[MC_prot_form_roll - 1][0][1])] # Selects the first hydrogen if two are added (which only affects N-terminus)
			
def MC_prot_change(networks, solution_pH):
	for network in networks:
		if network[0] == 'Solv': # Solvent accessibility means any residue can be protonated
			for residue in network[1]:
				MC_decide_solv(residue, solution_pH)
		else: # Otherwise the number of available protons is fixed 
			MC_decide_nosolv(network[1])
			
# Outputted pdb processing functions

def fix_atom_num_whitespace(atom_line, pdb_line_atom_num, atom_tally):
	old_atom_num = atom_line[2 * pdb_line_atom_num]
	new_atom_num = str(atom_tally)
	whitespace = atom_line[2 * pdb_line_atom_num - 1]

	if len(old_atom_num) != len(new_atom_num):
		if len(old_atom_num) > len(new_atom_num): # Add a whitespace character if the new atom number is shorter than the old one
			whitespace = whitespace + ' '
		else: # Remove a whitespace character if the new atom number is longer than the old one
			whitespace = whitespace[:-1]
	
	return atom_line
	
def new_atom_format(atom_line, residue, atom_name, pdb_line_atom_type, pdb_line_x_coord, pdb_line_element):
	atom_line[2 * pdb_line_atom_type] = atom_name # Atom name
	atom_line[2 * pdb_line_element] = 'H' # Element is hydrogen
	hydrogen_name_length = len(atom_line[2 * pdb_line_atom_type]) # Correct whitespace around atom name
	if hydrogen_name_length == 4:
		atom_line[2 * pdb_line_atom_type - 1] = ' '
		atom_line[2 * pdb_line_atom_type + 1] = ' '
	elif hydrogen_name_length == 3:
		atom_line[2 * pdb_line_atom_type - 1] = '  '
		atom_line[2 * pdb_line_atom_type + 1] = ' '
	elif hydrogen_name_length == 2:
		atom_line[2 * pdb_line_atom_type - 1] = '  '
		atom_line[2 * pdb_line_atom_type + 1] = '  '
	elif hydrogen_name_length == 1:
		atom_line[2 * pdb_line_atom_type - 1] = '  '
		atom_line[2 * pdb_line_atom_type + 1] = '   '
	
	for i in range(3): # Atom coordinates
		crop_length = 0 # Determine number of characters appropriate for this coordinate and correct whitespace
		if residue.atoms[-1][i + 1] >= 100.0:
			crop_length = 7
			atom_line[2 * (pdb_line_x_coord + i) - 1] = ' '
		elif residue.atoms[-1][i + 1] >= 10.0:
			crop_length = 6
			atom_line[2 * (pdb_line_x_coord + i) - 1] = '  '
		else:
			crop_length = 5
			atom_line[2 * (pdb_line_x_coord + i) - 1] = '   '
		if i == 0: # The x coordinate needs more whitespace before it
			atom_line[2 * (pdb_line_x_coord + i) - 1] += '    '
		
		atom_line[2 * (pdb_line_x_coord + i)] = str(residue.atoms[-1][i + 1])[:crop_length] # Because the last atom is the added proton, and there can only be one at most
	
	return atom_line

def generate_out_pdb(pdb_in_lines, all_titr_res):
	output = ''
	new_atom_tally = 1
	first_res = 'Start'
	current_res_num = ''
	atoms_changing = 'False'
	
	# assess pdb line format
	pdb_line_atom_num, pdb_line_atom_type, pdb_line_res_type, pdb_line_res_num, pdb_line_x_coord, pdb_line_element, pdb_line_chain = assess_pdb_line_format(pdb_in_lines)

	for line in pdb_in_lines:
		if len(line) > 20:
			if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
				split_line = re.split(r'(\s+)', line)
				if split_line[2 * pdb_line_res_type] in DATA.titr_amino_acids_3let or first_res == 'Start' or first_res == 'True': # first_res is for the N-terminus						
					if split_line[2 * pdb_line_res_num] != current_res_num:
						if first_res != 'False': # Obnoxious, but the N-terminus is the first residue, but not always consistently numbered so need something like this to account for that
							if first_res == 'Start':
								first_res = 'True'
							elif first_res == 'True':
								first_res = 'False'
						current_res_num = split_line[2 * pdb_line_res_num]
						if pdb_line_chain != 'none':
							current_res_num += split_line[2 * pdb_line_chain]
						atoms_changing = 'False'
						for residue in all_titr_res:
							if pdb_line_chain != 'none':
								if current_res_num == str(residue.res_num) + residue.chain and residue.change[0] != 'None':
									atoms_changing = 'True'
									atoms_changing_res = residue
									break
							else:
								if current_res_num == str(residue.res_num) and residue.change[0] != 'None':
									atoms_changing = 'True'
									atoms_changing_res = residue
									break
						
						if atoms_changing == 'True':
							if split_line[2 * pdb_line_atom_type] in atoms_changing_res.change[1]: # Remove the lost proton
								pass
							elif atoms_changing_res.change[0] == 'Add':
								if split_line[2 * pdb_line_atom_type] == atoms_changing_res.change_heteroatom[0]: # Add any new proton(s) after its corresponding heteroatom
									# First take care of the heteroatom
									split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
									new_atom_tally += 1
									for entry in split_line:
										output += entry
								
									# Now the new proton(s)
									# Unfortunately, for the N-terminus this procedure adds both new hydrogen on top of each other, but this doesn't matter for our DMD protocols (as it doesn't treat hydrogen atomistically)
									for new_prot_name in atoms_changing_res.change[1]:
										split_line = new_atom_format(split_line, atoms_changing_res, new_prot_name, pdb_line_atom_type, pdb_line_x_coord, pdb_line_element)
									        # Fix surrounding whitespace for new entries
										split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
										new_atom_tally += 1
										for entry in split_line:
											output += entry
								
								else:
									split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
									new_atom_tally += 1
									for entry in split_line:
										output += entry
								
							else:
								split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
								new_atom_tally += 1
								for entry in split_line:
									output += entry
								
						else:
							split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
							new_atom_tally += 1
							for entry in split_line:
								output += entry

						current_res_num = split_line[2 * pdb_line_res_num]
							
					else:
						if atoms_changing == 'True':
							if split_line[2 * pdb_line_atom_type] == atoms_changing_res.change[1][0]: # Remove the lost proton
								pass
							elif atoms_changing_res.change[0] == 'Add':
								if split_line[2 * pdb_line_atom_type] == atoms_changing_res.change_heteroatom[0]: # Add any new proton(s) after its corresponding heteroatom
									# First take care of the heteroatom
									split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
									new_atom_tally += 1
									for entry in split_line:
										output += entry
								
									for new_prot_name in atoms_changing_res.change[1]: # Now the new proton(s)
										split_line = new_atom_format(split_line, atoms_changing_res, new_prot_name, pdb_line_atom_type, pdb_line_x_coord, pdb_line_element)
									# Fix surrounding whitespace for new entries
										split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
										new_atom_tally += 1
										for entry in split_line:
											output += entry
								
								else:
									if first_res != 'False' and split_line[2 * pdb_line_atom_type] == 'H': # Change the naming of the N-terminal 'H' to 'HN1' to make it consistent with the naming scheme
										split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
										new_atom_tally += 1
										
										split_line[2 * pdb_line_atom_type] = 'HN1'
										split_line[2 * pdb_line_atom_type + 1] = ' '

										for entry in split_line:
											output += entry

									else:
										split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
										new_atom_tally += 1
										for entry in split_line:
											output += entry
								
							else:
								split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
								new_atom_tally += 1
								for entry in split_line:
									output += entry

						else:
							split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
							new_atom_tally += 1
							for entry in split_line:
								output += entry
							
				else:
					split_line = fix_atom_num_whitespace(split_line, pdb_line_atom_num, new_atom_tally)
					new_atom_tally += 1
					for entry in split_line:
						output += entry
			else:
				output += line
		else:
			output += line
			
	return output



# Developer switches
calc_pKa_type = 'full' # 'full' or 'byres': specifies whether to use a pKa prediction algorithm that runs on the full protein or on a single residue
pKa_calc_style = 'propka31' # 'propka31': specifies which pKa estimation program to use
solv_shell_style = 'propka31_prerun' # 'placevent', 'propka31_prerun': specifies which solvation program to use, 'propka31_prerun' assumes propka31 as the pKa_calc_style
propka_buried_cutoff = 0.85 # Propka buried percentage cutoff to consider a residue solvent inaccessible
proton_partner_cutoff = 3.5 # Distance in angstrom within which another molecule could possibly take a proton from another

if __name__ == "__main__":
	PATH = sys.argv[2] # Save path to working directory
	pdb_in_path = PATH + sys.argv[1] # Path to input pdb
	pdb_out_path = PATH + sys.argv[3] # Path to output pdb
	
	solution_pH = float(sys.argv[4]) # pH of the solution can vary with some calculations
	propka_buried_cutoff = float(sys.argv[5]) 
	propka_solv_access_prob = float(sys.argv[6])

	#output_control = sys.argv[4] # Controls a number of features, including whether residue names change with change in protonation state

	with open(pdb_in_path, 'r+') as pdb_in: # Open pdb for input
		pdb_in_lines = pdb_in.readlines()
		pdb_in.close()

	all_titr_res = process_in_pdb(pdb_in_lines) # Identify and save titratable residues
	define_connections(all_titr_res, proton_partner_cutoff) # Define connected amino acids
	if all_titr_res[0].chain == '': # Determine if there are multiple chains for the pKa and solvation scripts
		chains = 'False'
	else:
		chains = 'True'
	
	import run_pka # Calculate the pKa of each residue
	if calc_pKa_type == 'full': # If the pKa prediction program runs on the entire pdb
		calc_pKa_data = calc_pKa_total_pdb(pdb_in_path, all_titr_res, pKa_calc_style, chains)
		for titr_res in all_titr_res:
			titr_res.assign_pKa(calc_pKa_data, pKa_calc_style)
	elif calc_pKa_type == 'byres': # UNUSED: If the pKa prediction program runs on each residue individually
		for titr_res in all_titr_res:
			titr_res.calc_pKa()
			
	import run_solvation # Locate water positions or water accessibility of each residue #fix
	run_solv_data = find_solv_shell(PATH, sys.argv[1][:-4], solv_shell_style, chains)
	if solv_shell_style == 'propka31_prerun':
		titr_stack = [] # Construct the stack form of all_titr_res for use in find_solv_shell
		for res in all_titr_res:
			titr_stack += [res]
		all_networks = define_aa_networks(titr_stack) # Define amino acid networks
		all_networks = find_network_solvent_access(all_networks, run_solv_data, propka_buried_cutoff, propka_solv_access_prob, solv_shell_style) # Identify solvent accessible networks
	
	MC_prot_change(all_networks, solution_pH) # Decide the new protonation states by a Metropolis Monte Carlo algorithm: rolling dice
	for residue in all_titr_res: # Update the atom list in each residue
		residue.update_prots()
		
	output = generate_out_pdb(pdb_in_lines, all_titr_res) # Create output pdb
	with open(pdb_out_path, 'w+') as pdb_out: # And print it
		pdb_out.write(output)
		pdb_out.close()
	
	
	
	### Data Output ###
	
	if len(sys.argv) == 8:
		data_output_style = sys.argv[7]
	else:
		data_output_style = 'None'

	if data_output_style == 'None':
		print_networks = 'False'
		print_changes = 'False'
		print_residues = 'False'
	elif data_output_style == 'inConstr':
		print_changes = 'True'
		print_networks = 'False'
		print_residues = 'False'
	elif data_output_style == 'Test-Res':
		print_changes = 'False'
		print_networks = 'False'
		print_residues = 'True'
	
	if print_networks == 'True':
		for network in all_networks:
			print('\nNetwork:')
			print(network[0])
			for residue in network[1]:
				if residue.change[0] != 'None':
					print(residue.full_name + ' - ' + str(residue.pKa) + ' - ' + str(residue.change_prob) + ' - ' + str(residue.change_roll) + ' - ' + residue.change[0] + ' - ' + residue.change[1][0])
				else:
					print(residue.full_name + ' - ' + str(residue.pKa) + ' - ' + str(residue.change_prob) + ' - ' + str(residue.change_roll) + ' - ' + residue.change[0])

	if print_residues == 'True':
		for res in all_titr_res:
			print(res.full_name)
			print(res.prot_state)
			print(res.heteroatoms)
			print(res.change)
			print(res.change_prob)
			print(res.change_roll)
			print(res.atoms)

	if print_changes == 'True':
		output = ''
		for res in all_titr_res:
			output_line = ''
			# Output lines written in the format (De)protonate chain_#.res_#.modified_atom

			atom_change = res.change
			if atom_change[0] != 'None':
				if atom_change[0] == 'Add':
					output_line += 'Protonate '
					is_change = 'True'
				elif atom_change[0] == 'Remove':
					output_line += 'Deprotonate '
					is_change = 'True'

				if is_change == 'True':
					output_line += str(DATA.alphabet_order[res.chain]) + '.'
					output_line += str(res.res_num) + '.'
					output_line += res.change_heteroatom[0] + '\n'
					
					output += output_line
		print(output)

				
