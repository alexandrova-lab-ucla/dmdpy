# This script holds all of the functions that involve running water placement and solvation software

def find_solv_shell(PATH, base_name, method):
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
				if len(split_line) > 21:
					run_solv_data[split_line[0] + split_line[1]] = (float(split_line[4]) / 100) # Store the buried % as a decimal
				elif len(split_line) > 8:
					if split_line[0] == 'Coupled' and split_line[1] == 'residues': # And find location to stop reading
						end_read = True

	return run_solv_data


def find_partners(all_titr_res, run_solv_data, partner_cutoff, method):
	if method == 'propka31_prerun':
		
		for titr_res in all_titr_res:
			for partner_res in all_titr_res: # A nested for loop to identify amino acids of the correct protonation states close to each other to act as partners in a protonation/deprotonation reaction, this should be fine given how short the titr_res list is

				if titr_res.cps[0] == '-' and partner_res.cps[0] == '+': # Only a residue in a '-' protonation state can take protons from a residue in a '+' protonation state and don't let a residue reference its own protons
					temp_all_atom_pairs_list = []
					nearby_prots = 0
					for unbound_atom in titr_res.uac: # To preserve the order
						temp_atom_pairs_by_unbound_atom_list = []
						for proton in partner_res.hc:
							temp_mag = 0.0
							for i in range(len(proton)): # Calculate the distance
								temp_mag += (proton[i] - unbound_atom[i])**2
							temp_mag = math.sqrt(temp_mag)
	
							if temp_mag <= partner_cutoff:
								nearby_prots += 1
								temp_atom_pairs_by_unbound_atom_list += [[unbound_atom, proton]]

						temp_all_atom_pairs_list += [temp_atom_pairs_by_unbound_atom_list]

					if nearby_prots > 0: # Only add entry to self.pc if there are qualifying proton-unbound atom pairs
						current_res_list = [partner_res] + [[temp_all_atom_pairs_list]]
						partner_list = [titr_res] + [[temp_all_atom_pairs_list]]

						titr_res.pc += [[current_res_list]]
						bound_atom_partner.pc += [[partner_list]]


