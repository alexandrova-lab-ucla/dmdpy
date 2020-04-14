# This script holds all of the functions that involve running pKa estimation software

def yesimimporteddummy(message):
	print(message)

def calc_pKa_total_pdb(pdb_in_path, all_titr_res, method):
	if method == 'propka31': # Run propka3.1 on the input pdb and generate the dictionary of pKa by ionizable residue
		subprocess.call("propka31 " + pdb_in_path)
		
		propka_out_path = pdb_in_path[:-3] + '.pka' # Just change the extension for propka output
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
						calc_pKa_data[split_line[0] + split_line[1]] = float(split_line[3])

	return calc_pKa_data


