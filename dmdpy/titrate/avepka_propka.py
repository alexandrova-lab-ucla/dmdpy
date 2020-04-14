#!/usr/bin/env python

# This script calculate the propka pKas for each structure in the input trajectory file of pdbs and then averages the pKa values by residue

import sys
import os

def open_file_lines(file_path):
    in_file = open(file_path, 'r+')
    in_file_lines = in_file.readlines()
    in_file.close()

    return in_file_lines

def print_file(output, out_file_path):
    out_file = open(out_file_path, 'w')
    out_file.write(output)
    out_file.close()



def compile_pkas(pkas_list, propka_out):
    propka_lines = open_file_lines(propka_out)

    read_pkas = False
    for line in propka_lines:
        if 'Group      pKa' in line: # Only read the section labeled 'SUMMARY OF THIS PREDICTION'
            read_pkas = True
        elif read_pkas == True:
            if len(line) < 5:
                read_pkas = False
            else:
                split_line = line.split()
                curr_res = split_line[0] + ':' + split_line[1]
                if curr_res in pkas_list:
                    pkas_list[curr_res] += float(split_line[3])
                else:
                    pkas_list[curr_res] = float(split_line[3])

    return pkas_list

def update_propka_out(pkas_list, num_frames, propka_in_path, propka_out_path):
    propka_lines = open_file_lines(propka_in_path)

    read_pkas = False
    output = ''
    for line in propka_lines:
        if 'Group      pKa' in line:  
            read_pkas = True
            output += line
        elif read_pkas == True:
            if len(line) < 5:
                read_pkas = False
                output += line
            else:
                split_line = line.split()
                curr_res = split_line[0] + ':' + split_line[1]
                curr_pka = pkas_list[curr_res] / num_frames
                curr_pka_str = str(curr_pka)

                new_line = ''
                if curr_pka >= 10.0:
                    new_line = line[:16] + curr_pka_str[:5] + line[21:]
                else:
                    new_line = line[:17] + curr_pka_str[:4] + line[21:]

                output += new_line

        else:
            output += line

    print_file(output, propka_out_path)



def iter_movie_pdb_frames(movie_pdb_path, script_loc, propka_outname):
    movie_pdb_lines = open_file_lines(movie_pdb_path)

    pkas_list = {}
    frame_coords = ''
    num_frames = 0
    for line in movie_pdb_lines:
        if 'ENDMDL' in line or 'END' in line:
            # Process the last frame and add its pKas to the running tally
            print_file(frame_coords, 'avepka_temp.pdb')
            os.system(script_loc + '/../formatpdb/formatpdb.sh -i avepka_temp.pdb -o avepka_temp.pdb -f Standard')
            os.system(script_loc + '/propka31 avepka_temp.pdb > avepka_temp.stdout')
            pkas_list = compile_pkas(pkas_list, 'avepka_temp.stdout')
            frame_coords = ''
            num_frames += 1
        elif 'ATOM' in line or 'HETATM' in line:
            frame_coords += line

    update_propka_out(pkas_list, num_frames, 'avepka_temp.stdout', propka_outname)
    #os.system('rm avepka_temp.pdb avepka_temp.stdout')



pdb_trajectory = sys.argv[1]
propka_outname = sys.argv[2]
script_loc = sys.argv[3]

iter_movie_pdb_frames(pdb_trajectory, script_loc, propka_outname)
