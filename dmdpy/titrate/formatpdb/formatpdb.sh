#!/bin/bash

# This program requires python 3 to be loaded with the biopython package

while getopts hi:o:f: option
do
case "${option}"
in
h) help="True";;
i) input=${OPTARG};;
o) output=${OPTARG};;
f) format=$OPTARG;;
*)
esac
done

if [ "$help" == "True" ]; then
    echo 'FormatPDB'
    echo ''
    echo 'This script converts between PDB atom naming schemes'
    echo 'Call this script anywhere'
    echo 'Outputs a formatted pdb file'
    echo ''
    echo 'NOTE: this script requires python3'
    echo 'If on Hoffman2 and this is not in your path, issue the following command:'
    echo '   module load python/3.6.1_shared'
    echo 'Other versions can be found with:'
    echo '   module avail python'
    echo ''
    echo 'FormatPDB.sh [-h] [-i input.pdb] [-o output.pdb] [-f format]'
    echo '-h: help - display this information'
    echo '-i: input - path to pdb to be formatted'
    echo '-o: output - path where formatted pdb should be stored'
    echo '-f: format - naming convention to convert to, currently supported:'
    echo '           - Standard (used for titr-MD)'
    echo '           - OPLSAA'
    echo '           - DMD'
    echo '           - Amber03'
    echo '           - GROMOS96-54a7'
    echo ''
elif [ "$input" == "?" ] || [ "$output" == "?" ] || [ "$format" == "?" ]; then
    echo 'Must provide proper -i, -o, -f options'
    echo 'Call FormatPDB.sh -h to see the run information'
else
    echo 'Call FormatPDB.sh -h to see the run information'

    script_loc=$(dirname $(readlink -f $0))

    name_file_loc=$script_loc/resources/

    python3 $script_loc/fix4charpdb.py $input # Allow for 4 character residues... because biopython PDBparser doesn't allow you to remove the extra whitespace when you make a new 4 character residue, but can somehow only read them when there is the proper amount of whitespace
    python3 $script_loc/main.py $input $output $name_file_loc $format # Create the new file
    python3 $script_loc/fix4charpdb.py $output # And clean up after

fi
