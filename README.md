dmdpy
=
A python interface for [piDMD](http://www.turbomole.com/). Currently supports up the latest version of piDMD.

## Table of Contents

* [Description](#Description)
* [Technologies](#Technologies)
* [Installation](#Installation)
* [Define Options](#options)
* [Submission Keywords](#Submission-Keywords)
* [Usage](#Usage)
* [Contributing](#Contributing)
* [Authors](#Authors)
* [License](#License)

## Description
This project can be used as a set of executables or as python modules. The three main executable scripts are `setupturbomole.py`, `runturbomole.py`, and `submitturbomole.py`. `setupturbomole.py` sets up a turbomole job, `runturbomole.py` runs a turbomole job, and `submitturbomole.py` submits a turbomole job to a queueing system.

The modules `turbopy.calculation` and `turbopy.setupjob` can be imported to other python projects be accessed directly.

## Technologies
Project created with:
* [Python](https://www.python.org/downloads/) version: 3.6
* [TURBOMOLE](http://www.turbomole.com/) version: 6.6

## Installation
This code needs a python version greater than 3.7. Check your python version using
 ```bash
 python -V
 ```
 This code is still in development, so it is best to clone or download this git repository. Then you can install it in using pip! The `-e` allows for continual development.
```bash
git clone https://github.com/alexandrova-lab-ucla/dmdpy.git
cd turbopy
pip install --user -e .
```
To ensure that the executables are in your PATH, make sure to extend your path to include the `~/.local/bin` directory!

To use the executables, the codes program needs to know specific information from you about where piDMD is installed and your queuing system.
I would recommend running `setupdmd.py` as this will initialize dmdpy by creating `.phd` in your home directory. When you first execute this,
you will see a bunch of errors and warnings but this is normal. In the `~/.phd` directory, please edit the `phd_config.json` file to update 
the `DMD_DIR` to the absolute path to your installation of the piDMD binaries and the `parameters` to the absolute path to the directory with 
the piDMD parameters. If you will not be using the `submitdmd.py` script, you can ignore the `QUEUEING` section; however, if you are using it 
please update the value for each of the keys with the following information:

###### 1. max_nodes: 
Maximum number of nodes you want a job to be able to run on.
###### 2. node_types: 
The available node types for the system you are using with the entires being the number of cores available.
###### 3. max_time: 
Maximum time (hours) allowed for each of the node type listed above. Keep them in the same order.
###### 4. high_priority: 
Nodes that you have high priority access to.
###### 5. submit: 
the key word to submit a job to the queuing system.

After the `phd_config.json` file is updated, try rerunning `setupdmd.py`. This time, it should copy over a blank `dmdinput.json`. If it did this, 
then its installation was successful! If you are going to be using `submitdmd.py` then you also need to supply it with a submission template in 
the form jinja2 format. You can view a sample `submit.j2` in `/dmdpy/templates/submit.j2` for a UGE queuing system. The `submit.j2` is highly 
platform and system dependent. A list of keywords is listed in [Submission Keywords](#Submission-Keywords) section.

## DMD Options <a name="options"></a>
For a full list of what each of these parameters does, it is best to refer to the piDMD manual. This file is meant to make life easier to replicate
similar DMD simulations across various proteins. 
##### Thermostat
Theremostat to use. Default is `Anderson`

##### Initial Temperature
Initial temperature to start each DMD run at 

##### Final Temperature
Final temperature to end each DMD run at. During the course of each DMD simulation, DMD will gradually change the temperature from the initial 
to final temperature, allowing a gradient. If you are stringing multiple DMD simulations together (through the `commands` options), I would not 
recommend changing the temperature here (globally), but instead inside each of the command.

##### HEAT\_X\_C
Heat capacity as defined in the piDMD manual.

##### Echo File
File name to save all of the piDMD data to. Refer to the piDMD manual for explanation of the echo file.

##### Movie File
File name to save the DMD trajectory to. Refer to the piDMD manual for explanation of the movie file.

##### Restart File
File name to save the restart information to. This is the file that is used instead of the state file to restart a DMD simulation from.

##### dt
Interval of time steps in which to save data (Echo, Movie, Restart). 

##### Time
Total number of steps you wish each simulation to go for.

##### titr
Work in progress. Will allow for titratable DMD simulations in the future.

##### Freeze Non-Residues
This will hold static any residue, or groups of atoms, that are not part of any canonical amino acids. This is because DMD is not particularly 
good at modeling anything that isn't a normal amino acid or zinc.

##### Restrict Metal Ligands
This will search for any atoms that are within a certain cutoff distance of every metal and hold these atoms also static during the DMD simulation.
Additionally, all atoms bound to these now frozen metal ligand atom, will have their displacement restricted to 0.2 to prevent any explosions 
during the simulation.

##### Custom protonation states
List of residues specifiy any non-canonical protonation states. Should be listed as

`["A", 254, "protonate"]`

This will protonate residue 254 in chain A. Alternatively, you can directly specify the type of protonation, or deprotonation by using the index.

`["A", 254, "protonate", 2]`

If residue 254 is an ASP, then we are saying we want to actually protonate OD2, and not OD1. Note that if you do not specify and integer at the end, 
it assumes to act as if you put a 1 there instead. Hence, the initial example is equivalent to doing:

`["A", 254, "protonate", 1]`

##### Frozen atoms
Can specify to freeze a chain, residue, or atom in the various slots. For the chains, just list the chain letter. For the residues, specify `["A", 254]`
which would be chain A and residue 254. For an atom, just include the atom ID, for example `["A", 254, "CA"]` for the alpha carbon in residue 254.

##### Restrict Displacement
This will restrict the displacement between two atoms by a specified amount. Specify atom 1 and atom 2 (as in the Frozen atom section above), and 
finally include the displacement magnitude. For example

`[["A", 254, "SG"], ["A", 230, "SG"], 0.2]`

will restrict the displacement between these two sulfers by 0.2 angstroms. This is a common trick to help model disulfide bridges in proteins!

##### Commands
A way to create multiple runs, with some parameters tweaked, in a single run. For example, say you wanted to have 2 DMD simulations run back to back, 
with the first simulation only going 1000 time steps, and the second one going 200 time steps. And lets say afterwards, you want to gradually increase 
the temperature to 0.2 for 100 time steps. Then you would place the following into the commands section:

`
"1" : {

},
"2" : {
    "Time" : 200
},
"3" : {
    "Time" : 100,
    "Final Temperature" : 0.2
}
`

Notice how I created 3 internal dictionaries. Any keys within these dictionaries have to correspond to a key in the options outside of the commands section
(Time, Initial Temperature, Final Temperature, dt,...). Any new value inside will override what is previously stated. So we first have DMD run a 
simulation with eveything being the default in the dmdinput.json. Then we run for only 200 time steps. Finally, we do a simulation with the final 
temperature at 0.2 and go for 100 time steps. Instead of haveing to create 3 seperate dmdinput.json files and 3 seperate runs, we can string them all
together into one input file!

##### Remaining Commands
Not to be adjusted if you are unsure what you are doing. This is used by dmdpy to keep track of what still needs to be done in the original simulations
that you specified in the commands section.

## Submission Keywords
##### date
The date the script was made.

##### submit_dir
The directory that `submitturbomole.py` was called in.

##### node_type
The type of node to be used. These are listed in `turbopy_config.json` in the `node_types` list.

##### cores
The number of cores that your job will use.

##### time
The amount of time (in hours) that your job will run for.

##### high_priority
Will add `,highp` at this location if a high priority node is used.

##### user
This will be replaced with the user account of the person submitting the job.

##### job_name
This is where the job name will be placed.

##### run_script
This will auto-populate with the correct (absolute) path to the `rundmd.py` script. Following this call should be the 
number of cores, the time (if applicable), and scratch (if applicable).

## Usage
###### setupdmd.py
This script automates the set up of a piDMD job by reformating the pdb and creating all the necessary files. It uses `dmdinput.json` as a 
source of user parameters to generate any files and ensure all parameters are properly set. To run this script, simply call it in a directory 
with a `.pdb` file and a `dmdinput.json`. If there was no `dmdinput.json` in the directory, a default one will be copied over for you to edit.

###### rundmd.py
This script runs DMD simulations job on a user decided number of cores, time, and scratch directory. Call this script in the directory with 
a `.pdb` file and a `dmdinput.json` file. NOTE: it will not copy over a default `dmdinput.json` file if there is not one already present. 
Thus you should ensure that the `dmdinput.json` file has the proper parameters to run the job that you want. It decides which simulations to
perform by looking at `dmdinput.json` at the `commands` section. See [DMD Options](#options) for available options. 

An example usage for 16 cores for 24 hours on a scratch directory `/local/scratch/my_scratch` would be
```bash
rundmd.py -n 16 -t 24 -s /local/scratch/my_scratch
```

It is recommended that you run `setupdmd.py` prior to running `rundmd.py` as issues with piDMD can cause issues with setting up the job. 
Additionally, this script is mainly built to be run on a supercomputer cluster and partners with the `submitdmd.py` script; however it works 
perfectly fine on a local computer given the appropriate resources!

###### submitturbomole.py
This script submits a piDMD job to a queuing system. It uses a [jinja2](https://jinja.palletsprojects.com/en/2.10.x/) template that you provide
as a submission script. The `submit.j2` file is located in the `~/.phd/` direectory. For keywords, see the [Submission Keywords](#Submission-Keywords)
section on how to create the template. Users can then specify the amounts of nodes, cores/nodes, time, and if to submit at execution. 
It is highly recommended to build the job with `setupdmd.py` prior to submitting to the queue as aforementioned issues with piDMD can cause 
early termination. Ensure that a `.pdb` and `dmdinput.json` file are in the directory you are submitting from.

An example usage for 16 cores, for 24 hours with submission is:
```bash
submitdmd.py -n 16 -t 24 -sub
``` 
Note that the default number of nodes to use is 1, but that can be changed using the `-N` flag. 
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Authors
* Matthew Hennefarth - [Alexandrova Research Group](http://www.chem.ucla.edu/~ana/members.html), UCLA

## License
[GNU Public license](https://choosealicense.com/licenses/gpl-3.0/)
