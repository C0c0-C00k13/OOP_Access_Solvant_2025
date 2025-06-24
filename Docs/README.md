This repository contains the executable files* and data to calculate the Solvant Accessible Surface Area of a protein. 

\* NACCESS v2.1.1 has been downloaded from the [official site](https://www.bioinf.manchester.ac.uk/naccess/)

The project has been done in Linux OS and a virtual environment (conda version 24.9.1).

# Repertories
## Data :
Directory containing every files processed or necessary to run NACCESS and SASA.py.
- **environment.yml :** YAML files contaning the dependencies used to execute the files.
- **Insuline.pdb :** Insuline(60sec) and UV laser excited fluorescence.
- **2oe4.pdb :** High Pressure Psuedo Wild Type T4 Lysozyme.
- **2c8r.pdb :** Insuline(60sec) and UV laser excited fluorescence. (recommended for testing the functions) 
- **1bj5.pdb :** HUMAN SERUM ALBUMIN COMPLEXED WITH MYRISTIC ACID.
- **1c26.pdb :** CRYSTAL STRUCTURE OF P53 TETRAMERIZATION DOMAIN.
- **6pwf.pdb :** Cryo-EM structure of the ATPase domain of chromatin remodeling factor ISWI bound to the nucleosome 
- **standard.data :** File containing the standard MAX ASA of each amino acid with a 1.4 anstrom probe atom.
- **vdw.radii :** File containing the radius and polarity of every atoms of each organic molecule/ residue (amino acid, nucleotide, etc..).
 
## Results
Repertory containing the results from NACCESS v2.1.1 & SASA.py.

They are ordered in based on the file used to calculate the ASA (e.i. : *'Results/SASA/'*), the date of the creation of the file (*'Results/SASA/06-24-2025/'*) and the PDB accession ID of the protein (*'Results/SASA/06-24-2025/2c8r'*). 

## Analysis
Python scripts to evaluate and compare the NACCESS and SASA.py. 
- **comparison.py :** Execute the comparisons between naccess and SASA.py
- **atomAnalysis.py :** Contains functions related to the analysis of ASA of atoms.
- **residueAnalysis.py :** Contains functions related to the analysis of ASA of residues.
- **utils :** Contains functions usable for other context than atom or residue (e.i. timed_input function : Prompts the user to give an input within a specific duration).
- **durationAnalysis.py :** Contains the function used to analyze the how much time the the script has run.
- **Figs/ :** Repertory containing the plots shown by the functions from the analysis scripts.

## Root :
Contains the files used to calculate the 
 	- SASA.py : Execution file. Calculates the ASA of the protein (input) and returns 3 output files : <br>
		- **\<output>.asa** : Atom-level ASA.<br>
		- **\<output>.rsa** : Contains Residue-level ASA.<br>
		- **\<output>.log** : Contains details of the execution (e.i : date of exection, parameters etc.)


# Let's get started

## Run NACCESS
In a terminal, enter :

```bash
$ ./naccess ./Data/2c8r.pdb
```
You will get the following message, meaning that the program has successfully ended.
```
naccess: using vdw.radii in local directory
naccess: using STD FILE in local directory
```

## Run SASA.py
### Pre-requiered : if the virtual environment does not exist yet 
Check that the virtual environment `environment` exists:
```bash
$ conda env list
```
If not, create the virtual environment with the command:
```bash
$ conda create -f environment.yml
```

### Steps
1. Activate the environment:
```bash
$ conda activate environment
```

2. Then you can run `SASA.py`. The program will run only if a filename is provided to `pdb_file`. 
```bash
$ python SASA.py
```
If a file is not included in argument, you will get this message: 
```
06-23-2025 15:21:19
usage: SASA.py [-h] [-i {y,n}] [-n POINTS] [-o OUTPUT] [-p PROBE] [-r RADII] pdb_file
SASA.py: error: the following arguments are required: pdb_file
```

* `pdb_file` = PDB file to process.
* `-i`/ `--hetero` = Include heteratoms while calculating the ASA. Default set to n(o).
* `-n`/ `--points` = Number of points used to represent a sphere. Default set 92.
* `-o`/ `--output` = Name of the output files. Default set to 'output'
* `-p`/ `-probe`  = radius of the probe. Default value set to 1.4.
* `-r`/ `--radii`  = Include a customize radii file. Default is set to None.

```bash
$ python SASA.py ./Data/2c8r.pdb
```
Here is the message notifying of a successful run until the end.
```
06-23-2025 15:25:23
PDB file         : Data/2c8r.pdb
Include HETATM   : False
Probe radius     : 1.4
Custom radii file: None
The output files 'output.*' has been created in the directory: Results/SASA/06-23-2025/2c8r
Done
```