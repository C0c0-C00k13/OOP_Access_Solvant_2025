This repository contains the executable files* and data to calculate the Solvant Accessible Surface Area of a protein. 

\* NACCESS v2.1.1 has been downloaded from the [official site](https://www.bioinf.manchester.ac.uk/naccess/)

The project has been done in Linux OS and a virtual environment (conda version 24.9.1).

# Repertories
## Data :
- **environment.yml :** YAML files contaning the dependencies used to execute the files.
- **Insuline.pdb :** Insuline(60sec) and UV laser excited fluorescence.
- **2oe4.pdb :** High Pressure Psuedo Wild Type T4 Lysozyme.
- **2c8r.pdb :** Insuline(60sec) and UV laser excited fluorescence. (recommended for testing the functions) 
- **1bj5.pdb :** HUMAN SERUM ALBUMIN COMPLEXED WITH MYRISTIC ACID.
- **1c26.pdb :** CRYSTAL STRUCTURE OF P53 TETRAMERIZATION DOMAIN.
- **6pwf.pdb :** Cryo-EM structure of the ATPase domain of chromatin remodeling factor ISWI bound to the nucleosome 
 
## Results
Repertory containing the results from NACCESS v2.1.1 & SASA.py.

They are ordered in based on the date of the creation of the file (e.i. : *'Results/06-22-2025'*), the PDB accession of the protein (*'Results/06-22-25/2c8r')* and the file used to calculate the ASA (*'Results/06-22-2025/2c8r/SASA'*). 

## Analysis
Python scripts to evaluate and compare the NACCESS and SASA.py. 
- comparison.py : Execute the comparisons between naccess and SASA.py
- atomAnalysis.py : Contains functions related to the analysis of ASA of atoms.
- residueAnalysis.py : Contains functions related to the analysis of ASA of residues.
- utils : Contains functions usable for other context than atom or residue (e.i. timed_input function : Prompts the user to give an input within a specific duration).

- Execution :
 	- SASA.py : Execution file. Calculates the ASA of the protein (input) and returns 3 output files :
		- **\<output>.asa** : Atom-level ASA.
		- **\<output>.rsa** : Contains Residue-level ASA.
		- **\<output>.log** : Contains details of the execution (e.i : date of exection, parameters etc.)


# Let's get started

## Run NACCESS
In a terminal, enter :

```bash
$ ./naccess ./Data/2c8r.pdb
naccess pdb_file [-p probe_size] [-r vdw_file] [-s stdfile] [-z zslice] -[hwyfaclqb]
```

## Run SASA.py
Check that the virtual environment `environment` exists:
```bash
$ conda env list
```
If not, create the virtual environment with the line:
```bash
$ conda create -f environment.yml
```
Activate the environment:
```bash
$ conda activate environment
$ ./naccess projet-court-POO/Data/2c8r.pdb
naccess: using vdw.radii in local directory
naccess: using STD FILE in local directory
$ mkdir 2c8r
$ mv -i 2c8r.* 2c8r/
$ mv -i 2c8r/ Results/"$(date +"%m-%d-%y")"
```

```bash
$ python SASA.py [-h] [-i {y,n}] [-n POINTS] [-o OUTPUT] [-p PROBE] [-r RADII] pdb_file
```

* pdb_file = PDB file to process.
* -i/ --hetero = Include heteratoms while calculating the ASA. Default set to n(o).
* -n/ --points = Number of points used to represent a sphere. Default set 92.
* -o/ --output = Name of the output files. Default set to 'output'
* -p/ "-probe  = radius of the probe. Default value set to 1.4.
* -r/ --radii  = Include a customize radii file. Default is set to None.