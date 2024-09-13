Cerépertoire contient l'ensembledes scripts et documents utilisés pour réaliser le calcul de surface de protéine exposée au solvant.

Le projet a été réalisé dans un environnement conda sous un système Linux.

- Data :
	- environment.yml : fichier YAML contenant les dependencies utilisées pour le calcul
	- Insuline.pdb : fichier PDB test utilisé pour calculer la sur
    - 2oe4.pdb : ficheir PDB représentant une structure du lyzozyme
    - 2c8r.pdb : fichier PDB représentant une structure de l'insuline
 
- Results : Contient les résultats de calculs effectués avec le programme NACCESS v2.1.1

- Execution :
 	- Fonctions_utiles.py : Script rassembant toutes les fonctions personnalisées et non liées à une classe;
	- point_atome.py : script définissant le point d'un atome;  
    - Execution_file.py : Fichier d'exécution. C'est lui qui devrait (en théorie) réaliser l'ensemble des calculs.
	- Project_final.py : fichier d'exécution. Rassemble tous les scripts ci-dessus. Réalises également l'ensemble des calculs. 
	
	- Project_2.ipynb : notebook équivalent à Project_final.py. Contient l'ensemble du code pour une observation du script plus détaillée. 
 
