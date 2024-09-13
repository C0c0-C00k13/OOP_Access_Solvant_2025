"""Script d'exécution du calcul de surface de la protéine exposée au solvant.
"""

# Importation de modules
import numpy as np
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Atom
from Bio.PDB import NACCESS

import Fonctions_utiles
import point_atome


if __name__ == "__main__":
    
    protein1 = ("2C8Q","./Data/insuline.pdb")
    print("Début de lecture du fichier PDB.")
    list_atome = PDBRetrieve_Atoms(protein1[0],protein1[1])
    
    print("Fin de lecture du fichier PDB.")
    print("Début du calcul d'exposition dela protéine au solvant")
    Exposition_point_par_solvant(list_atome=list_atome)
    
    # Liste de tous les atomes de protéines. Non séparés par résidu
    TOTALE_ATOMS = []
    for res in list_atome:
        TOTALE_ATOMS = TOTALE_ATOMS + res
    TOTAL_POINTS = 92 * len(TOTALE_ATOMS)
    
    # Pourcentage de la protéine esposée au solvant
    solvated_region = 0
    for atome in TOTALE_ATOMS:
        solvated_region += len(atome.liste_points_solvant)
    
    # Pourcentage de la protéine esposée au solvant par résidu
    solvated_region_2 = []
    for res in list_atome:
        solvated_region_per_res = 0
        total_point_per_res = 92 * len(res)
        for atome in res:
            solvated_region_per_res += len(atome.liste_points_solvant)
        # Pourcentage duu résidu exposé au solvant
        tmp = solvated_region_per_res/ total_point_per_res * 100
        solvated_region_2.append(tmp)
            
             
    print(f"Proportion de protéine au solvant exposée :{solvated_region/(TOTAL_POINTS)*100}%.")
    print(f"Proportion exposées par résidu:")
    for idx, region in enumerate(solvated_region_2):
        print(f"Proportion exposée du résidu {idx} : {region}%")
    print("Fin d'éxecution.")