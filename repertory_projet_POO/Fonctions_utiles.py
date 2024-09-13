"""Script d'exécution du calcul de surface de la protéine exposée au solvant
"""

import numpy as np
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Atom
from Bio.PDB import NACCESS

import point_atome


def PDBRetrieve_Residus(id_prot, filename):
    """
    Fonction permettant de récupere les résidus d'une protéine
    ---
    Args
    id_prot : str
    filename : str
    ---
    Returns
    list
    """

    
    pdbparser = Bio.PDB.PDBParser(QUIET=True)
    struct = pdbparser.get_structure(id_prot, filename)

    # Recuperation des residus
    ensemble_res = struct.get_residues()
    
    list_res = [res for res in ensemble_res]
    # Recuperation des residus, exlusion des molecules d'eau
    return [res for res in list_res if not res.resname == "HOH" ]

# ------------------------
def PDBRetrieve_Atoms(id_prot, filename):
    """
    """

    # Recuperation des residus
    # print("Get residues")
    list_res = PDBRetrieve_Residus(id_prot, filename)
    # Recuperation des atomes + information par atome/ residu
    # print("Get atom generator")
    ensemble_atome = [res.get_atoms() for res in list_res]
    # print(f"Liste de generateur : {ensemble_atome}") # Visible
    # Enregistrement des atomes dans une liste
    list_atome = [[at for at in atome] for atome in ensemble_atome]
    # print(f"Nombre total d'atome : {len(list_atome)};\nListe d'atome : {list_atome}")
    return list_atome

# ------------------------
def saff_kuijlaars_points(N):
    """
    Génère N points quasi-uniformes sur une sphère unitaire
    à l'aide de l'algorithme de Saff et Kuijlaars.

    Args:
        N (int): Nombre de points à générer.

    Returns:
        points (ndarray): Un tableau (N, 3) avec les coordonnées x, y, z des points.
    """
    
    points = np.zeros((N, 3))
    
    for k in range(1, N + 1):
        h = -1 + 2 * (k - 1) / (N - 1)  # Hauteur du point
        theta = np.arccos(h)            # Colatitude
        phi = np.pi * (1 + np.sqrt(5)) * (k - 1)  # Longitude (angle d'or)
        
        # Coordonnées sphériques vers cartésiennes
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        
        points[k - 1] = np.array([x, y, z])
    
    return points

# ------------------------
def minifuction(atome,points):
    """Génère une liste de points représantant un atome de la proteine
    Parameters
    ------
    atome : Bio.PDB.Atom
        Atome auquel les points seront attachés.
    points : list
        Liste de points. 
    Returns
    ------
    list
        Liste des contenant les points rattachés à l'atome.
    """
    liste_point_coord = []
    for index,coord in enumerate(points):
        new_point = point_atome(atom_center=atome,x_pt=coord[0]+atome.coord[0], y_pt=coord[1]+atome.coord[1], z_pt= coord[2]+atome.coord[2])
        liste_point_coord.append(new_point)
        # print(index, new_point)
    return (liste_point_coord)
    # print(len(liste_point_coord), liste_point_coord)

# ------------------------
def comparaisonDistances(atome, pts_atome, totale_atoms):
    """Fonction renvoyant la liste des points exposés au solvant en fonction de leurs distances au reste des atomes
    Parameters
    ------
    atome : Bio.PDB.Atom
        Atome auquel les points sont reliés
    pts_atom : list
        Liste des points reliés à l'objet atome
    totale_atoms : list
        Ensemble des atomes contenus dans la protéine
    Returns
    ------
    list
        Liste des points exposés au solvant rattachés à atome. 
    """
    
    # Le seuil fixe correspond à la taille d'un atome d'oxygène
    liste_points_solvant = [] 
    SEUIL = 1.4
    # Comparer les distances
    for pt in pts_atome:
        # Au debut de la comparaison, un veariable vaut 0.
        condition_solvate = True
        for at_tot in totale_atoms:
            # Calcul des points pour chaque atome + renvoei d'une liste de points
            if at_tot == atome:
            # print(f"atome trouve en position {cpt_tot_at}")
                pass
            else:
                # print(f"Residu/atome {cpt_at_per_res+1, atome.element}; Distance {pt.calcul_distance(atome=at_tot)}")
                if pt.calcul_distance(atome=at_tot) < SEUIL:
                    condition_solvate = False
        if condition_solvate:
            liste_points_solvant.append(pt)
    return liste_points_solvant

# ------------------------
# fonction renvoyant les points exposés au solvant pour un residu
def Exposition_point_par_solvant(list_atome):
    """Associe à chaque atome l'ensemble de ses points exposés au solvant dans un nouvel attribut: liste_points_solvant.
    Parameters
    ------
    list_atome : list
        liste des atomes de la protéine. Les atomes sont groupés par les résidu.
    Returns
    ------
    None
    """
    # Liste de tous les atomes de protéines. Non séparés par résidu
    TOTALE_ATOMS = []
    for res in list_atome:
        TOTALE_ATOMS = TOTALE_ATOMS + res
    POINTS = saff_kuijlaars_points(92)
    
    for res in list_atome:
        for atome in res:
            atome.liste_points_solvant = []
            pts_atome = minifuction(atome=atome, points=POINTS)
            atome.liste_points_solvant = comparaisonDistances(atome, pts_atome, TOTALE_ATOMS)
    # for res in list_atome:
    #     print(res)
    #     for atome in res:
    #         print(atome,len(atome.liste_points_solvant))

if __name__ == "__main__":
    pass
