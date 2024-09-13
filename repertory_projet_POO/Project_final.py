"""Script d'exécution du calcul de surface de la protéine exposée au solvant
"""
import numpy as np
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Atom
from Bio.PDB import NACCESS


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

class point_atome:
    """
    """


    def __init__(self, atom_center, x_pt, y_pt, z_pt):
        self.atom_center = atom_center 
        self.x_pt = x_pt
        self.y_pt = y_pt
        self.z_pt = z_pt

    def __str__(self):
        return f"Atome {print(self.atom_center)},\t Point coords[{self.x_pt},{self.y_pt},{self.z_pt}]"

    def calcul_distance(self, atome):
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        if isinstance(atome,Bio.PDB.Atom.Atom):
            distX = abs(self.x_pt - atome.coord[0])
            distY = abs(self.y_pt - atome.coord[1])
            distZ = abs(self.z_pt - atome.coord[2])
            return pow( (pow((distX), 2) + pow((distY),2) + pow((distZ),2)), 0.5)

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
# p1 = point_atome(totale_atoms[0], 1.0,2.0,3.5)
# p1.calcul_distance(totale_atoms[4])
# totale_atoms[0] - totale_atoms[1]

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
    
    protein1 = ("","./Data/insuline.pdb")
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