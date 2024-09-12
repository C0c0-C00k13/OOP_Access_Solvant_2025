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


def PDBRetrieve_Atoms(id_prot, filename):
    """
    Récupère les atomes de chacun des résidus
    ---
    Args
    id_prot : str
    filename : str
    ---
    Returns
    list
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


def PDB_Direct_Retrieve_Atoms(id_prot, filename):
    """
    Renvoie les atomes présents dans la fichs PBD. Attention ! Les molécules d'eau sont incluses. 
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
    ensemble_atome = struct.get_atoms()

    # Enregistrement des atomes dans une liste
    list_atome = [atome for atome in ensemble_atome]
    # len(list_atome)
    return list_atome


if __name__ == "__main__": 
    protein1 = ("2C8Q","./Data/insuline.pdb")

    # Test sur residu
    # Observation de la classe residu
    list_res = PDBRetrieve_Residus(protein1[0],protein1[1])

    # Verification : tous les residus de la proteine ont ete recuperes
    print(f"Nombre total de residu : {len(list_res)}")
    print("--------")
    # Observation des attributs pour 1 seul résidu
    r1 = list_res[0]
    print(f"-- 1 seul residu --\nNom residu : {r1.resname},\nId res :{r1.id},\nFull id : {r1.full_id}\n\
    Nom residu avec fonction {r1.get_resname()}")
    print("--------")
    # Observation des attributs pour tous les résidus
    for r1 in list_res:
        print(f"-- Tous les residus --\nNom residu : {r1.resname},\nId res :{r1.id},\nFull id : {r1.full_id}\n\
        Nom residu avec fonction {r1.get_resname()}")

    # Test sur atom en passant par residue
    # Observation des attributs de la classe Calc_Atom à partir de la class 
    list_atome = PDBRetrieve_Atoms(protein1[0],protein1[1])
    # print(list_atome)

    # Consultation d'information sur les atomes
    new_liste_atom = []
    for grp_atom in list_atome:
        for atome in grp_atom:
            # print(f"Coordonnées : {atome.coord}; Full name : {atome.fullname}; Element : {atome.element}")
            new_atom = [atome.fullname, atome.element, atome.coord]
            new_liste_atom.append(new_atom)
    print(new_liste_atom)

    # Test sur direct atom
    # Recuperation de tous les atomes (H20 inclus)
    list_atome2 = PDB_Direct_Retrieve_Atoms(protein1[0],protein1[1])
    # Affichage de id + coord : fonctionnnel
    [[atome.id,atome.coord] for atome in list_atome2]

    # Verification : tous les atomes de la proteine ont ete recuperes
    print(len(list_atome2))
    for a1 in list_atome2:
        print(f"Coordonnées:{a1.coord};\nFull name : {a1.fullname};\nElement : {a1.element}")


# Inspiration source
# https://stackoverflow.com/questions/10324674/parsing-a-pdb-file-in-python?rq=3