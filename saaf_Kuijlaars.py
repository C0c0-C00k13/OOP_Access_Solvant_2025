import numpy as np
import Point

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


if __name__ == "__main__":
    protein1 = ("2C8Q","./Data/insuline.pdb")

    #Recupération des atomes
    list_atome = PDBRetrieve_Atoms(protein1[0],protein1[1])

    # Génération de la classe atome définie par Atom.py
    new_liste_atom = []
    for grp_atom in list_atome:
        for atome in grp_atom:
            # print(f"Coordonnées : {atome.coord}; Full name : {atome.fullname}; Element : {atome.element}")
            new_atom = Calc_Atom(atome.fullname, atome.element, atome.coord)
            new_liste_atom.append(new_atom)
    new_liste_atom[0]
    atome1, atome2 = new_liste_atom[0], new_liste_atom[1]

    # Génération de points à l'aide de l'algorithme de Saff & Kuijlaars
    points = saff_kuijlaars_points(92)
    
    # Test sur 1 atome
    print(f"Index\tCoordonnées initial\tCordonnées centrées")
    for index,coord in enumerate(points):
        print(f"{index}\t{coord[0:3]}\t{coord[0:3]+[atome1.x,atome1.y,atome1.z]}")

    # Test sur tout les atomes
    liste_point_coord = []
    for atome in new_liste_atom:
        for coord in points:
            new_point = point_atome(atom_center= atome, x_pt= coord[0]+atome.x, y_pt= coord[1]+atome.y, z_pt= coord[2]+atome.z)
            liste_point_coord.append(new_point)
            
    print(len(liste_point_coord),liste_point_coord)

