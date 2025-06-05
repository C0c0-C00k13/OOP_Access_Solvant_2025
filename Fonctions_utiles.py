"""Script d'exécution du calcul de surface de la protéine exposée au solvant
"""
import os
import sys
import time
import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Atom import Atom
from Get__Radii import get_radii


# Defaults Radii
VAN_DER_WAALS_RADII = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8
}

def get_atoms(file:str):
    """Returns list of atoms from a PDB file.
    
    Args
    ---
    file : str
    
    Return
    ---
    atoms : set
    """

    atoms = []
    with open(file, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                # Index
                index = int(line.strip().split()[1])
                # Coordinates
                coord_z = float(line.strip().split()[-6])
                coord_y = float(line.strip().split()[-5])
                coord_x = float(line.strip().split()[-4])
                position = (coord_x,coord_y,coord_z)
                # Element
                element = line.strip().split()[-1]
                atom = Atom(element=element,position=position,index=index)
                atoms.append(atom)

    return tuple(atoms)

def get_1st_atom(file:str):
    """Returns the 1st atom of the PDB file."""
    with open(file, "r") as pdb_file:
        line = pdb_file.readline()
        while not line.startswith("ATOM"):
            line = pdb_file.readline()
        # print(line.strip())
        # Index
        index = int(line.strip().split()[1])
        # Coordinates
        coord_z = float(line.strip().split()[-6])
        coord_y = float(line.strip().split()[-5])
        coord_x = float(line.strip().split()[-4])
        position = (coord_x,coord_y,coord_z)
        # Element
        element = line.strip().split()[-1]
        # print(f"Index:{index}; Position:{position}; Element:{element}")
        atom = Atom(element=element,position=position,index=index)
        # print(atom.__dict__)
        return atom

# ------------------------
def saff_kuijlaars_points(n, center=(0.0,0.0,0.0), radius=0.0):
    """
    Génère N points quasi-uniformes sur une sphère unitaire
    à l'aide de l'algorithme de Saff et Kuijlaars.

    Args
    ---
    n (int): Nombre de points à générer.

    Returns
    ---
    points (ndarray): Un tableau (n, 3) avec les coordonnées x, y, z des points.
    """

    points = np.zeros((n, 3))

    for k in range(1, n + 1):
        h = -1 + 2 * (k - 1) / (n - 1)  # Hauteur du point
        theta = math.acos(h)            # Colatitude
        if k == 1 or k == n:
            phi = 0.0
        else:
            phi += 3.6 / math.sqrt(n * (1 - h * h))

        # Coordonnées sphériques vers cartésiennes
        x = center[0] + math.sin(theta) * math.cos(phi) * radius
        y = center[1] + math.sin(theta) * math.sin(phi) * radius
        z = center[2] + math.cos(theta) * radius

        points[k - 1] = np.array([x, y, z])

    return points

def display_sphere(atom:Atom):
    """Displays the sphere representing an atom"""
    points = atom.points
    # Plotting
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=3, alpha=0.6)
    ax.set_title(f"Sphere of {len(atom.points)} Points\
                (center={atom.position}, radius={atom.radius})")
    ax.set_box_aspect([1, 1, 1])
    plt.show()

def distance(point_a, point_b)->float:
    """Calculates distance between 2 points"""
    return math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3)))

def comparaisonDistances(atom:Atom, list_points, totale_atoms):
    """Returns a list of points exposed to solvant based on
    distances au reste des atomes.

    Parameters
    ------
    atom : Atom
        Atome auquel les points sont reliés
    pts_atom : list
        Liste des points reliés à l'objet atom
    totale_atoms : list
        Ensemble des atomes contenus dans la protéine
    Returns
    ------
    list
        Liste des points exposés au solvant rattachés à atom. 
    """

    # Le seuil fixe correspond à la taille d'un atome d'oxygène
    liste_points_solvant = [] 
    SEUIL = 1.4
    # Comparer les distances
    for pt in list_points:
        # Au debut de la comparaison, un veariable vaut 0.
        condition_solvate = True
        for at_tot in totale_atoms:
            # Calcul des points pour chaque atome + renvoei d'une liste de points
            if at_tot == atom:
            # print(f"atom trouve en position {cpt_tot_at}")
                pass
            else:
                # print(f"Residu/atom {cpt_at_per_res+1, atom.element};\
                # Distance {pt.calcul_distance(atom=at_tot)}")
                if pt.calcul_distance(atom=at_tot) < SEUIL:
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
    points = saff_kuijlaars_points(92)

    for res in list_atome:
        for atome in res:
            atome.liste_points_solvant = []
            # list_points = minifuction(atome=atome, points=points)
            # atome.liste_points_solvant = comparaisonDistances(atome, list_points, TOTALE_ATOMS)
    # for res in list_atome:
    #     print(res)
    #     for atome in res:
    #         print(atome,len(atome.liste_points_solvant))


# ---------------------------------------------------------------------------------
# ------------------------------------
# import math

# Example input: list of atoms (x, y, z, element)
atoms = [
    (0.0, 0.0, 0.0, 'C'),
    (2.0, 0.0, 0.0, 'O'),
    (0.0, 2.0, 0.0, 'N'),
]

# van der Waals radii (Å)
vdw_radii = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80
}

probe_radius = 1.4  # Water probe radius
n_surface_points = 100  # More points = better accuracy

def generate_sphere_points(n):
    points = []
    for i in range(n):
        phi = math.acos(1 - 2*(i+0.5)/n)
        theta = math.pi * (1 + 5**0.5) * (i+0.5)
        x = math.sin(phi) * math.cos(theta)
        y = math.sin(phi) * math.sin(theta)
        z = math.cos(phi)
        points.append((x, y, z))
    return points

sphere_points = generate_sphere_points(n_surface_points)



def calculate_asa(atoms):
    asa_per_atom = []
    for i, (x, y, z, element) in enumerate(atoms):
        r = vdw_radii.get(element, 1.5) + probe_radius
        accessible_points = 0
        for dx, dy, dz in sphere_points:
            px, py, pz = x + r*dx, y + r*dy, z + r*dz
            exposed = True
            for j, (x2, y2, z2, element2) in enumerate(atoms):
                if i == j:
                    continue
                r2 = vdw_radii.get(element2, 1.5) + probe_radius
                if distance((px, py, pz), (x2, y2, z2)) < r2:
                    exposed = False
                    break
            if exposed:
                accessible_points += 1
        # Surface area of full sphere * exposed fraction
        sphere_area = 4 * math.pi * r**2
        asa = sphere_area * (accessible_points / n_surface_points)
        asa_per_atom.append((i, element, round(asa, 2)))
    return asa_per_atom

# --- EXECUTION
# Run ASA calculation
# asa_results = calculate_asa(atoms)

# Output results
# for i, elem, asa in asa_results:
#     print(f"Atom {i} ({elem}): ASA = {asa} Å²")


# Max ASA values (Tien et al. 2013)
max_asa = {
    'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0,
    'C': 148.0, 'Q': 214.0, 'E': 214.0, 'G': 97.0,
    'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
    'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0,
    'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0
}

# Example ASA data (residue index, residue name, ASA value)
asa_data = [
    (1, 'A', 55.0),
    (2, 'R', 120.0),
    (3, 'G', 45.0),
    (4, 'V', 150.0),
    (5, 'L', 100.0)
]

# Calculate RSA
rsa_data = []
for res_id, aa, asa in asa_data:
    max_val = max_asa.get(aa)
    if max_val:
        rsa = asa / max_val
        rsa_data.append((res_id, aa, asa, round(rsa, 3)))
    else:
        rsa_data.append((res_id, aa, asa, None))

# --- EXECUTION
# Print results
# for res_id, aa, asa, rsa in rsa_data:
#     print(f"Residue {res_id} ({aa}): ASA = {asa:.2f}, RSA = {rsa if rsa is not None else 'N/A'}")

# ------------------------------------
# ---------------------------------------------------------------------------------

if __name__ == "__main__":


    # Display the current date of run
    today = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    FILE_PDB = "Data/2c8r.pdb"
    # Checking if the file exists
    IS_EXIST = os.path.exists(FILE_PDB)
    print(f"Date of execution : {today}")
    if not IS_EXIST:
        print("This file does not exist. Please check the file path.\
            \nExit")
        sys.exit()

    # List radii
    FILE_RADIUS = "./Data/vdw.radii"

    IS_EXIST = os.path.exists(FILE_RADIUS)
    if IS_EXIST:
        # print('MAIN EXECUTION')
        # print(f"Radii references: '{FILE_RADIUS}' found...")
        time.sleep(1)
        VAN_DER_WAALS_RADII = get_radii(FILE_RADIUS)
    else:
        print(f"This file does not exist. Default values:\n{VAN_DER_WAALS_RADII}.")
    time.sleep(2)

    atom_1 = get_1st_atom(FILE_PDB)
    # print(atom_1)

    # print("Début de lecture du fichier PDB.")
    atoms = get_atoms(FILE_PDB)
    # print("Fin de lecture du fichier PDB.")
    atom = atoms[1]
    # print(atom)


    print("Generating Sphere points...")
    NUMBER_OF_POINTS = 92
    # time.sleep(2)
    # print(type(atom), atom.radius, type(atom.radius))
    atom.points = saff_kuijlaars_points(NUMBER_OF_POINTS,atom.position,atom.radius)
    # print(atom.points)

    for atom_i in atoms:
        if atom_i == atom:
            print("Same atom selected. Skip...")
            time.sleep(2)
        else:
            distance_i = distance(atom.position,atom_i.position)
            print(f"Distance between atom n°{atom.index} and n°{atom_i.index}: {distance_i}")


    # print("Début du calcul d'exposition dela protéine au solvant")
    # Exposition_point_par_solvant(list_atome=list_atome)

    # # Liste de tous les atomes de protéines. Non séparés par résidu
    # TOTALE_ATOMS = []
    # for res in list_atome:
    #     TOTALE_ATOMS = TOTALE_ATOMS + res
    # TOTAL_POINTS = 92 * len(TOTALE_ATOMS)

    # # Pourcentage de la protéine esposée au solvant
    # solvated_region = 0
    # for atome in TOTALE_ATOMS:
    #     solvated_region += len(atome.liste_points_solvant)

    # # Pourcentage de la protéine esposée au solvant par résidu
    # solvated_region_2 = []
    # for res in list_atome:
    #     solvated_region_per_res = 0
    #     total_point_per_res = 92 * len(res)
    #     for atome in res:
    #         solvated_region_per_res += len(atome.liste_points_solvant)
    #     # Pourcentage duu résidu exposé au solvant
    #     tmp = solvated_region_per_res/ total_point_per_res * 100
    #     solvated_region_2.append(tmp)

    # print(f"Proportion de protéine au solvant exposée :{solvated_region/(TOTAL_POINTS)*100}%.")
    # print(f"Proportion exposées par résidu:")
    # for idx, region in enumerate(solvated_region_2):
    #     print(f"Proportion exposée du résidu {idx} : {region}%")

    print("Fin d'éxecution.")
