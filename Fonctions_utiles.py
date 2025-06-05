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
    
    Parameters
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


def saff_kuijlaars_points(n, center=(0.0,0.0,0.0), radius=0.0):
    """
    Generates a n-points quasi-uniform sphere based on Saff and Kuijlaars algorithm.

    Parameters
    ---
    n (int): Number of points to generate.
    center (np.array) : Coordinates of the center of the sphere. Default value : (0.0,0.0,0.0).
    radius (float) : Radius of the atom. Default value : 0.0.

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


def distance(point_a, point_b)->float:
    """Calculates distance between 2 points
    
    Parameters
    ---
    point_a : Coordinates from either a point of a an atom or an atom.
    point_b : Coordinates from either a point of a an atom or an atom.
    
    Returns
    ---
    distance (float) : Distance from 2 points (a, b).
    """
    return math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3)))


def calculate_asa(current_atom, atoms_list, probe_radius=1.4)->float:
    """Calculates ASA of 1 atom.
    
    Parameters
    ---
    current_atom (Atom) : Atom currently processed. Its attribute 'points' \
    cannot be empty (type : (np.array)).
    atoms_list ([Atom]) : List of atoms from the PDB file.
    probe_radius (float) : Radius of the probe (Oxygen atom). Default value : 1.4 Å.

    Returns
    ---
    asa (float) : ASA of the current atom.
    """

    # Radius of current atom
    threshold_radius = current_atom.radius + probe_radius
    accessible_points = 0

    # Runs through the list of points of current atom
    for point in current_atom.points:
        exposed = True
        # Runs through the list of atoms
        for atom_i in atoms_list:
            if current_atom.index != atom_i.index:
                if distance(point, atom_i.position) < threshold_radius:
                    exposed = False
                    break
                if exposed:
                    accessible_points += 1
    # Surface area of full sphere * exposed fraction
    sphere_area = 4 * math.pi * threshold_radius**2
    return sphere_area * (accessible_points / len(current_atom.points))


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
# rsa_data = []
# for res_id, aa, asa in asa_data:
#     max_val = max_asa.get(aa)
#     if max_val:
#         rsa = asa / max_val
#         rsa_data.append((res_id, aa, asa, round(rsa, 3)))
#     else:
#         rsa_data.append((res_id, aa, asa, None))

# --- EXECUTION
# Print results
# for res_id, aa, asa, rsa in rsa_data:
#     print(f"Residue {res_id} ({aa}): ASA = {asa:.2f}, RSA = {rsa if rsa is not None else 'N/A'}")

# ------------------------------------
# ---------------------------------------------------------------------------------

if __name__ == "__main__":


    # Display the current date of run
    today = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    # -------------------------------------------
    # CHECK PDB FILE
    FILE_PDB = "Data/2c8r.pdb"
    # Checking if the file exists
    IS_EXIST = os.path.exists(FILE_PDB)
    print(f"Date of execution : {today}")
    if not IS_EXIST:
        print("This file does not exist. Please check the file path.\
            \nExit")
        sys.exit()

    # -------------------------------------------
    # CREATE RADII REFERENCE
    FILE_RADIUS = "./Data/vdw.radii"

    IS_EXIST = os.path.exists(FILE_RADIUS)
    if IS_EXIST:
        # print('MAIN EXECUTION')
        # print(f"Radii references: '{FILE_RADIUS}' found...")
        time.sleep(1)
        VAN_DER_WAALS_RADII = get_radii(FILE_RADIUS)
    else:
        print(f"This file does not exist. Default values:\n{VAN_DER_WAALS_RADII}.")

    # time.sleep(2)
    # -------------------------------------------
    # READ PDB FILE
    print("Reading PDB file...")
    atoms = get_atoms(FILE_PDB)
    # print(atoms)
    print("Reading PDB file - Done.")
    atom_1, atom = atoms[0], atoms[1]
    # print(atom)

    # time.sleep(2)
    # -------------------------------------------
    # CALCULATE ASA
    # -- GENERAL TEST
    NUMBER_OF_POINTS = 92
    print("Calculating ASA...")
    for atom_i in atoms[:10]:
        print(atom_i)
        # GENEREATE SPHERE
        # print("Generating Sphere points...")
        atom_i.points = saff_kuijlaars_points(NUMBER_OF_POINTS, atom_i.position, atom_i.radius)
        # print("Generating Sphere points - DONE")
        atom_i.asa = calculate_asa(current_atom=atom_i, atoms_list=atoms)
        # print(atom_i)
    print("Calculating ASA - Done")

    # -------------------------------------------
    # CALCULATE RSA
    rsa_data = []
    for res_id, aa, asa in asa_data:
        max_val = max_asa.get(aa)
        if max_val:
            rsa = asa / max_val
            rsa_data.append((res_id, aa, asa, round(rsa, 3)))
        else:
            rsa_data.append((res_id, aa, asa, None))

    # Print results
    for res_id, aa, asa, rsa in rsa_data:
        print(f"Residue {res_id} ({aa}): ASA = {asa:.2f},\
RSA = {rsa if rsa is not None else 'N/A'}")


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
