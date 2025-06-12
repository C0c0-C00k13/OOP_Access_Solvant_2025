"""Script d'exécution du calcul de surface de la protéine exposée au solvant
"""
import os
import sys
import time
import datetime
import math
import numpy as np

from Atom import Atom
from Get__Radii import get_radii, get_reference_total_asa


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
                # print(line.strip().split())
                # Index
                index = int(line.strip().split()[1])
                # Element
                element = line.strip().split()[2]
                atom_type = line.strip().split()[-1]
                # Residue
                residue = line.strip().split()[3]
                chain = line.strip().split()[4]
                id_res = int(line.strip().split()[5])
                # Coordinates
                coord_z = float(line.strip().split()[8])
                coord_y = float(line.strip().split()[7])
                coord_x = float(line.strip().split()[6])
                position = (coord_x,coord_y,coord_z)
                # print(f"Index:{index}; Position:{position}; Element:{element}")
                atom = Atom(element=element, atom_type=atom_type,chain=chain,
                            id_res=id_res,residue=residue,position=position,
                            index=index)
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
        # x = round((center[0] + math.sin(theta) * math.cos(phi) * radius), 3)
        # y = round((center[1] + math.sin(theta) * math.sin(phi) * radius), 3)
        # z = round((center[2] + math.cos(theta) * radius), 3)

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
    # return math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3)))
    return round(math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3))), 5)


def comparison_distance(central_atom, test_atom,central_point, probe_radius = 1.4):
    threshold_radius = central_atom.radius + probe_radius
    occluded =  False
    for pos_points in test_atom.points:
        distance_pt_at = distance(pos_points, central_atom.position)
        distance_pt_pt = distance(pos_points, central_point)
        if distance_pt_at < threshold_radius and distance_pt_pt < probe_radius:
            occluded = True
            break
    return occluded


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
                    # break
                if exposed:
                    accessible_points += 1
    print(current_atom,accessible_points)
    # Surface area of full sphere * exposed fraction
    sphere_area = 4 * math.pi * threshold_radius**2
    return sphere_area * (accessible_points / len(current_atom.points))


def read_max_asa(filename:str):
    """Returns max ASA for each residue."""
    # Max ASA values (Tien et al. 2013)
    max_asa = {
        'ALA': 121.0, 'ARG': 265.0, 'ASN': 187.0, 'ASP': 187.0,
        'CYS': 148.0, 'GLN': 214.0, 'GLU': 214.0, 'GLY': 97.0,
        'HIS': 216.0, 'ILE': 195.0, 'LEU': 191.0, 'LYS': 230.0,
        'MET': 203.0, 'PHE': 228.0, 'PRO': 154.0, 'SER': 143.0,
        'THR': 163.0, 'TRP': 264.0, 'TYR': 255.0, 'VAL': 165.0
    }
    IS_EXIST = os.path.exists(filename)
    if IS_EXIST:
        print(f"Max ASA references: '{filename}' found...")
        time.sleep(1)
        max_asa = {}
        with open(filename, 'r') as standard_file:

            for line in standard_file:
                if line.startswith("ATOM"):
                    tab_line = line.strip().split()
                    residue = tab_line[3]
                    max_asa.update({residue : tab_line[4:]})
    else:
        print(f"This file does not exist. Returns default values.")
    return max_asa


def calculate_residue_asa(list_atoms):
    """Calculates ASA of each residue"""
    residues_asa = {}
    # Run through atom list
    for atom_i in list_atoms:
        # Searches for the current residue
        if atom_i.residue in residues_asa:
            # Searches for the current residue index
            if atom_i.id_res in residues_asa[atom_i.residue]:
                new_asa_res = residues_asa[atom_i.residue][atom_i.id_res] + atom_i.asa
                residues_asa[atom_i.residue][atom_i.id_res] = new_asa_res
            # Creates a new emplacement for the residue index
            else:
                residues_asa[atom_i.residue].update({atom_i.id_res : atom_i.asa})
        # Creates a new emplacement for the residue name
        else:
            residues_asa[atom_i.residue] = {atom_i.id_res : atom_i.asa}

    # Rounding the ASA to 3 decimals
    for residue in residues_asa:
        for index in residues_asa[residue]:
            residues_asa[residue][index] = round(residues_asa[residue][index], 3)
    return residues_asa


def calculate_max_asa(list_atoms, probe_radius:float):
    """Calculates the max ASA of each residue."""
    max_asa = {}
    # Run through atom list
    for atom_i in list_atoms:
        max_asa_atom = 4 * math.pi * (atom_i.radius + probe_radius)**2
        # Searches for the current residue
        if atom_i.residue in max_asa:
            # Searches for the current residue index
            if atom_i.id_res in max_asa[atom_i.residue]:
                new_max_asa_res = max_asa[atom_i.residue][atom_i.id_res] + max_asa_atom
                max_asa[atom_i.residue][atom_i.id_res] = new_max_asa_res
            # Creates a new emplacement for the residue index
            else:
                max_asa[atom_i.residue].update({atom_i.id_res : max_asa_atom})
        # Creates a new emplacement for the residue name
        else:
            max_asa[atom_i.residue] = {atom_i.id_res : max_asa_atom}
    return max_asa


# Example ASA data (residue index, residue name, ASA value)
asa_data = [
    (1, 'A', 55.0),
    (2, 'R', 120.0),
    (3, 'G', 45.0),
    (4, 'V', 150.0),
    (5, 'L', 100.0)
]

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
    # CREATE RADII REFERENCES
    FILE_RADIUS = "./Data/vdw.radii"
    VAN_DER_WAALS_RADII = get_radii(FILE_RADIUS)

    # TOTAL ASA
    # --- READ
    # Reading Total ASA
    print("Reading Total ASA...")
    # filename_data = "./Data/standard.data"
    # max_asa = read_max_asa(filename=filename_data)
    print("Reading Total ASA file - Done.")

    # --- CALCULATING
    # Calculating Total ASA
    print("Calculating Total ASA...")
    PROBE_RADIUS = 1.4
    max_asa = get_reference_total_asa(filename=FILE_RADIUS, probe_radius=PROBE_RADIUS)
    print("Calculating Total ASA - Done.")
    # print(max_asa)

    # time.sleep(2)
    # -------------------------------------------
    # READ PDB FILE
    print("Reading PDB file...")
    atoms = get_atoms(FILE_PDB)
    # print(atoms)
    print("Reading PDB file - Done.")
    # atom_1, atom = atoms[0], atoms[1]
    # print(atom)
    # time.sleep(2)

    # -------------------------------------------
    # CALCULATE ASA
    # -- GENERAL TEST
    NUMBER_OF_POINTS = 92
    print("Calculating Atomic ASA...")
    for atom_i in atoms:
        # GENEREATE SPHERE
        # print("Generating Sphere points...")
        atom_i.points = saff_kuijlaars_points(NUMBER_OF_POINTS, atom_i.position, atom_i.radius)
        # print("Generating Sphere points - DONE")
        # time.sleep(2)

        # atom_i.asa = calculate_asa(current_atom=atom_i, atoms_list=atoms[:10])
        # print(atom_i)
        # Radius of current atom

    # ----------- TO REWRITE
    # List of atoms
    for central_atom in atoms[:8]:
        print(f"atom : {central_atom.element} , n°{central_atom.index}")
        # time.sleep(1)
        # Number of accessible points per atom
        accessible_points = 0
        # print(central_atom.points)

        # Runs through the list of points of Central atom
        for current_point in central_atom.points:
            # ---- TAG
            # print(f"current_point: {current_point}")
            # time.sleep(2)
            # --
            # Point defined as exposed by default
            exposed = True

            # Filter out the atom
            filtered_atoms = [atom for atom in atoms if atom != central_atom]
            # Runs through the list of atoms -- excluding the Central one
            for test_atom in filtered_atoms:
                # ---- TAG
                # print(f"Atom: {test_atom.index}")
                # time.sleep(1)
                # --

                # Compare distance of every points of Test atom to threshold
                if comparison_distance(central_atom, test_atom, current_point):
                    exposed = False
                    # print(distance(pos_points, test_atom.position))
                    # print(f"Current atom {central_atom.index}, point :{current_point}, occluded by :{test_atom.index}, point: {idx_point}")
                    # time.sleep(2)
                    # --
                    break
            if exposed:
                accessible_points += 1
            
        print(central_atom.element,central_atom.index,central_atom.residue,accessible_points)
        # Surface area of full sphere * exposed fraction
        sphere_area = 4 * math.pi * (central_atom.radius + PROBE_RADIUS)**2
        print (sphere_area * (accessible_points / len(central_atom.points)))
    print("Calculating Atomic ASA - Done.")
    # time.sleep(2)

    # # Calculating Residue ASA
    # print("Calculating Residue ASA...")
    # # residues_asa = calculate_residue_asa(atoms[:10])
    # residues_asa = {}
    # # Run through atom list
    # for atom_i in atoms[:5]:
    #     print(atom_i, )
        # # Searches for the current residue
        # if atom_i.residue in residues_asa:
        #     # Searches for the current residue index
        #     if atom_i.id_res not in residues_asa[atom_i.residue]:
        #         # Creates a new emplacement for the residue index
        #         residues_asa[atom_i.residue].update({atom_i.id_res : 0})
        # # Creates a new emplacement for the residue name
        # else:
        #     print(f"New residue: {atom_i.residue}|Current ASA: {residues_asa}")
        #     residues_asa[atom_i.residue] = {atom_i.id_res : 0}
        # residues_asa[atom_i.residue][atom_i.id_res] += atom_i.asa
        # print(f"Current loop: {atom_i.index}:{atom_i.element}:{atom_i.residue}|Current ASA: {residues_asa}")

    # Rounding the ASA to 3 decimals
    # for residue in residues_asa:
    #     for index in residues_asa[residue]:
    #         residues_asa[residue][index] = round(residues_asa[residue][index], 3)
    # print("Calculating Residue ASA - Done.")
    # print(residues_asa)
    time.sleep(2)
    # -------------------------------------------
#     # CALCULATE RSA
#     rsa_data = []
#     for res_id, aa, asa in asa_data:
#         max_val = max_asa.get(aa)
#         if max_val:
#             rsa = asa / max_val
#             rsa_data.append((res_id, aa, asa, round(rsa, 3)))
#         else:
#             rsa_data.append((res_id, aa, asa, None))

#     # Print results
#     for res_id, aa, asa, rsa in rsa_data:
#         print(f"Residue {res_id} ({aa}): ASA = {asa:.2f},\
# RSA = {rsa if rsa is not None else 'N/A'}")


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

    print("Done.")
