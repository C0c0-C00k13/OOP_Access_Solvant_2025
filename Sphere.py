"""Class representing sphere objects"""
import os
import sys
import time
import datetime
import numpy as np
from point_atome import PointAtom
from Atom import Atom


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
# p1 = PointAtom(totale_atoms[0], 1.0,2.0,3.5)
# p1.calcul_distance(totale_atoms[4])
# totale_atoms[0] - totale_atoms[1]

if __name__ == "__main__":

    # Display the current date of run
    today = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    FILE = "Data/2c8r.pdb"
    # Checking if the file exists
    IS_EXIST = os.path.exists(FILE)
    print(f"Date of execution : {today}")
    if not IS_EXIST:
        print("This file does not exist. Please check the file path.\
            \nExit")
        sys.exit()

    # In case the file is found
    print(f"File found: {FILE}. Opening now")
    time.sleep(2)
    with open(FILE, "r") as pdb_file:
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
        position = (coord_x,coord_y,coord_y)
        # Element
        element = line.strip().split()[-1]
        # print(f"Index:{index}; Position:{position}; Element:{element}")
        atom = Atom(element=element,position=position,index=index)
        # print(atom.__dict__)

        print("Generating sphere...")
    print("Done")
