"""
atom.py

Defines the Atom class used for representing atoms in a protein structure for 
solvent-accessible surface area (ASA) calculations. 

Each Atom includes:
- Element type (e.g., C, N, O)
- 3D position (numpy array)
- van der Waals radius (based on element)
- ASA (initialized to 0.0, calculated later)
- Optional atom index for tracking

The Atom class provides methods for computing distances and easy representation.
"""


import os
import sys
import time
import datetime
import numpy as np
import Get__Radii


# Dictionary of van der Waals radii (in Ångströms)
FILENAME = "./Data/vdw.radii"
# Checking if the file exists
IS_EXIST = os.path.exists(FILENAME)
if IS_EXIST:
    print(f"Radii references: '{FILENAME}' found...")
    time.sleep(2)
    VAN_DER_WAALS_RADII = Get__Radii.get_radii(FILENAME)
else:
    VAN_DER_WAALS_RADII = {
        'H': 1.2,
        'C': 1.7,
        'N': 1.55,
        'O': 1.52,
        'S': 1.8
    }
    print(f"This file does not exist. Default values:\n{VAN_DER_WAALS_RADII}.")
    time.sleep(2)


class Atom:
    """
    A class to represent a single atom in a molecular structure.

    Attributes:
        element (str): Chemical element symbol (e.g., 'C', 'O').
        position (np.ndarray): 3D coordinates of the atom.
        radius (float): van der Waals radius of the atom.
        asa (float): Solvent-accessible surface area (in Å²), default is 0.0.
        index (int): Atom index (optional, useful for tracking).
    """
    def __init__(self, element:str,  position, index=None):
        self.element = element
        self.position = np.array(position)
        self.index = index
        self.radius = VAN_DER_WAALS_RADII.get(element, 1.5)  # default if unknown
        self.asa = 0.0  # will be filled after ASA calculation

    def distance_to(self, other):
        """Compute Euclidean distance to another Atom."""
        return np.linalg.norm(self.position - other.position)

    def generate_sphere(self,n:int):
        """
        Generate a quasi-uniformes n points sphere based on
        Saff and Kuijlaars algorithm (1997).

        Args:
            n (int): Number of points to generate the spphere.

        Returns:
            points (ndarray): An array (n, 3) with coordinates x, y, z of every points.
        """

        points = np.zeros((n, 3))

        for k in range(1, n + 1):
            h = -1 + 2 * (k - 1) / (n - 1)  # Hauteur du point
            theta = np.arccos(h)            # Colatitude
            phi = np.pi * (1 + np.sqrt(5)) * (k - 1)  # Longitude (angle d'or)

            # Transfert Coordonnées sphériques vers cartésiennes
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)

            points[k - 1] = np.array([x, y, z])

        return points

    # def __repr__(self):
    #     return f"Atom({self.index}, {self.element}, ASA={self.asa:.2f})"

    def __str__(self):
        return f"Atom({self.index}, {self.element}, ASA = {self.asa:.2f})"


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
        for line in pdb_file:
            if line.startswith("ATOM"):
                # print(line.strip().split())
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
                print(atom.__dict__)
    print("Done")
