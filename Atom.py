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
        points (np.array): 3D coordinates of the points of the atom
        self.exposed_points (int): number of points exposed to solvant
        asa (float): Solvent-accessible surface area (in Å²), default is 0.0.
        index (int): Atom index (optional, useful for tracking).
    """
    def __init__(self, element:str,  position, index=None):
        self.element = element
        self.position = np.array(position)
        self.index = index
        self.radius = Get__Radii.get_radius(VAN_DER_WAALS_RADII, element)  # default if unknown
        self.points = None # Will be filled after Sphere generation
        self.exposed_points = None # Will be filled by distance calculation
        self.asa = 0.0  # will be filled after ASA calculation

    # def __repr__(self):
    #     return f"Atom({self.index}, {self.element}, ASA={self.asa:.2f})"

    def __str__(self):
        return f"Atom({self.index}, {self.element}, ASA = {self.asa:.2f} Å², radius={self.radius})"


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
                position = (coord_x,coord_y,coord_z)
                # Element
                element = line.strip().split()[-1]
                # print(f"Index:{index}; Position:{position}; Element:{element}")
                atom = Atom(element=element,position=position,index=index)
                # print(atom.__dict__)
                print(atom)

    print("Done")
