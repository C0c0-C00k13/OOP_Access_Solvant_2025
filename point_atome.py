"""Class Point of a sphere representing an Atom.
"""

import os
import sys
import time
import datetime
import numpy as np
from Atom import Atom
# import Fonctions_utiles


class PointAtom:
    """Represents the points forming a sphere around a central atom 
    """


    def __init__(self, atom_center:Atom, x_pt:float, y_pt:float, z_pt:float):
        self.atom_center = atom_center
        self.x_pt = x_pt
        self.y_pt = y_pt
        self.z_pt = z_pt

    def __str__(self):
        return f"Atome: {self.atom_center},\
 Point coords: [{self.x_pt},{self.y_pt},{self.z_pt}]"

    def calcul_distance(self, atome:Atom)->float:
        """Calculates distance between an atom and the point
        """
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        dist_x = abs(self.x_pt - atome.coord[0])
        dist_y = abs(self.y_pt - atome.coord[1])
        dist_z = abs(self.z_pt - atome.coord[2])
        return pow( (pow((dist_x), 2) + pow((dist_y),2) + pow((dist_z),2)), 0.5)

if __name__ == "__main__":
    # pass
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

        print("Generating point")
        point = PointAtom(atom_center=atom,x_pt=coord_x,y_pt=coord_y,z_pt=coord_z)
        print(point)

    print("Done")
