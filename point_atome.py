"""Class Point of a sphere representing an Atom.
"""

import os
import sys
import time
import datetime
import numpy as np
import Atom


import Fonctions_utiles


class point_atome:
    """
    """


    def __init__(self, atom_center:Atom, x_pt:float, y_pt:float, z_pt:float):
        self.atom_center = atom_center 
        self.x_pt = x_pt
        self.y_pt = y_pt
        self.z_pt = z_pt

    def __str__(self):
        return f"Atome {print(self.atom_center)},\t Point coords[{self.x_pt},{self.y_pt},{self.z_pt}]"

    def calcul_distance(self, atome:Atom):
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        distX = abs(self.x_pt - atome.coord[0])
        distY = abs(self.y_pt - atome.coord[1])
        distZ = abs(self.z_pt - atome.coord[2])
        return pow( (pow((distX), 2) + pow((distY),2) + pow((distZ),2)), 0.5)

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
