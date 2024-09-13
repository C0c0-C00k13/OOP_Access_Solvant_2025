"""Ceci est le fichier de la classe point_atome.
"""

import numpy as np
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Atom
from Bio.PDB import NACCESS

import Fonctions_utiles


class point_atome:
    """
    """


    def __init__(self, atom_center, x_pt, y_pt, z_pt):
        self.atom_center = atom_center 
        self.x_pt = x_pt
        self.y_pt = y_pt
        self.z_pt = z_pt

    def __str__(self):
        return f"Atome {print(self.atom_center)},\t Point coords[{self.x_pt},{self.y_pt},{self.z_pt}]"

    def calcul_distance(self, atome):
        help = "Methodes pemettant de calculer la distance entre 2 atomes"
        if isinstance(atome,Bio.PDB.Atom.Atom):
            distX = abs(self.x_pt - atome.coord[0])
            distY = abs(self.y_pt - atome.coord[1])
            distZ = abs(self.z_pt - atome.coord[2])
            return pow( (pow((distX), 2) + pow((distY),2) + pow((distZ),2)), 0.5)

if __name__ == "__main__":
    pass
