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
import numpy as np

# Dictionary of van der Waals radii (in Ångströms)
VAN_DER_WAALS_RADII = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8
}

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
    def __init__(self, element, position, index=None):
        self.element = element
        self.position = np.array(position)
        self.index = index
        self.radius = VAN_DER_WAALS_RADII.get(element, 1.5)  # default if unknown
        self.asa = 0.0  # will be filled after ASA calculation

    def distance_to(self, other):
        """Compute Euclidean distance to another Atom."""
        return np.linalg.norm(self.position - other.position)

    def __repr__(self):
        return f"Atom({self.index}, {self.element}, ASA={self.asa:.2f})"



