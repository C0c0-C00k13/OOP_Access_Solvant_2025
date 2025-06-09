"""Retrieves atoms radii from file 'vdw.radii'"""
import time
import os
import math


def get_radii(filename: str):
    """
    Returns the radii atoms based on each type of residue.
    ---
    Attributes
        filename:str
    ---
    Returns
        residues:Dict
    """

    # Checking if the file exists
    IS_EXIST = os.path.exists(filename)
    if IS_EXIST:
        print(f"Radii references: '{filename}' found...")
        time.sleep(2)

        atoms = {}
        with open(filename, 'r') as radii_file :
            is_heteroatom = False
            for integrate_line in radii_file:
                if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
                    # New residu
                    if integrate_line.startswith("RESIDUE"):
                        if integrate_line.strip().split()[1] == "HETATM":
                            is_heteroatom = True
                        else:
                            is_heteroatom = False
                    if integrate_line.startswith("ATOM"):
                        atom_i = integrate_line.strip().split()[1]
                        if is_heteroatom and atom_i == "N":
                            atom_i = "".join(integrate_line.strip().split()[1:3])
                            # print(atom_i)
                        if atom_i not in atoms:
                            # Atom : Radius
                            if is_heteroatom and atom_i[0] == "N":
                                atoms[atom_i] = float(integrate_line.strip().split()[3])
                            else:
                                atoms[atom_i] = float(integrate_line.strip().split()[2])
    else:
        print(f"This file does not exist. Returns default values.")
        # Defaults Radii
        atoms = {
            'H': 1.2,
            'C': 1.7,
            'N': 1.55,
            'O': 1.52,
            'S': 1.8
        }
    return atoms


def get_radius(radii_reference, element:str)->float:
    """Returns the radius of an element"""
    return radii_reference[element]


def get_reference_total_asa(filename:str):
    # Checking if the file exists
    IS_EXIST = os.path.exists(filename)
    if IS_EXIST:
        print(f"Radii references: '{filename}' found...")
        time.sleep(2)

        PROBE_RADIUS = 1.4
        max_asa = {}
        residue_name = ""
        with open(filename, 'r') as radii_file :
            is_heteroatom = False
            for integrate_line in radii_file:
                if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
                    # New residu
                    if integrate_line.startswith("RESIDUE"):
                        if integrate_line.strip().split()[1] == "HETATM":
                            is_heteroatom = True
                        else:
                            is_heteroatom = False
                        # Get residu name
                        # print(f'New residue: current dict:{max_asa}')
                        residue_name = integrate_line.strip().split()[2]
                        # print(integrate_line.strip())
                        # time.sleep(4)
                    elif integrate_line.startswith("ATOM"):
                        if integrate_line.strip().split()[1] == "N" and is_heteroatom:
                            radius = float(integrate_line.strip().split()[3]) + PROBE_RADIUS
                        else:
                            radius = float(integrate_line.strip().split()[2]) + PROBE_RADIUS
                        if residue_name in max_asa:
                            # print(f"previous value:{max_asa[residue_name]}")
                            max_asa[residue_name] += 4 * math.pi * radius**2
                            # print(f"New value:{max_asa[residue_name]}")
                        else:
                            max_asa[residue_name] = 4 * math.pi * radius**2
                            # print(f"New space:{max_asa[residue_name]}")
    else:
        print(f"This file does not exist. Returns default values.")
        # Max ASA values (Tien et al. 2013)
        max_asa = {
            'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0,
            'C': 148.0, 'Q': 214.0, 'E': 214.0, 'G': 97.0,
            'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0,
            'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0,
            'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0
        }
    return max_asa


if __name__ == "__main__":
    # Dictionary of van der Waals radii (in Ångströms)
    FILENAME = "./Data/vdw.radii"
    VAN_DER_WAALS_RADII = get_radii(FILENAME)
    time.sleep(2)

    for atom in VAN_DER_WAALS_RADII:
        radius = get_radius(VAN_DER_WAALS_RADII, atom)
        print(atom,":",radius)

    max_axa = get_reference_total_asa(FILENAME)
    print(max_axa)
