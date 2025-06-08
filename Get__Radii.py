"""Retrieves atoms radii from file 'vdw.radii'"""
import time
import os


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


if __name__ == "__main__":
    # Dictionary of van der Waals radii (in Ångströms)
    FILENAME = "./Data/vdw.radii"
    VAN_DER_WAALS_RADII = get_radii(FILENAME)
    time.sleep(2)

    for atom in VAN_DER_WAALS_RADII:
        radius = get_radius(VAN_DER_WAALS_RADII, atom)
        print(atom,":",radius)
