"""Retrieves atoms radii from file 'vdw.radii'"""
import time
import os

# Defaults Radii
VAN_DER_WAALS_RADII = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8
}


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

    atoms = {}
    with open(filename, 'r') as radii_file :
        for line in radii_file:
            integrate_line = radii_file.readline()

            if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
                # New residu

                if integrate_line.startswith("ATOM"):
                    atom_i = integrate_line.strip().split()[1]
                    if atom_i == "SD" or atom_i == "SG":
                        atom_i = "S"
                    if atom_i not in atoms:
                        # Atom : Radius
                        atoms[atom_i] = float(integrate_line.strip().split()[2])
    return atoms


def get_radius(radii_reference, element:str)->float:
    """Returns the radius of an element"""
    return radii_reference[element]


if __name__ == "__main__":
    # Dictionary of van der Waals radii (in Ångströms)
    FILENAME = "./Data/vdw.radii"

    # Checking if the file exists
    IS_EXIST = os.path.exists(FILENAME)
    if IS_EXIST:
        print(f"Radii references: '{FILENAME}' found...")
        time.sleep(2)
    #     VAN_DER_WAALS_RADII = get_radii(FILENAME)
    # else:
    #     print(f"This file does not exist. Default values:\n{VAN_DER_WAALS_RADII}.")
    # time.sleep(2)

    # for atom in VAN_DER_WAALS_RADII:
    #     radius = get_radius(VAN_DER_WAALS_RADII, atom)
    #     print(atom,":",radius)

    list_residues = []
    residues = {}
    atom = {}
    atoms = []
    with open(FILENAME, 'r') as radii_file :
        # ILE_STATEMENT = False
        residue, atom_name, radius, polarity = None, None, None, None
        for line in radii_file:
            # Ignore Header/ Footer and Empty line
            if not line.startswith("#") and not line.startswith("\n"):
                if line.startswith("RESIDUE"):
                    list_residues.append(residues)
                    residues = {residue : atoms}
                    print(f"{residues}\n")
                    time.sleep(1)
                    residue = line.strip().split()[2]
                    atom = {}
                    atoms = []
                    time.sleep(1)
                if line.startswith("ATOM"):
                    atom_name = line.strip().split()[1]
                    radius = float(line.strip().split()[2])
                    polarity = int(line.strip().split()[3])
                    atom[atom_name] = {'polarity' : polarity, 'radius' : radius}
                    atoms.append(atom)
