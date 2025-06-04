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

    # residues = {}
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
                    if not atom_i in atoms:
                        # Atom : Radius
                        atoms[atom_i] = integrate_line.strip().split()[2]
    return atoms

    #             if integrate_line.startswith("RESIDUE"):
    #                 residue = integrate_line.strip().split()[2]
    #                 residues[residue] = []
    #             elif integrate_line.startswith("ATOM"):
    #                 residues[residue].append({
    #                     'name' : integrate_line.strip().split()[1], 
    #                     'radius' : integrate_line.strip().split()[2],
    #                     'polar' : integrate_line.strip().split()[3]
    #                 })
    # return residues         

def get_radius(radii_reference, element:str):
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
        VAN_DER_WAALS_RADII = get_radii(FILENAME)
    else:
        print(f"This file does not exist. Default values:\n{VAN_DER_WAALS_RADII}.")
    time.sleep(2)

    # print(residues)
    # print(residues.keys())
    for atom in VAN_DER_WAALS_RADII.keys():
        radius = get_radius(VAN_DER_WAALS_RADII, atom)
        print(atom,":",radius)