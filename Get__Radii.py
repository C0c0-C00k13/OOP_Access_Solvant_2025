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

def set_dict_atoms(filename:str):
    VAN_DER_WAALS_RADII = {}
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

                    residue = integrate_line.strip().split()[-2]
                    nb_atoms = integrate_line.strip().split()[-1]
                    # print(VAN_DER_WAALS_RADII)
                    # print(residue, nb_atoms)
                    VAN_DER_WAALS_RADII[residue] = {'nb_atoms' : nb_atoms}
                    # print(VAN_DER_WAALS_RADII)
                    # time.sleep(1)

                # Atoms
                elif integrate_line.startswith("ATOM"):
                    # print(integrate_line.strip().split())
                    if is_heteroatom and integrate_line.strip().split()[1] == "N":
                        atom_i = "".join(integrate_line.strip().split()[1:3])
                        radius = float(integrate_line.strip().split()[3])
                        polarity = int(integrate_line.strip().split()[4])
                        # print(atom_i, radius, polarity)
                    else:
                        atom_i = integrate_line.strip().split()[1]
                        radius = float(integrate_line.strip().split()[2])
                        polarity = int(integrate_line.strip().split()[3])

                    values = ({
                            'radius' : radius,
                            'polarity' : polarity 
                        })
                    VAN_DER_WAALS_RADII[residue].update({atom_i : values})
    return VAN_DER_WAALS_RADII

def get_reference_total_asa(filename:str, probe_radius:float=1.4):
    # Checking if the file exists
    IS_EXIST = os.path.exists(filename)
    if IS_EXIST:
        print(f"Radii references: '{filename}' found...")
        time.sleep(2)

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
                            radius = float(integrate_line.strip().split()[3]) + probe_radius
                        else:
                            radius = float(integrate_line.strip().split()[2]) + probe_radius
                        if residue_name in max_asa:
                            # print(f"previous value:{max_asa[residue_name]}")
                            max_asa[residue_name] += 4 * math.pi * radius**2
                            # print(f"New value:{max_asa[residue_name]}")
                        else:
                            max_asa[residue_name] = 4 * math.pi * radius**2
                            # print(f"New space:{max_asa[residue_name]}")
            # Rounding the ASA to 3 decimals
        for residue in max_asa:
            max_asa[residue] = round(max_asa[residue], 3)
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

def max_sasa(radius, probe_radius=1.4):
    """Calculate the max solvent-accessible surface area of an atom."""
    return 4 * math.pi * (radius + probe_radius) ** 2

def parse_and_compute_sasa(file_path):
    results = []
    current_residue = None
    current_atoms = []

    with open(file_path, 'r') as f:
        is_het = False
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            elif line.startswith('RESIDUE'):
                if "HETATM" in line:
                    is_het = True
                else:
                    is_het = False
                # Save previous residue
                if current_residue and current_atoms:
                    total_sasa = 0.0
                    polar_sasa = 0.0
                    nonpolar_sasa = 0.0

                    for radius, polarity in current_atoms:
                        sasa = max_sasa(radius)
                        total_sasa += sasa
                        if polarity == 1:
                            polar_sasa += sasa
                        else:
                            nonpolar_sasa += sasa

                    results.append((current_residue, total_sasa, polar_sasa, nonpolar_sasa))

                # Start new residue
                parts = line.split()
                resname = parts[2]
                resid = parts[3]
                current_residue = f"{resname} {resid}"
                current_atoms = []

            elif line.startswith('ATOM'):
                if is_het and line.split()[1] == "N":
                    parts = line.split()
                    radius = float(parts[3])
                    polarity = int(parts[4])
                    current_atoms.append((radius, polarity))
                else:
                    parts = line.split()
                    radius = float(parts[2])
                    polarity = int(parts[3])
                    current_atoms.append((radius, polarity))

        # Handle last residue
        if current_residue and current_atoms:
            total_sasa = 0.0
            polar_sasa = 0.0
            nonpolar_sasa = 0.0

            for radius, polarity in current_atoms:
                sasa = max_sasa(radius)
                total_sasa += sasa
                if polarity == 1:
                    polar_sasa += sasa
                else:
                    nonpolar_sasa += sasa

            results.append((current_residue, total_sasa, polar_sasa, nonpolar_sasa))

    return results



if __name__ == "__main__":
    # Dictionary of van der Waals radii (in Ångströms)
    FILENAME = "./Data/vdw.radii"
    VAN_DER_WAALS_RADII = get_radii(FILENAME)
    time.sleep(2)

    for atom in VAN_DER_WAALS_RADII:
        radius = get_radius(VAN_DER_WAALS_RADII, atom)
        print(atom,":",radius)

    max_axa = get_reference_total_asa(FILENAME)
    for residue in max_axa:
        max_axa[residue] = round(max_axa[residue], 3)
    print(max_axa)

    time.sleep(2)

    sasa_per_residue = parse_and_compute_sasa(FILENAME)

    # Print results
    print(f"{'Residue':<10} {'Total_SASA(Å²)':>15} \
{'Polar_SASA(Å²)':>15} {'Nonpolar_SASA(Å²)':>20}")
    for res, total, polar, nonpolar in sasa_per_residue:
        print(f"{res:<10} {total:>15.2f} {polar:>15.2f} {nonpolar:>20.2f}")

    dict_atoms = set_dict_atoms(FILENAME)
    # print(dict_atoms)
