import math
import os
import argparse
from collections import defaultdict

# ===========================
# CONFIGURATION
# ===========================

FILENAME_RADIUS = "./Data/vdw.radii"
PROBE_RADIUS = 1.4  # Radius of water probe (O radius in Å)
POINTS_PER_SPHERE = 92  # More points = higher accuracy, slower
DEFAULT_INCLUDE_HETATM = False
VDW_RADII = {
    "H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8,
}
# Max ASA values for residues (Å²), for RSA computation
MAX_ASA = {
    "ALA": 113, "ARG": 241, "ASN": 158, "ASP": 151, "CYS": 140,
    "GLN": 189, "GLU": 183, "GLY": 85, "HIS": 194, "ILE": 182,
    "LEU": 180, "LYS": 211, "MET": 204, "PHE": 218, "PRO": 143,
    "SER": 122, "THR": 146, "TRP": 259, "TYR": 229, "VAL": 160
}

MAIN_CHAIN_ATOMS = {"N", "CA", "C", "O", "OXT", "H", "HA"}
POLAR_ELEMENTS = {"N", "O", "S"}


# ===========================
# ARGUMENT PARSER
# ===========================

def parse_args():
    """Parses the arguments provided in the command line.

    Returns
    ---
    parser.parse_args() (Namespace) : 
    """

    parser = argparse.ArgumentParser(description="Compute ASA/RSA from PDB.")

    parser.add_argument("pdb_file", type=str, help="PDB file to process")
    parser.add_argument("-i", "--hetero", choices=["y", "n"], default="n",
                        help="Include HETATM records (y/n), default: n")
    parser.add_argument("-n", "--points", type=int, default=92,
                        help="Number of points representing the sphere, default: 92")
    parser.add_argument("-p", "--probe", type=float, default=PROBE_RADIUS,
                        help="Probe radius (default: 1.4)")
    parser.add_argument("-r", "--radii", type=str, default=None,
                        help="Custom radii file (default: None)")

    return parser.parse_args()


# ===========================
# STEP 1: READ PDB FILE
# ===========================

def read_pdb(filename, include_hetatm):
    """Returns the atoms from the PDB file.
    
    Parameters
    ---
    filename (str) : Name (or Path) of the PDB file.

    Returns
    ---
    atoms ([Dict]) : List of atoms. 
    """

    atoms = []
    heteroatoms = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atoms.append(parse_pdb_line(line=line))
            if include_hetatm:
                if line.startswith("HETATM"):
                    heteroatoms.append(parse_pdb_line(line=line))
    return atoms, heteroatoms


def parse_pdb_line(line):
    """Parse a single HETATM line from a PDB file.
    
    Parameters
    ---
    line (str) : Description of an atom. A line in the format of PDB is recommended.
    
    Returns
    ---
    (Dict) : Summary of the description in the line.
    """

    return {
        "atom_serial": int(line[6:11].strip()),
        "atom_name": line[12:16].strip(),
        "res_name": line[17:20].strip(),
        "chain_id": line[21].strip(),
        "res_id": int(line[22:26].strip()),
        "x": float(line[30:38].strip()),
        "y": float(line[38:46].strip()),
        "z": float(line[46:54].strip()),
        "element": line[76:78].strip()
    }


def set_description(filename):
    """Returns the descriptions contained in the file.
    
    Parameters
    ---
    filename (str) : Name of the file access by the function.
    
    Returns
    ---
    atoms (Dict) : List of the radius of every atoms.
    polar_atoms (set) : List of the polar atoms.
    max_residue_asa (Dict) : List of the max asa of every residue.
    """

    atoms = {}
    polar_atoms = set()
    is_heteroatom = False
    max_residue_asa = {}

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("RESIDUE"):
                type = line.strip().split()[1] 
                if type == "ATOM" or type == "NUCL":
                    residue_name = line.strip().split()[-2]
                    max_residue_asa.update({residue_name:0})
                is_heteroatom = line.strip().split()[1] == "HETATM"

            elif line.startswith("ATOM"):
                atom_name, radius, polarity = line_atom(line, is_heteroatom)
                max_residue_asa[residue_name] += calculate_max_atom_surface(radius)
                if polarity == 1:
                    polar_atoms.add(atom_name)
                if atom_name in atoms:
                    continue
                atoms.update({atom_name : radius})

    return atoms, polar_atoms, max_residue_asa


def line_atom(line, is_hetatm):
    """Returns the description of an atom.
    
    Parameters
    ---
    line (str) : Description of an atom.
    is_hetatm (bool) : Determines if the current atom is a heteroatom.

    Returns
    ---
    atom_name (str) : Name of the atom described in the line.
    radius (float) : Radius of the atom described in the line.
    polarity (int) : Polarity of the atom described in the line.
    """

    if is_hetatm and line.strip().split()[1] == "N":
        atom_name = "_".join(line.strip().split()[1:3])
    else:
        atom_name = line.strip().split()[1]

    radius = float(line.strip().split()[-2])
    polarity = int(line.strip().split()[-1])
    return atom_name, radius, polarity


def calculate_max_atom_surface(radius)->float:
    """Returns the sphere surface area of an atom. 

    Parameter
    ---
    radius (float) : Radius of an atom.

    Returns
    ---
    (float) : Surface area of the atom.
    """

    return 4 * math.pi * radius**2

# ===========================
# STEP 2: SAFF & KUIJLAARS SPHERE
# ===========================

def generate_sphere_points(n):
    """
    Generates a n-points quasi-uniform sphere based on Saff and Kuijlaars algorithm.

    Parameters
    ---
    n (int) : Number of points to generate.

    Returns
    ---
    points ([tuple]): A list with an atom in the pdb file and  x, y, z of a point.
    """

    points = []
    offset = 2.0 / n
    increment = math.pi * (3.0 - math.sqrt(5))
    for k in range(n):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        points.append((x, y, z))
    return points

# ===========================
# STEP 3: ACCESSIBLE POINT TEST
# ===========================

def is_point_exposed(px, py, pz, atoms, this_atom, probe, hetatoms):
    """ Returns a False if there is an overlap between to point s of different atoms.
    
    Parameters
    ---
    px (float) : Point coordinate on the x-axis.
    py (float) : Point coordinate on the y-axis.
    pz (float) : Point coordinate on the z-axis.
    atoms (List) : List of of atoms.
    this_atom (Dict): Atom linked to the current point. This atom will ba skipped when comparing distances. 
    probe (float) : Radius of probe.

    Returns
    ---
    True/ False (boolean)    
    """

    for atom in atoms:
        if atom is this_atom:
            continue
        ex, ey, ez = atom["x"], atom["y"], atom["z"]
        r = VDW_RADII.get(atom["element"], 1.7) + probe
        dx, dy, dz = px - ex, py - ey, pz - ez
        if dx*dx + dy*dy + dz*dz < r * r:
            return False
    return True

# ===========================
# STEP 4: CALCULATE ASA/RSA
# ===========================

def calculate_asa(atoms, hetatoms, probe=PROBE_RADIUS):
    """ Returns the ASA (Accessible Solvant Area) of each residue and chain.

    Parameters
    ---
    atoms ([Dict]) : List of atoms to process to calculate ASA.
    probe (float) : Radius of the probe. Default value is set to constant PROBE_RADIUS.

    Returns
    ---
    res_asa (Dict) : List of ASA for every residue of the list.
    chain_stat (Dict) : List of ASA of each chain of the protein.
    """

    sphere = generate_sphere_points(POINTS_PER_SPHERE)
    point_area = 4 * math.pi / POINTS_PER_SPHERE
    res_asa = {}
    chain_stats = defaultdict(lambda: {
        "main": 0.0, "side": 0.0,
        "polar": 0.0, "apolar": 0.0,
        "total": 0.0
    })

    for atom in atoms:
        x, y, z = atom["x"], atom["y"], atom["z"]
        element = atom["element"]
        atom_name = atom["atom_name"]
        res = atom["res_name"]
        chain = atom["chain_id"]
        res_id = atom["res_id"]
        key = (chain, res_id, res)

        r = VDW_RADII.get(element, 1.7) + probe
        exposed_points = 0

        for dx, dy, dz in sphere:
            px = x + r * dx
            py = y + r * dy
            pz = z + r * dz
            if is_point_exposed(px=px, py=py, pz=pz,
                                atoms=atoms, this_atom=atom, probe=probe, hetatoms=hetatoms):
                exposed_points += 1

        atom_asa = exposed_points * point_area * (r ** 2)
        if key not in res_asa:
            res_asa[key] = {"total": 0.0, "polar": 0.0, "apolar": 0.0}
        res_asa[key]["total"] += atom_asa

        chain_stats[chain]["total"] += atom_asa

        # Discriminates MAIN chain an SIDE chain
        if atom_name in MAIN_CHAIN_ATOMS:
            chain_stats[chain]["main"] += atom_asa
        else:
            chain_stats[chain]["side"] += atom_asa

        # Discriminates POLAR elements ande NON-POLAR elements
        if element in POLAR_ELEMENTS:
            res_asa[key]["polar"] += atom_asa
            chain_stats[chain]["polar"] += atom_asa
        else:
            res_asa[key]["apolar"] += atom_asa
            chain_stats[chain]["apolar"] += atom_asa

    return res_asa, chain_stats

# ===========================
# MAIN
# ===========================

def main():
    args = parse_args()

    pdb_file = args.pdb_file
    include_hetatm = args.hetero == "y"
    probe_radius = args.probe
    custom_radii_file = args.radii

    print(f"PDB file         : {pdb_file}")
    print(f"Include HETATM   : {include_hetatm}")
    print(f"Probe radius     : {probe_radius}")
    print(f"Custom radii file: {custom_radii_file}")
    
    # if custom_radii_file is not None:
    #     if os.path.isfile(custom_radii_file):
    #         VDW_RADII, POLAR_ELEMENTS, MAX_ASA = set_description(custom_radii_file)

    atoms, hetatoms = read_pdb(filename=pdb_file, include_hetatm=include_hetatm)
    
    # Residues, Chain
    res_asa, chain_stats = calculate_asa(atoms=atoms,hetatoms=hetatoms, probe=PROBE_RADIUS)
    print("\nResidue ASA and RSA:")
    print(f"{'Chain':<5} {'ResID':<6} {'ResName':<7} "
          f"{'TotalASA':<10} {'RSA(%)':<8} "
          f"{'PolarASA':<10} {'PolarRSA':<10} "
          f"{'ApolarASA':<11} {'ApolarRSA':<10}")

    for (chain, res_id, res_name), data in sorted(res_asa.items()):
        max_ref = MAX_ASA.get(res_name, 200)
        total = data["total"]
        polar = data["polar"]
        apolar = data["apolar"]

        rsa_total = (total / max_ref) * 100
        rsa_polar = (polar / max_ref) * 100
        rsa_apolar = (apolar / max_ref) * 100

        print(f"{chain:<5} {res_id:<6} {res_name:<7} "
              f"{total:<10.2f} {rsa_total:<8.2f} "
              f"{polar:<10.2f} {rsa_polar:<10.2f} "
              f"{apolar:<11.2f} {rsa_apolar:<10.2f}")

    print("\nPer-Chain ASA Summary:")
    print(f"{'Chain':<5} {'Main ASA':<12} {'Side ASA':<12} {'Polar ASA':<12} {'Apolar ASA':<12} {'Total ASA':<12}")
    for chain, stats in sorted(chain_stats.items()):
        print(f"{chain:<5} {stats['main']:<12.2f} {stats['side']:<12.2f} "
              f"{stats['polar']:<12.2f} {stats['apolar']:<12.2f} {stats['total']:<12.2f}")

if __name__ == "__main__":
    
    main()
