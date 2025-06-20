import math
import os
import argparse
import datetime
from collections import defaultdict
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
# create logger with '__name__'
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('./Logs/main.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

# ===========================
# CONFIGURATION
# ===========================

PROBE_RADIUS = 1.4  # Radius of water probe (O radius in Å)
POINTS_PER_SPHERE = 92  # More points = higher accuracy, slower
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
    parser.add_argument("-n", "--points", type=int, default=POINTS_PER_SPHERE,
                        help=f"Number of points representing the sphere, default: {PROBE_RADIUS}")
    parser.add_argument("-o", "--output", type=str, default="output",
                        help="Name of the output file, default : output")
    parser.add_argument("-p", "--probe", type=float, default=PROBE_RADIUS,
                        help="Probe radius (default: 1.4)")
    parser.add_argument("-r", "--radii", type=str, default=None,
                        help="Custom radii file (default: None)")

    return parser.parse_args()


# ===========================
# STEP 1: READ PDB FILE
# ===========================

def read_pdb(filename, add_hetatm):
    """Returns the atoms from the PDB file.
    
    Parameters
    ---
    filename (str) : Name (or Path) of the PDB file.
    include_hetatm (bool) : Determines if the function returns heteroatoms.

    Returns
    ---
    atoms ([Dict]) : List of atoms.
    heteroatoms ([Dict]) : List of atoms.
    """

    atoms = []
    heteroatoms = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atoms.append(parse_pdb_line(line=line))
            if add_hetatm:
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
                atom_name, radius, polarity = parse_line_radii(line, is_heteroatom)
                max_residue_asa[residue_name] += calculate_max_atom_surface(radius)
                if polarity == 1:
                    polar_atoms.add(atom_name)
                if atom_name in atoms:
                    continue
                atoms.update({atom_name : radius})

    return atoms, polar_atoms, max_residue_asa


def parse_line_radii(line, is_hetatm):
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

def is_point_exposed(px, py, pz, atoms, this_atom,
                    probe=PROBE_RADIUS, vdw_radii=VDW_RADII):
    """ Returns a False if there is an overlap between to points of different atoms.
    
    Parameters
    ---
    px (float) : Point coordinate on the x-axis.
    py (float) : Point coordinate on the y-axis.
    pz (float) : Point coordinate on the z-axis.
    atoms (List) : List of of atoms.
    this_atom (Dict): Atom linked to the current point. This atom will be skipped when comparing distances. 
    probe (float) : Radius of probe.
    vdw_radii (Dict) : Listt of radius of every atom.

    Returns
    ---
    True/ False (boolean)    
    """

    for atom in atoms:
        if atom is this_atom:
            continue
        ex, ey, ez = atom["x"], atom["y"], atom["z"]
        r = vdw_radii.get(atom["atom_name"], 1.7) + probe
        dx, dy, dz = px - ex, py - ey, pz - ez
        distance = dx*dx + dy*dy + dz*dz
        r_square = r**2
        logger.debug(msg=f"Radius of current 'other' atom : {r} ; Square {r_square} ; Distance {distance}")
        if distance <r_square:
            logger.debug(msg=f"False for atom {this_atom}, point {(px, py, pz)} compared with atom  {atom}")
            return False
    
    return True

# ===========================
# STEP 4: CALCULATE ASA/RSA
# ===========================

def calculate_asa(atoms, hetatoms, probe=PROBE_RADIUS,
                point_per_sphere=POINTS_PER_SPHERE, vdw_radii=VDW_RADII,
                polar_list=POLAR_ELEMENTS, main_chain_elements=MAIN_CHAIN_ATOMS):
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

    sphere = generate_sphere_points(point_per_sphere)
    point_area = 4 * math.pi / point_per_sphere
    res_asa = {}
    chain_stats = defaultdict(lambda: {
        "main": 0.0, "side": 0.0,
        "polar": 0.0, "apolar": 0.0,
        "total": 0.0
    })

    total_atom = atoms + hetatoms
    for atom in atoms:
        x, y, z = atom["x"], atom["y"], atom["z"]
        element = atom["element"]
        atom_name = atom["atom_name"]
        res = atom["res_name"]
        chain = atom["chain_id"]
        res_id = atom["res_id"]
        key = (chain, res_id, res)

        r = vdw_radii.get(atom_name, 1.7) + probe
        exposed_points = 0

        for dx, dy, dz in sphere:
            px = x + r * dx
            py = y + r * dy
            pz = z + r * dz
            if is_point_exposed(px=px, py=py, pz=pz,atoms=total_atom,
                                this_atom=atom, probe=probe, vdw_radii=vdw_radii):
                exposed_points += 1

        atom_asa = exposed_points * point_area * (r ** 2)
        if key not in res_asa:
            res_asa[key] = {
                "total": 0.0,
                "polar": 0.0, "apolar": 0.0,
                "main": 0.0, "side": 0.0
            }
        res_asa[key]["total"] += atom_asa

        chain_stats[chain]["total"] += atom_asa

        # Discriminates MAIN chain an SIDE chain
        if atom_name in main_chain_elements:
            chain_stats[chain]["main"] += atom_asa
            res_asa[key]["main"] += atom_asa
        else:
            chain_stats[chain]["side"] += atom_asa
            res_asa[key]["side"] += atom_asa

        # Discriminates POLAR elements ande NON-POLAR elements
        if element in polar_list:
            res_asa[key]["polar"] += atom_asa
            chain_stats[chain]["polar"] += atom_asa
        else:
            res_asa[key]["apolar"] += atom_asa
            chain_stats[chain]["apolar"] += atom_asa

    return res_asa, chain_stats



# ===========================
# STEP 5: WRITE OUTPUT
# ===========================

def write_output(filename, res_asa, chain_stats, max_asa=MAX_ASA):
    # with open(filename, 'w') as f:
    print("\nResidue ASA and RSA:\n")
    print(f"{'Chain':<5} {'ResID':<6} {'ResName':<7} "
            f"{'TotalASA':<10} {'RSA(%)':<8} "
            f"{'PolarASA':<10} {'PolarRSA':<10} "
            f"{'ApolarASA':<11} {'ApolarRSA':<10}\n")
    for (chain, res_id, res_name), data in sorted(res_asa.items()):
        max_ref = max_asa.get(res_name, 200)
        total = data["total"]
        # main_chain = data["main"]
        # side_chain = data["side"]
        polar = data["polar"]
        apolar = data["apolar"]
        rsa_total = (total / max_ref) * 100
        # rsa_main = (main_chain / max_ref) * 100
        # rsa_side = (side_chain / max_ref) * 100
        rsa_polar = (polar / max_ref) * 100
        rsa_apolar = (apolar / max_ref) * 100
        print(f"{chain:<5} {res_id:<6} {res_name:<7} "
                f"{total:<10.2f} {rsa_total:<8.2f} "
                # f"{main_chain:<10.2f} {rsa_main:<10.2f} "
                # f"{side_chain:<10.2f} {rsa_side:<10.2f} "
                f"{polar:<10.2f} {rsa_polar:<10.2f} "
                f"{apolar:<11.2f} {rsa_apolar:<10.2f}\n")
    print("\nPer-Chain ASA Summary:\n")
    print(f"{'Chain':<5} {'Main ASA':<12} {'Side ASA':<12} {'Polar ASA':<12} {'Apolar ASA':<12} {'Total ASA':<12}\n")
    for chain, stats in sorted(chain_stats.items()):
        print(f"{chain:<5} {stats['main']:<12.2f} {stats['side']:<12.2f} "
                f"{stats['polar']:<12.2f} {stats['apolar']:<12.2f} {stats['total']:<12.2f}\n")

# ===========================
# MAIN
# ===========================

def main():
    logging.basicConfig(filename='./Logs/SASA.log', level=logging.INFO, filemode='w')
    logger.info(msg='Started')
    today = datetime.datetime.now().strftime("%m-%d-%Y %H:%M:%S")
    print(today)
    args = parse_args()

    pdb_file = args.pdb_file
    include_hetatm = args.hetero == "y"
    point_per_sphere = args.points
    probe_radius = args.probe
    custom_radii_file = args.radii
    output = args.output

    print(f"PDB file         : {pdb_file}")
    print(f"Include HETATM   : {include_hetatm}")
    print(f"Probe radius     : {probe_radius}")
    print(f"Custom radii file: {custom_radii_file}")
    
    description_radii, polar_elements, max_axa = VDW_RADII, POLAR_ELEMENTS, MAX_ASA
    if custom_radii_file is not None:
        if os.path.exists(custom_radii_file):
            description_radii, polar_elements, max_axa = set_description(custom_radii_file)

    atoms, hetatoms = read_pdb(filename=pdb_file, add_hetatm=include_hetatm)
    
    # Residues, Chain
    res_asa, chain_stats = calculate_asa(atoms=atoms,hetatoms=hetatoms, probe=probe_radius,
                                        polar_list=polar_elements, vdw_radii=description_radii,
                                        point_per_sphere=point_per_sphere, main_chain_elements=MAIN_CHAIN_ATOMS)


    path_result_directory = f"./Results/{today.split()[0]}/{pdb_file.strip().split('/')[-1][:-4]}/SASA"
    nested_directory_path = Path(path_result_directory)
    nested_directory_path.mkdir(parents=True, exist_ok=True)
    output_filename = f"{path_result_directory}/{output}.asa"
    write_output(filename=output_filename, res_asa=res_asa, chain_stats=chain_stats, max_asa=max_axa)
    print(f"The output file '{output}.asa' has been created in the directory: {nested_directory_path}")
    print('Done')


if __name__ == "__main__":

    main()
