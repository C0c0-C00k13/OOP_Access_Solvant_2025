import math
import time

FILENAME_RADIUS = "./Data/vdw.radii"
FILENAME_PROTEIN = "./Data/2c8r.pdb"
RADIUS_PROBE = 1.4

# ======================================
# STEP 1 : READ FILES
# ======================================

# Read PDB file
def atoms_descriptions(filename_pdb, filename_raddi):
    """"""
    radii_dict = simplified_radii(filename_raddi)
    list_atoms = {}
    with open(filename_pdb,"r") as f_in:
        for line in f_in:
            if line.startswith("ATOM"):
                id_atom = int(line.strip().split()[1])
                [x, y, z] = line.strip().split()[-6:-3]
                list_atoms[id_atom] = {
                    "element" : line.strip().split()[2],
                    "residue" : line.strip().split()[3],
                    "num_residue" : int(line.strip().split()[5]),
                    "chain" : line.strip().split()[4],
                    "coordinates" : (float(x),float(y),float(z)),
                    "radius": radii_dict[line.strip().split()[2]]['radius'],
                    "polarity" : radii_dict[line.strip().split()[2]]['polarity']
                }
    return list_atoms


# Read Radius file
def read_radii(filename):
    atoms = {}
    with open(filename, 'r') as radii_file :
        is_heteroatom = False
        for integrate_line in radii_file:
            # New residu
            if integrate_line.startswith("RESIDUE"):
                if integrate_line.strip().split()[1] == "HETATM":
                    is_heteroatom = True
                else:
                    is_heteroatom = False
                residue = integrate_line.strip().split()[2]
                atoms[residue] = {}

            # New atom
            if integrate_line.startswith("ATOM"):
                atom, radius, polarity = line_atom(integrate_line, is_heteroatom)
                atom_desc = {'radius' : radius, 'polarity' : polarity}
                atoms[residue].update({atom : atom_desc})
    return atoms


# Simplified version
def simplified_radii(filename):
    atoms = {}
    with open(filename, 'r') as radii_file :
        is_heteroatom = False
        for integrate_line in radii_file:
            # New residu
            if integrate_line.startswith("RESIDUE"):
                if integrate_line.strip().split()[1] == "HETATM":
                    is_heteroatom = True
                else:
                    is_heteroatom = False

            # New atom
            if integrate_line.startswith("ATOM"):
                atom, radius, polarity = line_atom(integrate_line, is_heteroatom)
                if atom in atoms:
                    continue
                atoms[atom] =  {'radius' : radius, 'polarity' : polarity}
    return atoms


def line_atom(line, is_hetatm):
    if is_hetatm and line.strip().split()[1] == "N":
        atom_i = "_".join(line.strip().split()[1:3])
    else:
        atom_i = line.strip().split()[1]

    radius = float(line.strip().split()[-2])
    polarity = int(line.strip().split()[-1])
    return atom_i, radius, polarity

# ======================================
# STEP 2 : CALCULATE SPHERE + LIST OF POINTS
# ======================================

def saff_kuijlaars_points(n, id_atom, center=(0.0,0.0,0.0), radius=1.0):
    """
    Generates a n-points quasi-uniform sphere based on Saff and Kuijlaars algorithm.

    Parameters
    ---
    n (int) : Number of points to generate.
    id_atom : Position of the atom in the PDB file.
    center (tuple) : Coordinates of the center of the sphere. Default value : (0.0,0.0,0.0).
    radius (float) : Radius of the atom. Default value : 0.0.

    Returns
    ---
    points ([tuple]): A list with an atom in the pdb file and  x, y, z of a point.
    """

    points = []

    for k in range(1, n + 1):
        h = -1 + 2 * (k - 1) / (n - 1)  # Hauteur du point
        theta = math.acos(h)            # Colatitude
        if k == 1 or k == n:
            phi = 0.0
        else:
            phi += 3.6 / math.sqrt(n * (1 - h * h))

        # Coordonnées sphériques vers cartésiennes
        x = round((center[0] + math.sin(theta) * math.cos(phi) * radius), 3)
        y = round((center[1] + math.sin(theta) * math.sin(phi) * radius), 3)
        z = round((center[2] + math.cos(theta) * radius), 3)

        points.append( (id_atom, (x, y, z)) )

    return points

# ======================================
# ACCESSIBLE POINT TEST
# ======================================

# ======================
# DISTANCE
# ======================

def distance(point_a, point_b)->float:
    """Calculates distance between 2 points
    
    Parameters
    ---
    point_a : Coordinates from either a point of a an atom or an atom.
    point_b : Coordinates from either a point of a an atom or an atom.
    
    Returns
    ---
    distance (float) : Distance from 2 points (a, b).
    """
    # return math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3)))
    return round(math.sqrt(sum((point_a[i]-point_b[i])**2 for i in range(3))), 3)

# ======================
# ACCESSIBLE POINT TEST
# ======================

def point_exposed(point, list_point, probe):
    for current_point in list_point:
        if current_point[0] is point[0]:
            continue
        if distance(point[1],current_point[1]) < 2 * probe:
            # print(point, current_point, distance(point[1],current_point[1]))
            # time.sleep(2)
            return False
    # print(point, current_point, distance(point[1],current_point[1]))
    # time.sleep(2)
    return True

def calculate_asa(atom, nb_point, list_exposed):
    """"""
    # -- ASA
    # Angstrom value
    total_sphere_area = 4 * math.pi * atom['radius']**2
    # Sphere surface per point
    unit_sphere_area = (4 * math.pi * atom['radius']**2)/ nb_point

    # -- RSA
    # Angstrom value
    exposed_surface = unit_sphere_area * list_exposed[atom]
    # Percentage
    relative_surface = 100 * (4 * math.pi * atom['radius']**2)/ nb_point
    
    # return list_exposed

# ======================
# MAIN
# ======================


print("READING FILES...")
extracted_atoms = atoms_descriptions(FILENAME_PROTEIN, FILENAME_RADIUS)
print("READING FILES - DONE")
list_points = []

for item in extracted_atoms.items():
    # print(item)
    list_points += saff_kuijlaars_points(n=92, id_atom=item[0], center=item[1]['coordinates'],radius=item[1]['radius'])
# print(list_points)
print("Done")
list_exposed = {}
start = time.time()
for central_point in list_points:
    if central_point[0] not in list_exposed:
        list_exposed[central_point[0]] = 0
    # print(central_point)
    exposed = point_exposed(point=central_point,list_point=list_points,probe=1.4)
    # print(exposed)
    # time.sleep(2)
    if exposed:
        list_exposed[central_point[0]] += 1
        print(list_exposed)
        time.sleep(2)
# time.sleep(4)
print(len(list_points),len(list_exposed))
print(list_exposed)

end = time.time()
print(f"calculus duration = {(end - start)} ; {math.floor((end-start)/60)}' {round((end - start)%60,3)}")
