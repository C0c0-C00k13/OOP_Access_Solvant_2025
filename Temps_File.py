import math

FILENAME_RADIUS = "./Data/vdw.radii"
FILENAME_PROTEIN = "./Data/2c8r.pdb"
RADIUS_PROBE = 1.4

# ======================
# READ FILES
# ======================

# Read PDB file
def read_pdb(filename_pdb):
    list_atoms = {}
    with open(filename_pdb,"r") as f_in:
        for line in f_in:
            if line.startswith("ATOM"):
                # print(line.strip())
                # print(line.strip().split()[-6:-3])
                id_atom = int(line.strip().split()[1])
                [x, y, z] = line.strip().split()[-6:-3]
                list_atoms[id_atom] = {
                    "element" : line.strip().split()[2],
                    "residue" : line.strip().split()[3],
                    "num_residue" : int(line.strip().split()[5]),
                    "chain" : line.strip().split()[4],
                    "coordinates" : (float(x),float(y),float(z))
                }
    return list_atoms


# ======================
# CALCULATE SPHERE
# ======================

def saff_kuijlaars_points(n, id_atom, center=(0.0,0.0,0.0), radius=0.0):
    """
    Generates a n-points quasi-uniform sphere based on Saff and Kuijlaars algorithm.

    Parameters
    ---
    n (int): Number of points to generate.
    center (np.array) : Coordinates of the center of the sphere. Default value : (0.0,0.0,0.0).
    radius (float) : Radius of the atom. Default value : 0.0.

    Returns
    ---
    points (ndarray): Un tableau (n, 3) avec les coordonnées x, y, z des points.
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

        points.append(id_atom, (x, y, z))

    return points

# ======================
# MAIN
# ======================

# ======================
# MAIN
# ======================

file_atoms = read_pdb(FILENAME_PROTEIN)

# print()