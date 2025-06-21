"""Every static or specific function not included in an existing python library."""

import math
import logging
logger = logging.getLogger(__name__)

def generate_sphere_template(n):
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

    logger.debug(msg=f"Generates a {n} points quasi-uniform sphere based on Saff and Kuijlaars algorithm.")
    points = []
    offset = 2.0 / n
    increment = math.pi * (3.0 - math.sqrt(5))

    for k in range(n):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * increment
        x =  math.cos(phi) * r
        z =  math.sin(phi) * r
        points.append((x, y, z))

    # logger.debug(msg=f"Points of sphere : {points}.")
    return points


def distance_between_two_3d_coordinates(coord_a, coord_b)->float:
    """Calculates distance between 2 points
    
    Parameters
    ---
    coord_a : Coordinates from either a point of a an atom or an atom.
    coord_b : Coordinates from either a point of a an atom or an atom.
    
    Returns
    ---
    distance (float) : Distance from 2 points (a, b).
    """

    # logger.debug(msg=f"Running 'distance_between_two_3d_coordinates'...")

    dx = coord_a[0] - coord_b[0]
    dy = coord_a[1] - coord_b[1]
    dz = coord_a[2] - coord_b[2]
    return dx*dx + dy*dy + dz*dz


def calculate_sphere_surface(radius)->float:
    """Returns the sphere surface area of a sphere. 

    Parameter
    ---
    radius (float) : Radius of a sphere.

    Returns
    ---
    (float) : Surface area of the sphere.
    """

    # logger.debug(msg=f"Running 'calculate_sphere_surface'...")
    return 4 * math.pi * radius**2


def calculate_point_surface(sphere_surface:float, point_per_sphere:int)->float:
    """Returns the fraction a sphere surface occupied by 1 point.

    Parameters
    ---
    sphere_surface (float) :
    point_per_sphere (int) :

    Returns
    ---
    (float) : Surface occupied by 1 point of a sphere.
    
    """
    return sphere_surface / point_per_sphere

if __name__ == "main":
    pass