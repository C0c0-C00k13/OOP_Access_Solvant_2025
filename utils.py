"""Every static or specific function not included in an existing python library."""

import math
import logging
logger = logging.getLogger(__name__)

# Utils

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

    logger.debug(msg=f"Points of sphere : {points}.")
    return points


# def distance_between_two_3d_coordinates():

if __name__ == "main":
    pass