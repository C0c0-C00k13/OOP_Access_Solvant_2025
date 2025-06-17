"""Define Atom object"""

import utils
import logging
logger = logging.getLogger(__name__)


# class Atom():
# Attributes
# Methods
# determine_the_atom_exposed_surface


def generating_atom_with_caracteristic_list(simple_atom_list, atom_caracteristics):
    """Returns a list

    Parameters
    ---
    simple_atom_list ([Dict]) : List of atoms.
    atom_caracteristics ([Dict]) : List of supplementary caracteristics.
    
    Returns
    ---
    complete_atom_list (List) :
    """ 

    # Lookup dictionary for atom_caracteristics by element
    logger.debug(msg=f"Mapping of the list of atom caracteristic ")
    characteristics_dict = {d['atom_name']: d for d in atom_caracteristics}
    
    # Usage of filter and map to generate the merged list
    logger.debug(msg=f"Concatenation of dictionaries with the same value 'atom_name'.")
    complete_atom_list = list(
        map(
            lambda atom: {**atom, **characteristics_dict[atom['atom_name']]} 
            if atom['atom_name'] in characteristics_dict else None,
            simple_atom_list
        )
    )

    return complete_atom_list

def generate_spheres_from_atom_list(atom_list, n):
    """Returns a list of spheres representing every atom of the list.
    
    Parameters
    ---
    atom () : List of atoms/ Atom represented by the sphere.s.
    n (int) : Number of points used to generate a sphere.

    Returns
    ---
    sphere_from_atom_list : List of spheres representing every atom of the list.
    """

    sphere_points_template = utils.generate_sphere_template(n)
    sphere_from_atom_list = {}

    for atom in atom_list:
        logger.debug(msg=f" Atom/sphere {atom}.")
        sphere_points = []
        for point in sphere_points_template:

            x = float(atom['x']) + float(atom['radius']) * point[0]
            y = float(atom['y']) + float(atom['radius']) * point[1]
            z = float(atom['z']) + float(atom['radius']) * point[2]
            # logger.debug(msg=f" Atom/sphere {atom["atom_serial"]} : sphere points {point} - Current point [{x,y,z}].")
            sphere_points.append( (x,y,z) )

        sphere_from_atom_list[atom["atom_serial"]] = sphere_points

    return sphere_from_atom_list


if __name__ == "__main__":
    pass