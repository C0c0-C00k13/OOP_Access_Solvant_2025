"""List of functions specific to an Atom object"""

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
    # logger.debug(msg=f"Mapping of the list of atom caracteristic ")
    characteristics_dict = {d['atom_name']: d for d in atom_caracteristics}
    
    # Usage of filter and map to generate the merged list
    # logger.debug(msg=f"Concatenation of dictionaries with the same value 'atom_name'.")
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

    # logger.debug(msg=f"Running 'generate_spheres_from_atom_list' function.")
    sphere_points_template = utils.generate_sphere_template(n)
    sphere_from_atom_list = {}

    for atom in atom_list:
        # logger.debug(msg=f" Atom/sphere {atom}.")
        sphere_points = []
        for point in sphere_points_template:

            x = float(atom['x']) + float(atom['radius']) * point[0]
            y = float(atom['y']) + float(atom['radius']) * point[1]
            z = float(atom['z']) + float(atom['radius']) * point[2]
            # logger.debug(msg=f" Atom/sphere {atom["atom_serial"]} : sphere points {point} - Current point [{x,y,z}].")
            sphere_points.append( (x,y,z) )

        sphere_from_atom_list[atom["atom_serial"]] = sphere_points
    # logger.debug(msg=f"'generate_spheres_from_atom_list' function : done.")
    return sphere_from_atom_list

def calculating_atom_max_exposed_surface(atom_list):
    """Returns the sphere surface area of each atom of the atom_list. 

    Parameter
    ---
    radius ([Dict]) : List of atoms.

    Returns
    ---
    ([Dict]) : List of atoms updated with the surface area of each atom.
    """

    # logger.debug(msg=f"Running 'calculating_atom_max_asa'")
    for atom in atom_list:
        logger.debug(msg=f"Current atom : {atom}")
        max_asa = utils.calculate_sphere_surface(float(atom['radius']))
        # logger.debug(msg=f"Updating current atom with {max_asa}")
        atom.update({'max_asa' : max_asa})

    return atom_list


def nb_of_points_exposed_per_atom(spheres_list, atom_list, probe):
    """
    """

    list_nb_points_exposed_per_atom = {}
    logger.debug(msg=f"1 : {spheres_list}.")
    for point_list_couple in spheres_list.items():
        key = point_list_couple[0]
        list_points = point_list_couple[1]
        logger.debug(msg=f"2 : {key, list_points}.")
        nb_points_exposed = 0
        logger.debug(msg=f"Current atom serial : {key}.")

        for point in list_points:
            if is_point_exposed(px=point[0], py= point[1], pz= point[2], atoms=atom_list, probe=probe, current_atom_serial = key):
                nb_points_exposed += 1
        list_nb_points_exposed_per_atom[key] = nb_points_exposed

    return list_nb_points_exposed_per_atom


def is_point_exposed(px, py, pz, atoms, current_atom_serial, probe):
    """ Returns a False if there is an overlap between to points of different atoms.
    
    Parameters
    ---
    px (float) : Point coordinate on the x-axis.
    py (float) : Point coordinate on the y-axis.
    pz (float) : Point coordinate on the z-axis.
    atoms (List) : List of of atoms.
    current_atom_serial (int): Serial ID of atom.  
    probe (float) : Radius of probe.
    vdw_radii (Dict) : List of radius of every atom.

    Returns
    ---
    True/ False (boolean)    
    """

    for atom in atoms:
        if atom['atom_serial'] == current_atom_serial:
            continue

        logger.debug(msg=f"Current 'other' atom : {atom}")
        atom_coordinates = (float(atom["x"]), float(atom["y"]), float(atom["z"]))
        threshold = float(atom["radius"]) + probe

        threshold_square =  threshold **2
        distance = utils.distance_between_two_3d_coordinates((px,py,pz), atom_coordinates)
        logger.debug(msg=f"Radius of current 'other' atom : {threshold} ; Square {threshold**2} ; Distance {distance}")

        if distance < threshold_square:
            # logger.debug(msg=f"False for atom {current_atom_serial} compared with atom  {atom}")
            return False
    logger.debug(msg=f"True for atom {current_atom_serial} point {px,py,pz}")
    return True

if __name__ == "__main__":
    pass