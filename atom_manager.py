"""List of functions specific to an Atom object"""

from collections import defaultdict
import logging
import utils
logger = logging.getLogger(__name__)


def generating_atom_with_caracteristic_list(simple_atom_list, atom_caracteristics):
    """Returns a list

    Parameters
    ---
    simple_atom_list ([Dict]) : List of atoms.
    atom_caracteristics ([Dict]) : List of supplementary caracteristics.
    
    Returns
    ---
    complete_atom_list (List[Dict]) :
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
        # logger.debug(msg=f"Current atom : {atom}")
        max_asa = utils.calculate_sphere_surface(float(atom['radius']))
        # logger.debug(msg=f"Updating current atom with {max_asa}")
        atom.update({'max_asa' : max_asa})

    return atom_list


def nb_of_points_exposed_per_atom(spheres_list, atom_list, probe):
    """
    """

    list_nb_points_exposed_per_atom = {}
    # logger.debug(msg=f"1 : {spheres_list}.")
    for point_list_couple in spheres_list.items():
        key = point_list_couple[0]
        list_points = point_list_couple[1]
        # logger.debug(msg=f"2 : {key, list_points}.")
        nb_points_exposed = 0
        # logger.debug(msg=f"Current atom serial : {key}.")

        for point in list_points:
            if is_point_exposed(point, atom_list, key, probe):
                nb_points_exposed += 1
        # logger.debug(msg=f"Atom {key} : currently {nb_points_exposed} exposed point(s).")
        list_nb_points_exposed_per_atom[key] = nb_points_exposed

    return list_nb_points_exposed_per_atom


def is_point_exposed(current_point, atoms, current_atom_serial, probe):
    """ Returns a False if there is an overlap between to points of different atoms.
    
    Parameters
    ---
    current_point () : Point .
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

        # logger.debug(msg=f"Current 'other' atom : {atom}")
        atom_coordinates = (float(atom["x"]), float(atom["y"]), float(atom["z"]))
        threshold = float(atom["radius"]) + probe
        threshold_square =  threshold **2

        distance = utils.distance_between_two_3d_coordinates(current_point, atom_coordinates)
        # logger.debug(msg=f"Radius of current 'other' atom : {threshold} ; Square {threshold**2} ; Distance {distance}")

        if distance < threshold_square:
            # logger.debug(msg=f"False for atom {current_atom_serial} compared with atom  {atom}")
            return False
    # logger.debug(msg=f"Atom {current_atom_serial} - is_exposed : True -point exposed : {current_point}.")
    return True


def calculate_atom_asa(atoms_list, exposed_points, nb_points_per_sphere):
    """Returns the ASA (Accessible Solvent Area) of each atom and each residue.

    Parameters
    ---
    atoms_list ([Dict]) : List of atoms to process to calculate ASA.
    exposed_points (Dict) : Dict of exposed points keyed by atom serial number.
    nb_points_per_sphere (int) : Number of points used for the sphere surface.

    Returns
    ---
    atom_asa_list (Dict) : ASA per atom (angstrom^2 and percent).
    residue_asa (Dict)   : Total ASA per residue.
    """

    main_chain_elements = {"N", "CA", "C", "O", "OXT", "H", "HA"}
    polar_list = {"N", "O", "S"}
    chain_stats = defaultdict(lambda: {
        "main": 0.0, "side": 0.0,
        "polar": 0.0, "apolar": 0.0,
        "total": 0.0
    })
    atom_asa_list = {}
    residue_asa = {}

    for atom in atoms_list:
        logger.debug(msg=f"ASA Calculus - Atom about to about to be processed : {atom}")
        atom_max_asa = atom["max_asa"]
        point_area = utils.calculate_point_surface(atom_max_asa, nb_points_per_sphere)
        atom_key = atom["atom_serial"]
        atom_angstrom_asa = exposed_points[atom_key] * point_area
        atom_percent_asa = (atom_angstrom_asa / atom_max_asa) * 100

        # Add to atom-level dictionary
        atom_asa_list[atom_key] = {
            'angstrom_asa': atom_angstrom_asa,
            'percent_asa': atom_percent_asa
        }
        logger.debug(msg=f"Integration in Atom Dict of - Atom n°{atom_key}, {atom["atom_name"]} | Current value of de 'atom_asa_list' {atom_asa_list}.")

        # Build residue key
        residue_key = (atom["chain"], atom["num_residue"], atom["residue"])
        chain = atom["chain"]

        if residue_key not in residue_asa:
            residue_asa[residue_key] = {
                "total": 0.0,
                "polar": 0.0, "apolar": 0.0,
                "main": 0.0, "side": 0.0
            }
        residue_asa[residue_key]["total"] += atom_angstrom_asa

        chain_stats[chain]["total"] += atom_angstrom_asa

        # Discriminates MAIN chain an SIDE chain
        if atom['atom_name'] in main_chain_elements:
            chain_stats[chain]["main"] += atom_angstrom_asa
            residue_asa[residue_key]["main"] += atom_angstrom_asa
        else:
            chain_stats[chain]["side"] += atom_angstrom_asa
            residue_asa[residue_key]["side"] += atom_angstrom_asa

        # Discriminates POLAR elements ande NON-POLAR elements
        if atom['element'] in polar_list:
            residue_asa[residue_key]["polar"] += atom_angstrom_asa
            chain_stats[chain]["polar"] += atom_angstrom_asa
        else:
            residue_asa[residue_key]["apolar"] += atom_angstrom_asa
            chain_stats[chain]["apolar"] += atom_angstrom_asa

        logger.debug(f"Atom {atom_key} ({atom['atom_name']}): ASA = {atom_angstrom_asa:.2f} Å², ASA (%) = {atom_percent_asa:.2f} %")
        logger.debug(f"Updated residue ASA for {residue_key}: MAIN : {residue_asa[residue_key]["main"]:.2f}\n\
SIDE : {residue_asa[residue_key]["side"]:.2f} \nPOLAR : {residue_asa[residue_key]["polar"]:.2f} \nAPOLAR : {residue_asa[residue_key]["apolar"]:.2f} Å²")

    return atom_asa_list, residue_asa, chain_stats

# def print_result(atom_list, atom_asa_list):
#     for atom, atom_asa in zip(atom_list, atom_asa_list):
#         record = ["record"]
#         atom_serial = atom["atom_serial"]
#         atom_name = atom["atom_name"]
#         residue = atom["residue"]
#         chain = atom["chain"]
#         num_residue = atom["num_residue"]
#         x, y, z = atom["x"], atom["y"], atom["z"]
#         radius = atom["radius"]

if __name__ == "__main__":
    pass