"""Protein Execution file to analyze the solvant accessible surface area of protein"""
import logging
import atom_manager
import file_manager

FILENAME_RADIUS = "./Data/vdw.radii"
FILENAME_PROTEIN = "./Data/2c8r.pdb"
RADIUS_PROBE = 1.4
ATOM_FILE_HEADER = ["record", "atom_serial", "atom_name", "residue", "chain", "num_residue", "x", "y", "z", "occupancy", "b-factor", "element"]
ATOM_CARACTERISTIC_FILE_HEADER = ["record", "atom_name", "radius", "polarity"]
SPHERE_NB_POINTS = 6
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

# functions to calculate RSA and ASA
# Protein Analysis


if __name__ == "__main__":
    logging.basicConfig(filename='./Logs/main.log', level=logging.DEBUG)
    logger.info(msg='Started')

    try:
        logger.info(msg="Loading atoms list from file...")
        atom_list = file_manager.load_data_from_file(FILENAME_PROTEIN, ATOM_FILE_HEADER)
        logger.debug(msg=f"Values contained in 'atom_list' variable : {atom_list}")
        atom_caracteristics = file_manager.load_data_from_file(FILENAME_RADIUS, ATOM_CARACTERISTIC_FILE_HEADER)
        logger.debug(msg=f"Values contained in 'atom_caracteristics' variable : {atom_caracteristics}")
        logger.info(msg="Loading atoms from list - DONE")

        logger.info(msg="Generating atoms list with complete caracteristics...")
        complete_atom_list = atom_manager.generating_atom_with_caracteristic_list(atom_list, atom_caracteristics)
        logger.debug(msg=f"Values contained in 'complete_atom_list' variable : {complete_atom_list}")
        logger.info(msg="Generating atoms list with complete caracteristics - DONE")

        logger.info(msg="Generating list of spheres represented by atoms...")
        spheres_list = atom_manager.generate_spheres_from_atom_list(complete_atom_list, SPHERE_NB_POINTS)
        logger.debug(msg=f"Values contained in 'spheres_list' variable : {spheres_list}")
        logger.info(msg="Generating list of spheres represented by atoms - DONE")

        logger.info(msg="Calculating Maximum Accessible Surface of each atom...")
        new_complete_atom_list = atom_manager.calculating_atom_max_asa(complete_atom_list)
        logger.debug(msg=f"Values contained in 'new_complete_atom_list' variable : {new_complete_atom_list}")
        logger.info(msg="Calculating Maximum Accessible Surface of each atom - DONE")

    except:
        logger.error(msg="Something went wrong during the process.")
    logger.info(msg='End.')
    # pass