"""Define Atom object"""

import logging
logger = logging.getLogger(__name__)


# class Atom():
# Attributes
# Methods
# determine_the_atom_exposed_surface


def generating_atom_with_caracteristic_list(simple_atom_list, atom_characteristics):
    logger.info(msg=f"Generate atom with caracteristics list\n{simple_atom_list}")
    complete_atom_list = []
    for atom in simple_atom_list:
        logger.debug(msg=f"Looking for the right atom caracteristic for : {atom}")
        atom_characteristic = filter(lambda x: x['atom_name'] == atom['atom_name'], atom_characteristics)
        complete_atom = atom | atom_characteristic
        logger.debug(msg=f"Complete atom: {complete_atom}")
        complete_atom_list.append(complete_atom)
    return complete_atom_list

if __name__ == "__main__":
    pass