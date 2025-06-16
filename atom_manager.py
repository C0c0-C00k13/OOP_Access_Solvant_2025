"""Define Atom object"""

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
    simple_atom_list ([Dict]) :
    atom_caracteristics ([Dict]) :
    
    Returns
    ---
    complete_atom_list (List) :
    """ 
    # Lookup dictionary for atom_caracteristics by element
    logger.debug(msg=f"Mapping of the list of atom caracteristic ")
    characteristics_dict = {d['element']: d for d in atom_caracteristics}
    
    # Usage of filter and map to generate the merged list
    logger.debug(msg=f"Concatenation of dictionaries with the same value 'atom_name'.")
    complete_atom_list = list(
        filter(
            None,
            map(
                lambda atom: {**atom, **characteristics_dict[atom['atom_name']]} 
                if atom['atom_name'] in characteristics_dict else None,
                simple_atom_list
            )
        )
    )

    # logger.info(msg=f"Generate atom with caracteristics list\n{simple_atom_list}")
    # complete_atom_list = []

    # for atom in simple_atom_list:
    #     logger.debug(msg=f"Looking for the right atom caracteristic for : {atom}")
    #     logger.debug(msg=f"First item of atom_caracteristic list : {atom_caracteristics[0]}")
    #     print(atom_caracteristics[0])
        
    #     atom_caracteristic = filter(lambda x: x['atom_name'] == atom['atom_name'], atom_caracteristics)
    #     complete_atom = atom | atom_caracteristic
    #     logger.debug(msg=f"Complete atom: {complete_atom}")
    #     complete_atom_list.append(complete_atom)
    return complete_atom_list

if __name__ == "__main__":
    pass