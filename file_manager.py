"""Functions to load atoms from any files and write files."""

import logging
logger = logging.getLogger(__name__)


def load_data_from_file(filename, header):
    """Returns a list of data extracted from the file accessed.
     
    Parameters
    ---
    filename (str) : Name of the file accessed.
    header (List) : List of each column describing the data retrieved from the file. 

    Returns
    ---
    data_list ([Dict]) :
    """
    logger.info(msg=f"Loading data from '{filename}' ...")
    data_list = []
    with open(filename,"r") as f_in:
        for line in f_in:
            if line.startswith("ATOM"):
                data = generate_data_from_line(line, header)
                if data is not None:
                    data_list += [data]
    logger.debug(msg=f"Data list before return : {data_list}")
    return data_list

def generate_data_from_line(line, header):
    """Returns structured data generated out of the provided line.

    Parameters
    ---
    line (str) :
    header (List) :

    Returns
    ---
    data (Dict) :
    """
    data_from_line = line.strip().split()
    if len(data_from_line) != len(header):
        logger.warning(f"Invalid structure. {line}")
        return None
    data = dict(zip(header, data_from_line))
    logger.debug(msg=f"{data}")
    return data

def generate_output_rsa_file(output_filename):
    """Generates and write an output file at the format of a ASA/RSA file.
    
    Parameters
    ---
    output_filename (str) :

    Returns
    ---
    None
    """
    with open(file= output_filename, mode='w', encoding='utf-8') as output:
        pass



if __name__ == "__main__":
    pass