"""Functions to load atoms from any files."""
import logging
logger = logging.getLogger(__name__)

def load_data_from_file(filename, header):
    """Returns a list of atoms """
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
    data_from_line = line.strip().split()
    if len(data_from_line) != len(header):
        logger.warning(f"Invalid structure. {line}")
        return None
    data = dict(zip(header, data_from_line))
    logger.debug(msg=f"{data}")
    return data



if __name__ == "__main__":
    pass