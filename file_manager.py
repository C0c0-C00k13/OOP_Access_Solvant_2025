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
    # logger.debug(msg=f"Data list before return : {data_list}")
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
    # logger.debug(msg=f"Data : {data}")
    return data

def write_residue_asa_output(filename, res_asa, chain_stats):
    """Write the output file of ASA of residue.
    
    Parameters
    ---
    filename (str) : Name of the output file.
    res_asa (Dict) : ASA of residues.
    chain_stats (Dict) : ASA of chains.
    
    Returns
    ---
    None
    
    """

    max_asa = {
    "ALA": 113, "ARG": 241, "ASN": 158, "ASP": 151, "CYS": 140,
    "GLN": 189, "GLU": 183, "GLY": 85, "HIS": 194, "ILE": 182,
    "LEU": 180, "LYS": 211, "MET": 204, "PHE": 218, "PRO": 143,
    "SER": 122, "THR": 146, "TRP": 259, "TYR": 229, "VAL": 160
    }
    with open(filename, 'w') as f:
        f.write("\nResidue ASA and RSA:\n")
        f.write(f"{'Chain':<5} {'ResID':<6} {'ResName':<7} "
                f"{'TotalASA':<10} {'RSA(%)':<8} "
                f"{'MainASA':<10} {'MainRSA(%)':<8} "
                f"{'SideASA':<10} {'SideRSA(%)':<8} "
                f"{'PolarASA':<10} {'PolarRSA':<10} "
                f"{'ApolarASA':<11} {'ApolarRSA':<10}\n")
        for (chain, res_id, res_name), data in sorted(res_asa.items()):
            max_ref = max_asa.get(res_name, 200)
            total = data["total"]
            main_chain = data["main"]
            side_chain = data["side"]
            polar = data["polar"]
            apolar = data["apolar"]
            rsa_total = (total / max_ref) * 100
            rsa_main = (main_chain / max_ref) * 100
            rsa_side = (side_chain / max_ref) * 100
            rsa_polar = (polar / max_ref) * 100
            rsa_apolar = (apolar / max_ref) * 100
            f.write(f"{chain:<5} {res_id:<6} {res_name:<7} "
                    f"{total:<10.2f} {rsa_total:<8.2f} "
                    f"{main_chain:<10.2f} {rsa_main:<10.2f} "
                    f"{side_chain:<10.2f} {rsa_side:<10.2f} "
                    f"{polar:<10.2f} {rsa_polar:<10.2f} "
                    f"{apolar:<11.2f} {rsa_apolar:<10.2f}\n")
        f.write("\nPer-Chain ASA Summary:\n")
        f.write(f"{'Chain':<5} {'Main ASA':<12} {'Side ASA':<12} {'Polar ASA':<12} {'Apolar ASA':<12} {'Total ASA':<12}\n")
        for chain, stats in sorted(chain_stats.items()):
            f.write(f"{chain:<5} {stats['main']:<12.2f} {stats['side']:<12.2f} "
                    f"{stats['polar']:<12.2f} {stats['apolar']:<12.2f} {stats['total']:<12.2f}\n")


def write_log(filename, date, pdb_file, points, probe, radii_file, rsa_filename, asa_filename):
    """Write the log file of ASA of residue.

    Parameters
    ---
    filename (str) : Name of the log file.
    date (str) : Date of Execution.
    rsa_filename (str) : Name of the ASA per residue file.
    asa_filename (str) : Name of the ASA per atom file.

    Returns
    ---
    None
    
    """
    with open(filename, 'w') as f_o:
        f_o.write(f"{'DATE':<30}: {date}\n")
        f_o.write("\n" + "*" * 40 + "\n")

        f_o.write("INPUT\n")
        f_o.write(f"{'PDB FILE':<30}: {pdb_file}\n")
        f_o.write("*" * 40 + "\n" * 2)

        f_o.write(f"{'INCLUDE HETATM':<30}: NO\n")
        f_o.write(f"{'POINTS PER SPHERE':<30}: {points}\n")
        f_o.write(f"{'PROBE RADIUS':<30}: {probe}\n")
        f_o.write(f"{'Custom radii file':<30}: {radii_file}\n")

        f_o.write("\n" + "*" * 40 + "\n")
        f_o.write("OUTPUT\n")
        f_o.write("*" * 40 + "\n" * 2)
        f_o.write(f"{'ACCESSIBLE SURFACE RESIDUE':<30}: '{rsa_filename}'\n")
        f_o.write(f"{'ACCESSIBLE SURFACE ATOM':<30}: '{asa_filename}'\n")


def write_atom_asa_output(filename, atom_asa_list):
    """Write the output file of ASA of atom.

    Parameters
    ---
    filename (str) : Name of the output file.
    atom_asa_list ([Dist]) : ASA of atoms.

    Returns
    ---
    None
    """
    with open(filename, 'w') as f:
        for atom in atom_asa_list:
            f.write(f"ATOM\t{atom['atom_serial']}\t{atom['atom_name']}\t"
                    f"{atom['res_name']}\t{atom['chain']}\t{atom['res_id']}\t"
                    f"{atom['x']:.3f}\t{atom['y']:.3f}\t{atom['z']:.3f}\t"
                    f"{atom['asa']:.3f}\t{atom['radius']:.3f}\n")





if __name__ == "__main__":
    pass