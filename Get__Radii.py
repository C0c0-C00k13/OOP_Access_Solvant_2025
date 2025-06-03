"""Retrieves atoms radii from file 'vdw.radii'"""
import time


def get_radii(filename: str):
    """
    Returns the radii atoms based on each type of residue.
    ---
    Attributes
        filename:str
    ---
    Returns
        residues:Dict
    """

    output = ""
    residues = {}
    atoms = {}
    with open(filename, 'r') as radii_file :
        for line in radii_file:
            integrate_line = radii_file.readline()
            
            if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
                # New residu
    
                if integrate_line.startswith("ATOM"):
                    atom_i = integrate_line.strip().split()[1]
                    atoms[atom_i] = []
                    atoms[atom_i].append({
                        'radius' : integrate_line.strip().split()[2],
                        'polar' : integrate_line.strip().split()[3]
                    })
    return atoms

    #             if integrate_line.startswith("RESIDUE"):
    #                 residue = integrate_line.strip().split()[2]
    #                 residues[residue] = []
    #             elif integrate_line.startswith("ATOM"):
    #                 residues[residue].append({
    #                     'name' : integrate_line.strip().split()[1], 
    #                     'radius' : integrate_line.strip().split()[2],
    #                     'polar' : integrate_line.strip().split()[3]
    #                 })
    # return residues         

if __name__ == "__main__":
    # Get every radii for each atom of each residue 
    FILENAME = "./Data/vdw.radii"
    residues = get_radii(FILENAME)
    # print(residues)
    # print(residues.keys())
    print(residues['N'][0]['radius'])