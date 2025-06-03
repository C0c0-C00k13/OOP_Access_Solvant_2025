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
    with open(filename, 'r') as radii_file :
        for line in radii_file:
            integrate_line = radii_file.readline()
            
            if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
                # New residu
                if integrate_line.startswith("RESIDUE"):
                    residue = integrate_line.strip().split()[2]
                    residues[residue] = []
                elif integrate_line.startswith("ATOM"):
                    residues[residue].append({
                        'name' : integrate_line.strip().split()[1], 
                        'radius' : integrate_line.strip().split()[2],
                        'polar' : integrate_line.strip().split()[3]
                    })
    return residues         

if __name__ == "__main__":
    FILENAME = "./Data/vdw.radii"
    residues = get_radii(FILENAME)
    print(residues)