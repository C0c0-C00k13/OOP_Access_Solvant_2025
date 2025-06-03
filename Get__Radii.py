"""Retrieves atoms radii from file 'vdw.radii'"""
import time

FILENAME = "./Data/vdw.radii"
output = ""
residues = {}

with open(FILENAME, 'r') as radii_file :
    for line in radii_file:
        integrate_line = radii_file.readline()
        
        if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
            # New residu
            if integrate_line.startswith("RESIDUE"):
                # print(integrate_line.strip().split()[1:])
                residue = integrate_line.strip().split()[2]
                residues[residue] = []
            elif integrate_line.startswith("ATOM"):
                residues[residue].append({
                    'name' : integrate_line.strip().split()[1], 
                    'radius' : integrate_line.strip().split()[2],
                    'polar' : integrate_line.strip().split()[3]
                })
                
                # residue['type'],residue['name'],residue['nb_atoms'] = 
            # print(integrate_line)
            # output = output+integrate_line

print(residues)

# print("OUTPUT")
# time.sleep(2)
# print (output)