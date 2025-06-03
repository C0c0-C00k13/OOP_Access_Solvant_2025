"""Retrieves atoms radii from file 'vdw.radii'"""
import time

FILENAME = "./Data/vdw.radii"
output = ""

with open(FILENAME, 'r') as radii_file :
    for line in radii_file:
        integrate_line = radii_file.readline()
        if not integrate_line.startswith('#') and not integrate_line.startswith("\n"):
            # print(integrate_line)
            output = output+integrate_line

print("OUTPUT")
time.sleep(2)
print (output)