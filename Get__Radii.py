import time

FILENAME = "./Data/vdw.radii"
file_content=""
output = ""
with open(FILENAME, 'r') as radii_file :
    for line in radii_file:
        integrate_line = radii_file.readline()
        file_content = file_content+integrate_line
        if not integrate_line.startswith('#'):
            # print(integrate_line)
            output = output+integrate_line

print("INITIAL FILE")
time.sleep(2)
print(file_content)
time.sleep(3)
print("OUTPUT")
time.sleep(2)
print (output)