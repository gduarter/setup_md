import sys
import argparse
from optparse import OptionParser 


parser = OptionParser()
parser.add_option("-f", "--file", help="input MOL2 file", dest="mol2file", metavar="FILE")
(options, args) = parser.parse_args()

# Open MOL2 file
mol2filename = options.mol2file
mol2file = open(mol2filename, "r+")

# Make appropriate changes
lines = mol2file.readlines()
newlines = []
for line in lines:
    if "Cl" in line:
        newline = line.replace("Cl", "CL")
        newlines.append(newline)
    elif "Br" in line:
        newline = line.replace("Br", "BR")
        newlines.append(newline)
    else:
        newlines.append(line)

# Create new MOL2 file
newmol2filename = mol2filename.split(".")[0] + "_fixed.mol2"
lines = "".join(newlines)
with open(newmol2filename, "w+") as f:
    f.write(lines)

# Close original MOL2 file
mol2file.close()

