import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-m", "--mol2", help="input MOL2 file", dest="mol2file", metavar="FILE")
parser.add_option("-c", "--chg", help="CHG file created by Multiwfn", dest="chgfile", metavar="FILE")

(options, args) = parser.parse_args()

chgfile = options.chgfile
charges = np.genfromtxt(chgfile, usecols=4, dtype=float)

filename = options.mol2file

line_id1 = '@<TRIPOS>ATOM'
line_id2 = '@<TRIPOS>BOND'
f = open(filename, "r+")
lines = f.readlines()
for idx, line in enumerate(lines):
    if line_id1 in line:
        start = idx
    if line_id2 in line:
        end = idx
f.close()

worklines = lines[start+1: end]
newlines = []
for i, line in enumerate(worklines):
    tmplist = line.split()[:-1]
    tmplist.append(charges[i])
    tmp = f'{tmplist[0]:>7} {tmplist[1]:<8}{tmplist[2]:>11}{tmplist[3]:>11}{tmplist[4]:>11} {tmplist[5]:<4}{tmplist[6]:>8} {tmplist[7]:>3}{round(tmplist[8],6):>15}\n'
    newlines.append(tmp)
    #print(tmp)
    
new_mol2 = lines[:start+1] + newlines + lines[end:]
#print(np.array(new_mol2))

with open(f'{filename[:-5]}_charged.mol2', 'w+') as g:
    for line in new_mol2:
        g.write(line)

print(f"QM charges added to {filename[:-5]}")

