import parmed
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p","--parm",dest="parmfile", help="Amber topology file",metavar="FILE")
parser.add_option("-c","--cord",dest="coordfile", help="Amber Coordinate File",metavar="FILE")

(options, args) = parser.parse_args()

gromacs_topology = parmed.load_file(options.parmfile, options.coordfile) 

gromacs_topology.save(options.parmfile.split('.')[0]+'.top', format='gromacs')
gromacs_topology.save(options.parmfile.split('.')[0]+'.gro')

print("Conversion from Amber to GMX done!")
