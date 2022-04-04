import parmed
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-m","--mol2",dest="mol2file", help="mol2 file with atomic charges", metavar="FILE")
(options, args) = parser.parse_args()
mol2_filename = options.mol2file

#Correct the mol2 file partial atom charges to have a total net integer molecule charge                  
mol2f = parmed.formats.Mol2File
mol2f.write(parmed.load_file(mol2_filename).fix_charges(), mol2_filename)

