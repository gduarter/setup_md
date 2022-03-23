# Import general packages
import os
import argparse

# Import modules for dealing with chemical information
from rdkit import Chem


### Auxiliary functions

def mol2_mol_supplier(file):
    ''' This function extracts all the molecules in the multi-molecule
        MOL2 file `file` and returns a list of rdkit.Chem.rdchem.mol 
        object.
        
        Variable         I/O          dtype           default value?
        ------------------------------------------------------------
        file              I           string                  None
        mols              O           list                    N/A
        mols[i]           O           rdkit.Chem.rdchem.mol   N/A
        
    '''
    
    # with open statements make sure you do not forget to close the
    # file you are using.
    with open(file, 'r') as f:
        
        # f.readline() reads each line separately
        line = f.readline()
        
        # Loop goes line by line to find situations where conditions
        # are met
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                
                # All the data for a single molecule will be stored 
                # in mol.
                mol = []
                mol.append(line)
                line = f.readline()
                
                # "@<TRIPOS>MOLECULE" marks the beginning of a new 
                # molecule in the MOL2 file. If the code finds it,
                # it skips the next block of lines within the block
                # of the while loop.
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                
                # Makes final adjustments to the data. It must look
                # like the MOL2 file of a single molecule.
                mol[-1] = mol[-1].rstrip() # removes blank line at file end
                block = ",".join(mol).replace(',','')
                
                # Converts the data of a single molecule to a 
                # rdkit.Chem.rdchem.mol object.
                m=Chem.MolFromMol2Block(block,
                                        sanitize=False,
                                        cleanupSubstructures=False)
            
    return m 


###### Main

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--mol2", dest="mol2file", help="input mol2 file", nargs=1, required=True)
args = parser.parse_args()
if len(args.mol2file) != 1:
    parser.print_help()
    sys.exit()

mol = mol2_mol_supplier( args.mol2file[0] )
#print(mol.GetNumAtoms())
mol = Chem.rdmolops.RemoveHs( mol )
smi = Chem.MolToSmiles( mol, canonical=True )
print( smi )
















