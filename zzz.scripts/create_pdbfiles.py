
import sys
import argparse
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem


# Define arguments read by this script
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csvfile", help="input CSV file",
                    type=str, nargs=1, required=True)
args = parser.parse_args()
if len(args.csvfile) != 1:
    parser.print_help()
    sys.exit()

filename = args.csvfile[0]

# Open CSV
df = pd.read_csv(filename, comment=";")
# Create mol objects
df["RDMol"] = df["SMILES"].apply(lambda smi: Chem.MolFromSmiles(smi))
df["RDMol_H"] = df["RDMol"].apply(lambda mol: Chem.AddHs(mol))
df["Fail"] = df["RDMol_H"].apply(lambda mol: AllChem.EmbedMolecule(mol))

#print(df.head())
#print(df.columns)
# Loop through dataframe and write files
for idx, row in df.iterrows():
    print(Chem.MolToPDBBlock(row["RDMol_H"]),
          file=open(f'{row["Name"]}.pdb', 'w+'))

          
