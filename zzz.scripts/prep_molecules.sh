#!/bin/sh

# Author: Guilherme Duarte R. Matos
# Date: March 07, 2022


############# File check ###############
helpFunction(){
    echo ""
    echo "Usage: $0 -l ligand.mol2"
    echo -e "\t-l ligand's filename"
    exit 1 # Exits script after printing help
}

# Assign typed in arguments to variables
while getopts "l:i:" opt
do
    case $opt in
        l ) LIG_NAME="$OPTARG";;
        i ) LIG_CODE="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFunction in case there's only one parameter or
# parameters are missing


if [ -z "${LIG_NAME}" ] || [ -z "${LIG_CODE}" ]
then
    echo "Some parameters are missing";
    helpFunction
fi

########################################

# Make LIG_CODE all uppercase and take only first three letters
tmp=${LIG_CODE::3}
CODE=$(echo $tmp | tr '[:lower:]' '[:upper:]')

# Create string variable without ".mol2" in the end
#name=${LIG_NAME::-4} # works in bash 4.x 
name=${LIG_NAME%.*} # works in bash 3.x (Mac included)

# Create and enter working directory
root=$(pwd)
mkdir -p ${root}/${name}
cd ${root}/${name}

## Prep for AMBER simulations
# You can check the Amber manual for more details
# use antechamber to assign force field parameters for the ligand
echo "Starting Antechamber"
# run antechamber assuming ligand is neutral ("-nc 0" flag)
antechamber -i ${root}/${name}.mol2 -fi mol2 -o ${name}_antechamber.mol2 -fo mol2 -at gaff2 -c bcc -rn ${CODE} -nc 0
# 'gaff2' is the force field
# 'bcc' indicates the method which generated the partial charges (AM1-BCC)

# check missing force field parameters with parmchk2
echo "Starting parmchk2"
parmchk2 -i ${name}_antechamber.mol2 -f mol2 -o ${name}_antechamber.frcmod

cat <<EOF > tl.in
source leaprc.gaff2
unk = loadmol2 ${name}_antechamber.mol2
loadamberparams ${name}_antechamber.frcmod
check unk
saveamberparm unk molecule.prmtop molecule.inpcrd
quit
EOF

# Run tleap
tleap -f tl.in

# Run amber_to_gmx.py
python ../amber_to_gmx.py -p molecule.prmtop -c molecule.inpcrd

# Go back to root
cd ${root}









