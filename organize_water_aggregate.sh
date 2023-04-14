#!/usr/bin/env bash

# Author: Guilherme Duarte Ramos Matos
# Date: March 2022

###### Help Function ########
helpFunction(){
    echo -e "\tUsage: $0 -u solutes_csv -t temperature -p pressure -n number_of_solute_mols"
    exit 1
}

# Assign typed arguments to variables
while getopts "u:s:t:p:n:" opt
do
    case $opt in
        u ) SOLUTE="$OPTARG";;
        t ) temperature="$OPTARG";;
        p ) pressure="$OPTARG";;
        n ) NSOLU="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${SOLUTE}" ] || [ -z "${temperature}" ] || [ -z "${pressure}" ] || [ -z "${NSOLU}" ] 
then
    echo "You are misusing this script"
    helpFunction
fi

# define paths
root=$(pwd)
scriptdir=${root}/zzz.scripts


array=(${SOLUTE})
for val in "${array[@]}"
do

    # Create PDB files for each molecule
    echo "Creating PDB files from SMILES"
    pdbdir=${root}/001.molecules

    mkdir -p ${pdbdir}
    cd ${pdbdir}
    python3 ${scriptdir}/create_pdbfiles.py -c ${root}/$val

    #exit 0
    # Remove lines containing 'CONECT'
    for elem in *.pdb
    do
        sed '/CONECT/d' $elem > tmp.pdb
        mv tmp.pdb $elem
    done
    #exit 0

    echo "Parameterizing each generated molecule"
    for elem in *.pdb
    do

        echo ""
        echo ${elem}
        echo "All molecules are assumed to be neutral"
        molname=${elem%.*}
        # Define 3-letter code
        tmp_code=${molname::3}
        code=$(echo $tmp_code | tr '[:lower:]' '[:upper:]')
        echo "Starting Antechamber with ${molname}, code ${code}"
        # Fixing residue names
        antechamber -i ${molname}.pdb -fi pdb -o ${molname}_renamed.pdb -fo pdb -rn ${code}

        #exit 0
        rm ${molname}.pdb
        mv ${molname}_renamed.pdb ${molname}.pdb
        # Creating simulation files
        antechamber -i ${molname}.pdb -fi pdb -o ${molname}.mol2 -fo mol2 -at gaff2 -c bcc -rn ${code} -nc 0
        # Fix excess charges
        python3 ${scriptdir}/fix_charges.py -m ${molname}.mol2
        echo "Starting parmchk2"
        parmchk2 -i ${molname}.mol2 -f mol2 -o ${molname}.frcmod
        rm sqm.*

        # Transform Cl and Br in mol2 files in CL and BR.
        # tleap will fail if you don't do it
        echo "Checking if there are halogens with inappropriate names"
        python3 ${scriptdir}/fixChlorineBromine.py -f ${molname}.mol2
        if [ -f "${molname}_fixed.mol2" ]
        then
            mv ${molname}.mol2 ${molname}.mol2.bckp
            mv ${molname}_fixed.mol2 ${molname}.mol2
        fi

    done

    cd ${root}

done

# Create another simulation directory
boxes=${root}/002.aggregates_in_water
mkdir -p ${boxes}
cd ${boxes}

# For each template, generate solvated files with monomers
echo "Generate solvated boxes"
for mol1 in ${root}/001.molecules/*.pdb
do
        # Create name of molecule 1
        tmp_mol1=${mol1%.pdb}
        name1=${tmp_mol1##*/}

        # Create directory and copy files
        mkdir -p ${boxes}/${name1}_in_water
        cd ${boxes}/${name1}_in_water
        echo ${mol1}
        cp ${mol1%.*}.* ${boxes}/${name1}_in_water

        # Prep packmol files
        echo "Prep packmol input files"
        cat <<EOF > input.inp
tolerance 3.0 # tolerance distance
output ${name1}_no_water.pdb # output file name
filetype pdb # output file type
#
# Create a box of ${name1} molecules
#
structure ${name1}.pdb
number ${NSOLU} # Number of molecules
resnumbers 3 # Sequential numbering
#fixed 30. 30. 30. 0. 0. 0.
inside cube 0. 0. 0. 40. # x, y, z coordinates of box, and length of box in Angstroms
add_amber_ter
end structure
#
EOF
        # Run packmol
        echo "Run packmol for ${name1}"
        packmol < input.inp

        #Define 3-letter code for each molecule
        tmp_code1=${name1::3}
        code1=$(echo $tmp_code1 | tr '[:lower:]' '[:upper:]')

        # Create tleap input
        echo "Create tleap input"
        cat <<EOF > tl.in
source leaprc.gaff2
source leaprc.water.tip3p
${code1} = loadmol2 ${name1}.mol2
loadamberparams ${name1}.frcmod
fullbox = loadPdB ${name1}_no_water.pdb
solvatebox fullbox TIP3PBOX 20.0
setbox fullbox centers 0.0

saveAmberParm fullbox ${name1}_in_water.prmtop ${name1}_in_water.inpcrd
quit
EOF
        echo "Running tleap"
        tleap -f tl.in

        # Create gromacs gro e top files
        python3 ${scriptdir}/amber_to_gmx.py -p ${name1}_in_water.prmtop -c ${name1}_in_water.inpcrd

        # Create minimization MDP file
        cat<<EOF > 01.minimization.mdp
; RUN CONTROL PARAMETERS
integrator               = steep
; Start time and timestep in ps
nsteps                   = 50000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; ENERGY MINIMIZATION OPTIONS
; Force tolerance and initial step-size
emtol                    = 100
emstep                   = 0.01
; Max number of iterations in relax_shells
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints
fcstep                   = 0
; Frequency of steepest descents steps when doing CG
nstcgsteep               = 1000
nbfgscorr                = 10

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 0
; Output frequency for energies to log file and energy file
nstlog                   = 0
nstenergy                = 0
; Output frequency and precision for xtc file
nstxout-compressed       = 0
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
cutoff-scheme            = Verlet
; nblist update frequency
nstlist                  = 10
; ns algorithm (simple or grid)
ns_type                  = grid
; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; nblist cut-off
rlist                    = 1.2

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = pme
rcoulomb-switch          = 0
rcoulomb                 = 1.2
; Dielectric constant (DC) for cut-off or DC of reaction field
; Method for doing Van der Waals
vdw-type                 = Cut-off
; cut-off lengths
rvdw                     = 1.0
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = AllEnerPres
; Extension of the potential lookup tables beyond the cut-off
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.10
; FFT grid size, when a value is 0 fourierspacing will be used
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = yes

; OPTIONS FOR BONDS
constraints              = ;hbonds ;all-bonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
unconstrained-start      = no
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 1e-04
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no
EOF

        # Create NVT equilibration MDP file
        cat<<EOF > 02.equil_nvt.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 0
; Output frequency for energies to log file and energy file
nstlog                   = 0
nstenergy                = 0
; Output frequency and precision for xtc file
nstxout-compressed       = 0
compressed-x-precision   = 1000

; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz

; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16

; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = AllEnerPres

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
; Groups to couple separately
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
; REMEMBER THIS SIMULATION USES A LANGEVIN INTEGRATOR
tau_t                    = 2.0
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = no

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = yes
gen_temp                 = ${temperature}
gen_seed                 = 2022
ld_seed                  = -1

; OPTIONS FOR BONDS
constraints              = all-bonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = no
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 1e-04
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no
EOF

        cat<<EOF > 03.equil_npt.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 0
; Output frequency for energies to log file and energy file
nstlog                   = 0
nstenergy                = 0
; Output frequency and precision for xtc file
nstxout-compressed       = 0
compressed-x-precision   = 1000

; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz

; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16

; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = AllEnerPres

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 2.0
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = all-bonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = yes
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 1e-04
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = no
ld_seed                  = -1
EOF
        cat<<EOF > 04.equil_npt2.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 0
; Output frequency for energies to log file and energy file
nstlog                   = 0
nstenergy                = 0
; Output frequency and precision for xtc file
nstxout-compressed       = 0
compressed-x-precision   = 1000

; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz

; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16

; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = AllEnerPres

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 2.0
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 5
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = all-bonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = yes
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 1e-04
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = no
ld_seed                  = -1
EOF

        cat<<EOF > 05.prod.mdp
;gromacs stuff
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 15000000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 100

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 1000
nstvout                  = 1000
nstfout                  = 1000
; Checkpointing helps you continue after crashes
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxout-compressed       = 1000
compressed-x-precision   = 1000

; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16

; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = AllEnerPres

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 2.0
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = Parrinello-Rahman
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 5
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = all-bonds
; Type of constraint algorithm
constraint-algorithm     = Lincs
; Do not constrain the start configuration
continuation             = yes
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 1e-04
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = no
ld_seed                  = -1
EOF

        cat<<EOF > gromacs.md.slurm
#!/bin/bash
#SBATCH --partition=int_medium
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=${name1}

## Prepare o ambiente
module purge
module load gnu8/8.3.0
module add openmpi3/3.1.4
## GROMACS para ser rodado em CPUs
module load gromacs-2020.3-gcc-8.3.0-djklvgf

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR
echo -e "Nome do job:   \$SLURM_JOB_NAME"
echo -e "Nos alocados:  \$SLURM_JOB_NODELIST"

# Mensagem
echo "Molecular dynamics simulation to evaluate ligand stability in the binding"
echo "in the binding site by calculating the RMSD from the production trajectory."

## Execute o programa
## Minimizacao da energia (optimizacao da geometria)
echo "Starting minimization"
gmx grompp -f 01.minimization.mdp -p ${name1}_in_water.top -c ${name1}_in_water.gro -o 01.minimization.tpr
gmx mdrun -ntomp \${SLURM_CPUS_PER_TASK} -deffnm 01.minimization
echo "Minimization ended"

## Equilibracao NVT
echo "Starting NVT equilibration"
gmx grompp -f 02.equil_nvt.mdp -p ${name1}_in_water.top -c 01.minimization.gro -o 02.equil_nvt.tpr
gmx mdrun -ntomp \${SLURM_CPUS_PER_TASK} -deffnm 02.equil_nvt
echo "NVT equilibration ended"

## Equilibracao NPT com barostato de Berendsen
echo "Starting NPT equilibration (Berendsen)"
gmx grompp -f 03.equil_npt.mdp -p ${name1}_in_water.top -c 02.equil_nvt.gro -o 03.equil_npt.tpr
gmx mdrun -ntomp \${SLURM_CPUS_PER_TASK} -deffnm 03.equil_npt
echo "NPT equilibration (Berendsen) ended"

## Equilibracao NPT com barostato de Parrinello-Rahman
echo "Starting NPT equilibration (Parrinello-Rahman)"
gmx grompp -f 04.equil_npt2.mdp -p ${name1}_in_water.top -c 03.equil_npt.gro -o 04.equil_npt2.tpr
gmx mdrun -ntomp \${SLURM_CPUS_PER_TASK} -deffnm 04.equil_npt2
echo "NPT equilibration (Parrinello-Rahman) ended"

## Producao
echo "Starting production"
gmx grompp -f 05.prod.mdp -p ${name1}_in_water.top -c 04.equil_npt2.gro -o 05.prod.tpr
gmx mdrun -ntomp \${SLURM_CPUS_PER_TASK} -deffnm 05.prod
echo "Production ended"

echo "Simulation done!"
EOF

done

cd ${root}
