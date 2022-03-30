#!/usr/bin/env bash

# Author: Guilherme Duarte Ramos Matos
# Date: March 2022

###### Help Function ########
helpFunction(){
    echo -e "\tUsage: $0 -u solutes_csv -s solvent_csv"
    exit 1
}

# Assign typed arguments to variables
while getopts "u:s:" opt
do
    case $opt in
        u ) SOLUTE="$OPTARG";;
        s ) SOLVENT="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${SOLUTE}" ] || [ -z "${SOLVENT}" ]
then
    echo "You are misusing this script"
    helpFunction
fi

# define paths
root=$(pwd)
scriptdir=${root}/zzz.scripts
pdbdir=${root}/001.single_mol_files
packdir=${root}/004.packmol_files
amberdir=${root}/005.amber_files

# Define relevant numbers for packmol
NSOLU=1
NSOLV=150

# Define relevant numbers for gromacs
temperature=298.15 #kelvin
pressure=1.01325 #bar


array=(${SOLUTE} ${SOLVENT})
for val in "${array[@]}"
do

    # Create PDB files for each molecule
    echo "Creating PDB files from SMILES"
    if [ $val == $SOLUTE ]
        then
        pdbdir=${root}/001.solutes
    else
        pdbdir=${root}/002.solvents
    fi

    mkdir -p ${pdbdir}
    cd ${pdbdir}
    python3 ${scriptdir}/create_pdbfiles.py -c ${root}/$val
    # Remove lines containing 'CONECT'
    for elem in *.pdb
    do
        sed '/CONECT/d' $elem > tmp.pdb
        mv tmp.pdb $elem
    done

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
        rm ${molname}.pdb
        mv ${molname}_renamed.pdb ${molname}.pdb
        # Creating simulation files
        antechamber -i ${molname}.pdb -fi pdb -o ${molname}.mol2 -fo mol2 -at gaff2 -c bcc -rn ${code} -nc 0
        echo "Starting parmchk2"
        parmchk2 -i ${molname}.mol2 -f mol2 -o ${molname}.frcmod
        rm sqm.*

        # Transform Cl and Br in mol2 files in CL and BR.
        # tleap will fail if don't do it
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
boxes=${root}/003.initial_boxes
mkdir -p ${boxes}
cd ${boxes}

# For each template, generate solvated files with monomers
echo "Generate solvated boxes"
for mol1 in ${root}/001.solutes/*.pdb
do
    for mol2 in ${root}/002.solvents/*.pdb
    do
        # Create name of molecule 1
        tmp_mol1=${mol1%.pdb}
        name1=${tmp_mol1##*/}
        # Create name of molecule 2
        tmp_mol2=${mol2%.pdb}
        name2=${tmp_mol2##*/}

        # Create directory and copy files
        mkdir -p ${boxes}/${name1}_in_${name2}
        cd ${boxes}/${name1}_in_${name2}
        echo ${mol1}
        echo ${mol2}
        cp ${mol1%.*}.* ${boxes}/${name1}_in_${name2}
        cp ${mol2%.*}.* ${boxes}/${name1}_in_${name2}

        # Prep packmol files
        echo "Prep packmol input files"
        cat <<EOF > input.inp
tolerance 1.5 # tolerance distance
output ${name1}_in_${name2}.pdb # output file name
filetype pdb # output file type
#
# Create a box of ${name1} in ${name2} molecules
#
structure ${name1}.pdb
number ${NSOLU} # Number of molecules
resnumbers 3 # Sequential numbering
inside cube 0. 0. 0. 30. # x, y, z coordinates of box, and length of box in Angstroms
add_amber_ter
end structure
#
structure ${name2}.pdb
number ${NSOLV} # Number of molecules
resnumbers 3 # Sequential numbering
inside cube 0. 0. 0. 30.
add_amber_ter
end structure
EOF
        # Run packmol
        echo "Run packmol for ${name1}_in_${name2}"
        packmol < input.inp

        #Define 3-letter code for each molecule
        tmp_code1=${name1::3}
        code1=$(echo $tmp_code1 | tr '[:lower:]' '[:upper:]')
        tmp_code2=${name2::3}
        code2=$(echo $tmp_code2 | tr '[:lower:]' '[:upper:]')

        # Create tleap input
        echo "Create tleap input"
        cat <<EOF > tl.in
source leaprc.gaff2

${code1} = loadmol2 ${name1}.mol2
loadamberparams ${name1}.frcmod
${code2} = loadmol2 ${name2}.mol2
loadamberparams ${name2}.frcmod

fullbox = loadPdB ${name1}_in_${name2}.pdb

setbox fullbox centers

saveAmberParm fullbox ${name1}_in_${name2}.prmtop ${name1}_in_${name2}.inpcrd
quit
EOF
        echo "Running tleap"
        tleap -f tl.in

        # Create gromacs gro e top files
        python3 ${scriptdir}/amber_to_gmx.py -p ${name1}_in_${name2}.prmtop -c ${name1}_in_${name2}.inpcrd

        # Create minimization MDP file
        cat<<EOF > minimization.mdp
; RUN CONTROL PARAMETERS
integrator               = steep
; Start time and timestep in ps
nsteps                   = 2500
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1

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
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxtcout-compressed     = 0
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
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
constraints              = hbonds
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
        cat<<EOF > equil_nvt.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 25000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 10

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxtcout-compressed     = 0
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
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
constraints              = hbonds
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

        # Create NPT equilibration MDP file
        cat<<EOF > equil_npt.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 25000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxtcout-compressed     = 0
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
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

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 0.1
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = hbonds
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

        # Create second NPT equilibration MDP file
        cat<<EOF > equil_npt2.mdp
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 25000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxtcout-compressed     = 0
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
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

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 0.1
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = hbonds
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

        # Create production stage MDP file
        cat<<EOF > prod.mdp
;gromacs stuff
; RUN CONTROL PARAMETERS
integrator               = sd
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 2500000
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Checkpointing helps you continue after crashes
nstcheckpoint            = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
nstxtcout-compressed     = 1000
compressed-x-precision   = 1000

; NEIGHBORSEARCHING PARAMETERS
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

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling
Tcoupl                   = no
tc-grps                  = System
; Time constant (ps) and reference temperature (K)
tau_t                    = 0.1
ref_t                    = ${temperature}
; Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = isotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 1
compressibility          = 4.5e-5
ref_p                    = ${pressure}

; OPTIONS FOR BONDS
constraints              = hbonds
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
        cd ${boxes}
    done
done

cd ${root}
