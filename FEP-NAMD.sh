#!/bin/bash

VAR=8fza_chainA_h
STIME=50000

source ~/anaconda3/etc/profile.d/conda.sh
conda activate AmberTools23


cat >>vmd.tcl <<EOF
#run vmd to separate the chains, molecules
set pdb [lindex \$argv 0]
set resid [lindex \$argv 1]
mol load pdb \$pdb.pdb
set rna [atomselect top "nucleic"]
set water [atomselect top "water"]
set ions [atomselect top "ions"]
set lig [atomselect top "not water and not ions and not nucleic"]
\$rna writepdb \$pdb-rna.pdb
\$lig writepdb \$pdb-ligand.pdb
\$water writepdb \$pdb-water.pdb
\$ions writepdb \$pdb-ions.pdb
exit
EOF
vmd -dispdev text -eofexit <vmd.tcl -args $VAR

CHG=-1
sed -i "1d" $VAR-ligand.pdb
LIG=$(awk '{print $4}' $VAR-ligand.pdb | head -n1)

echo $CHG $LIG 
antechamber -i $VAR-ligand.pdb -fi pdb -nc $CHG -o $VAR-ligand.mol2 -fo mol2 -c bcc -s 1 -pf yes -at gaff -m 1 
parmchk2 -i $VAR-ligand.mol2 -f mol2 -o $VAR-ligand.frcmod

cat $VAR-rna.pdb $VAR-ligand.pdb >$VAR-comp.pdb
sed -i "s/END/TER/" $VAR-comp.pdb

cat >>test.sh <<EOF
source leaprc.RNA.OL3
source leaprc.DNA.bsc1
source leaprc.gaff
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
loadamberparams $VAR-ligand.frcmod
$LIG = loadmol2 $VAR-ligand.mol2
list
check $LIG
saveoff $LIG $VAR-ligand.lib
saveamberparm $LIG $VAR-ligand.prmtop $VAR-ligand.rst7

#loadamberparams $VAR-ligand.frcmod
loadoff $VAR-ligand.lib
comp = loadpdb $VAR-comp.pdb
saveamberparm comp $VAR-complex.prmtop $VAR-complex.rst7
savepdb comp $VAR-complex.pdb

source leaprc.water.tip3p
solvatebox comp TIP3PBOX 12.0
addionsrand comp Na+ 0
saveamberparm comp $VAR-solvate.prmtop $VAR-solvate.rst7
savepdb comp $VAR-solvate.pdb
quit
EOF
tleap -f test.sh
rm test.sh vmd.tcl leap.log sqm*

cat >>namd-min.conf <<EOF
#Example NAMD input for Minimization
#Temperature & pressure control using Langevin
#set inp lindex 
#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber 		on
parmfile	$VAR-solvate.prmtop	
ambercoor 	$VAR-solvate.rst7 

seed 35462
temperature	300

######################################################
# NONBOND CUTOFF
######################################################
switching       off
cutoff		998.0 #This_is_for_VDW
switchdist	996.0
pairlistdist	1000.0 #need at least 1.5 between cutoff and pairlist
margin		1 #default=0, once equilibrated margin can be reduced

exclude         scaled1-4
1-4scaling      0.833333
scnb 		2.0
outputname 	$VAR-min
restartfreq     250     ;# 500steps = every 1ps
dcdfreq         250
xstFreq         250
binaryoutput	yes

minimize 	1000
EOF

namd3 +p16 namd-min.conf >namd-min.log
##########################################
cat >>box-dimensions.tcl <<EOF
#Box dimensions and constraining heavy atoms for equilibration
mol load parm7 $VAR-solvate.prmtop namdbin $VAR-min.coor
set all [atomselect top all]
set prot [atomselect top "noh and (nucleic or resname $LIG)"]

\$all set beta 0.0
\$prot set beta 1.0

\$all writepdb $VAR-const-min.pdb

#finding box dimensions 
set out [open "box-dimensions.txt" w]
set all [atomselect top "all"]
set minmax [measure minmax \$all]
set vec [vecsub [lindex \$minmax 1] [lindex \$minmax 0]]
set center [measure center \$all]
puts \$out "xdim [lindex \$vec 0]
ydim [lindex \$vec 1]
zdim [lindex \$vec 2]
center1 [lindex \$center 0]
center2 [lindex \$center 1]
center3 [lindex \$center 2]"

exit
EOF
vmd -dispdev text -eofexit <box-dimensions.tcl
xdim=`grep xdim box-dimensions.txt | sed 's/xdim/ /g'`
ydim=`grep ydim box-dimensions.txt | sed 's/ydim/ /g'`
zdim=`grep zdim box-dimensions.txt | sed 's/zdim/ /g'`
center1=`grep center1 box-dimensions.txt | sed 's/center1/ /g'`
center2=`grep center2 box-dimensions.txt | sed 's/center2/ /g'`
center3=`grep center3 box-dimensions.txt | sed 's/center3/ /g'`
echo $xdim
echo $ydim
echo $zdim
echo $center1
echo $center2
echo $center3
##########################################

cat >>namd-equil-amber.conf <<EOF
#############################################################
## JOB DESCRIPTION ##
#############################################################

#############################################################
## ADJUSTABLE PARAMETERS ##
#############################################################
set temperature 310

firsttimestep 0

#############################################################
## SIMULATION PARAMETERS ##
#############################################################

seed	145872
amber           on
parmfile       $VAR-solvate.prmtop ;#input prmtop
ambercoor      $VAR-solvate.rst7 ;#input inpcrd

set outputname  $VAR-equil     ;#output name

bincoordinates $VAR-min.restart.coor
binvelocities  $VAR-min.restart.vel
extendedSystem $VAR-min.restart.xsc

#temperature \$temperature

# Force-Field Parameters
exclude scaled1-4
1-4scaling 0.833333333
cutoff 12.
switching off
switchdist 10.
pairlistdist 13.5

# Integrator Parameters
timestep 2.0 ;# 2fs/step
rigidBonds all ;# needed for 2fs steps
nonbondedFreq 1
fullElectFrequency 2
stepspercycle 10

# Constant Temperature Control
langevin on ;# do langevin dynamics
langevinDamping 5 ;# damping coefficient (gamma) of 5/ps
langevinTemp \$temperature
langevinHydrogen off ;# don't couple langevin bath to hydrogens

# Periodic Boundary Conditions
cellBasisVector1    $xdim  0.   0.
cellBasisVector2     0.    $ydim 0.
cellBasisVector3     0.      0.  $zdim
cellOrigin     0.   0.   0. 
#cellOrigin     $center1 $center2 $center3

wrapAll on

# PME (for full-system periodic electrostatics)
PME yes
PMETolerance 1.0e-6
PMEInterpOrder 4
PMEGridSpacing      1.0
PMEPencils          1

useGroupPressure yes ;# needed for rigidBonds
useFlexibleCell no
useConstantArea no

#Note the lack of langevinpiston, which is the pressure control
fixedAtoms	on
fixedAtomsFile	$VAR-const-min.pdb
fixedAtomsCol	B

# Output
outputName \$outputname

restartfreq 500 ;# 500steps = every 1ps
dcdfreq 2500
xstFreq 2500
outputEnergies 1000
outputPressure 1000

#############################################################
## EXECUTION SCRIPT ##
#############################################################

# Minimization
minimize 100
reinitvels \$temperature

run 50000 ;# 
EOF
namd3 +auto-provision namd-equil-amber.conf>namd-equil-amber.log
###########################################

cat >>namd-prod-amber.conf <<EOF 
#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       $VAR-solvate.prmtop ;#input prmtop
ambercoor      $VAR-solvate.rst7 ;#input inpcrd

set outputname  $VAR-prod

bincoordinates $VAR-equil.restart.coor
binvelocities  $VAR-equil.restart.vel
extendedSystem $VAR-equil.restart.xsc

seed 35462
#temperature	310
timestep	2.0 #timestep_2_fs
stepspercycle	20 #number_of_steps_before_making_new_pairlist


######################################################
# OUTPUT PARAMETERS
######################################################
outputEnergies  500 #display every 1ps
outputTiming    1000 #display every 2ps
xstFreq         5000 #save every 10ps
dcdFreq         1000 #save every 2ps
restartfreq     2500 #update restart file every 5ps
wrapAll         on

rigidBonds all #SHAKE_H-bonds
rigidTolerance 0.00001

MTSAlgorithm impulse   # Integrator
fullElectFrequency 1 
nonbondedFreq  1
longSplitting c1

######################################################
#Temperature & Pressure Control
######################################################
Langevin		on
LangevinDamping		1
LangevinTemp		310
LangevinHydrogen	on
LangevinPiston 		on
LangevinPistonTarget 	1.01325 #target pressure in (bar)
LangevinPistonPeriod 	200 #oscilation period in (fs)
LangevinPistonDecay 	100 #dumping time scale (fs)
LangevinPistonTemp 	310 # barostat noise T = target T

######################################################
# NONBOND CUTOFF
######################################################
switching       on
cutoff		12.0 #This_is_for_VDW
switchdist	10.0
pairlistdist	14.0 #need at least 1.5 between cutoff and pairlist
margin		1 #default=0, once equilibrated margin can be reduced

exclude         scaled1-4
1-4scaling      0.833333333
dielectric      1.0

######################################################
# CELL PARAMETERS
######################################################

useFlexibleCell	yes # Allow cell dimension to change along the Z-axis
useConstantArea yes
cellBasisVector1    $xdim  0.   0.
cellBasisVector2     0.    $ydim 0.
cellBasisVector3     0.      0.  $zdim
cellOrigin  0  0  0 
######################################################
# PMEwald
######################################################
PME                 yes
PMETolerance 	    1.0e-6
PMEInterpOrder 	    4
PMEGridSpacing      1.0
PMEPencils          1

######################################################
# OUTPUT FILES
######################################################

binaryoutput    on
outputname     \$outputname
restartname     \$outputname.restart
binaryrestart   on

run $STIME
EOF
namd3 +auto-provision namd-prod-amber.conf>namd-prod-amber.log

########################################################
cat >>fep-write-vmd.tcl <<EOF
#Box dimensions and constraining heavy atoms for equilibration
mol load parm7 $VAR-solvate.prmtop namdbin $VAR-prod.coor
set all [atomselect top all]
set lig [atomselect top "resname $LIG"]

\$all set beta 0.0
\$lig set beta -1.0

\$all writepdb $VAR-solvate.fep
EOF
vmd -dispdev text -eofexit <fep-write-vmd.tcl

#########################################################

cat >>namd-fepforward-amber.conf <<EOF 
#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       $VAR-solvate.prmtop ;#input prmtop
ambercoor      $VAR-solvate.rst7 ;#input inpcrd

set outputname  $VAR-prod

bincoordinates $VAR-prod.restart.coor
binvelocities  $VAR-prod.restart.vel
extendedSystem $VAR-prod.restart.xsc

seed 35462
#temperature	310
timestep	2.0 #timestep_2_fs
stepspercycle	20 #number_of_steps_before_making_new_pairlist


######################################################
# OUTPUT PARAMETERS
######################################################
outputEnergies  500 #display every 1ps
outputTiming    1000 #display every 2ps
xstFreq         5000 #save every 10ps
dcdFreq         1000 #save every 2ps
restartfreq     2500 #update restart file every 5ps
wrapAll         on

rigidBonds all #SHAKE_H-bonds
rigidTolerance 0.00001

MTSAlgorithm impulse   # Integrator
fullElectFrequency 1 
nonbondedFreq  1
longSplitting c1

######################################################
#Temperature & Pressure Control
######################################################
Langevin		on
LangevinDamping		1
LangevinTemp		310
LangevinHydrogen	on
LangevinPiston 		on
LangevinPistonTarget 	1.01325 #target pressure in (bar)
LangevinPistonPeriod 	200 #oscilation period in (fs)
LangevinPistonDecay 	100 #dumping time scale (fs)
LangevinPistonTemp 	310 # barostat noise T = target T

######################################################
# NONBOND CUTOFF
######################################################
switching       on
cutoff		12.0 #This_is_for_VDW
switchdist	10.0
pairlistdist	14.0 #need at least 1.5 between cutoff and pairlist
margin		1 #default=0, once equilibrated margin can be reduced

exclude         scaled1-4
1-4scaling      0.833333333
dielectric      1.0

######################################################
# CELL PARAMETERS
######################################################

useFlexibleCell	yes # Allow cell dimension to change along the Z-axis
useConstantArea yes
cellBasisVector1    $xdim  0.   0.
cellBasisVector2     0.    $ydim 0.
cellBasisVector3     0.      0.  $zdim
cellOrigin  0  0  0 
######################################################
# PMEwald
######################################################
PME                 yes
PMETolerance 	    1.0e-6
PMEInterpOrder 	    4
PMEGridSpacing      1.0
PMEPencils          1

######################################################
# OUTPUT FILES
######################################################
# OUTPUT AND RESTART

outputname              alchemy
#restartname             alchemy

binaryoutput            no
binaryrestart           yes

# FEP

source                  fep.tcl

alch                    on
alchType                fep
alchFile                $VAR-solvate.fep
alchCol                  B
alchOutFile             forward-alchemy.fepout
alchOutFreq              10

alchVdwLambdaEnd        1.0
alchElecLambdaStart     0.5
alchVdwShiftCoeff       5.0
alchDecouple            no

# LOOP OVER LAMBDA-STATES

alchEquilSteps           250000
set nSteps               750000

set Lambda0    0.0
set dLambda    0.1

while {\$Lambda0 < 0.99} {
 alchLambda \$Lambda0
 set Lambda0 [expr \$Lambda0 + \$dLambda]
 alchLambda2 \$Lambda0
 run  \$nSteps
}
EOF
namd3 +auto-provision namd-fepforward-amber.conf>namd-fepforward-amber.log


#########################################################

cat >>namd-fepback-amber.conf <<EOF 
#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       $VAR-solvate.prmtop ;#input prmtop
ambercoor      $VAR-solvate.rst7 ;#input inpcrd

set outputname  $VAR-prod

bincoordinates $VAR-prod.restart.coor
binvelocities  $VAR-prod.restart.vel
extendedSystem $VAR-prod.restart.xsc

seed 35462
#temperature	310
timestep	2.0 #timestep_2_fs
stepspercycle	20 #number_of_steps_before_making_new_pairlist


######################################################
# OUTPUT PARAMETERS
######################################################
outputEnergies  500 #display every 1ps
outputTiming    1000 #display every 2ps
xstFreq         5000 #save every 10ps
dcdFreq         1000 #save every 2ps
restartfreq     2500 #update restart file every 5ps
wrapAll         on

rigidBonds all #SHAKE_H-bonds
rigidTolerance 0.00001

MTSAlgorithm impulse   # Integrator
fullElectFrequency 1 
nonbondedFreq  1
longSplitting c1

######################################################
#Temperature & Pressure Control
######################################################
Langevin		on
LangevinDamping		1
LangevinTemp		310
LangevinHydrogen	on
LangevinPiston 		on
LangevinPistonTarget 	1.01325 #target pressure in (bar)
LangevinPistonPeriod 	200 #oscilation period in (fs)
LangevinPistonDecay 	100 #dumping time scale (fs)
LangevinPistonTemp 	310 # barostat noise T = target T

######################################################
# NONBOND CUTOFF
######################################################
switching       on
cutoff		12.0 #This_is_for_VDW
switchdist	10.0
pairlistdist	14.0 #need at least 1.5 between cutoff and pairlist
margin		1 #default=0, once equilibrated margin can be reduced

exclude         scaled1-4
1-4scaling      0.833333333
dielectric      1.0

######################################################
# CELL PARAMETERS
######################################################

useFlexibleCell	yes # Allow cell dimension to change along the Z-axis
useConstantArea yes
cellBasisVector1    $xdim  0.   0.
cellBasisVector2     0.    $ydim 0.
cellBasisVector3     0.      0.  $zdim
cellOrigin  0  0  0 
######################################################
# PMEwald
######################################################
PME                 yes
PMETolerance 	    1.0e-6
PMEInterpOrder 	    4
PMEGridSpacing      1.0
PMEPencils          1

######################################################
# OUTPUT FILES
######################################################
# OUTPUT AND RESTART

outputname              alchemy
#restartname             alchemy

binaryoutput            no
binaryrestart           yes

# FEP

source                  fep.tcl

alch                    on
alchType                fep
alchFile                $VAR-solvate.fep
alchCol                  B
alchOutFile             backward-alchemy.fepout
alchOutFreq              10

alchVdwLambdaEnd        1.0
alchElecLambdaStart     0.5
alchVdwShiftCoeff       5.0
alchDecouple            no

# LOOP OVER LAMBDA-STATES

alchEquilSteps           250000
set nSteps               750000

set Lambda0    1.0
set dLambda    0.1

while {\$Lambda0 > 0.01} {
 alchLambda \$Lambda0
 set Lambda0 [expr \$Lambda0 - \$dLambda]
 alchLambda2 \$Lambda0
 run  \$nSteps
}
EOF
namd3 +auto-provision namd-fepback-amber.conf>namd-fepback-amber.log

######################################################
#ParseFEP analysis
######################################################
cat >> parsefep.tcl <<EOF 
package require parsefep
parsefep -forward forward-alchemy.fepout -entropy -gc 3 -gauss -backward backward-alchemy.fepout -bar
EOF

vmd -dispdev text -eofexit <parsefep.tcl >final-fep-out.txt 
