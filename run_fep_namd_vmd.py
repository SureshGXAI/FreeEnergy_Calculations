#!/usr/bin/env python

import os, sys
from prody import *
from rdkit import Chem
from rdkit.Chem import AllChem
from io import StringIO
import pypdb
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import GetFormalCharge


if not os.path.exists("fep.tcl"):
    with open("fep.tcl", 'w') as w:
        w.write("""##############################################################
# FEP SCRIPT
# Jerome Henin <jhenin@ifr88.cnrs-mrs.fr>
#
# Changes:
# 2010-04-24: added runFEPmin
# 2009-11-17: changed for NAMD 2.7 keywords
# 2008-06-25: added TI routines
# 2007-11-01: fixed runFEP to handle backwards transformations
#             (i.e. dLambda < 0)
##############################################################

##############################################################
# Example NAMD input:
#
# source fep.tcl
#
# alch                  on
# alchFile              system.fep
# alchCol               B
# alchOutFreq           10
# alchOutFile           system.fepout
# alchEquilSteps        500
#
# set nSteps      5000
# set init {0 0.05 0.1}
# set end {0.9 0.95 1.0}
#
# runFEPlist $init $nSteps
# runFEP 0.1 0.9 0.1 $nSteps
# runFEPlist $end $nSteps
##############################################################

##############################################################
# proc runFEPlist { lambdaList nSteps }
#
# Run n FEP windows joining (n + 1) lambda-points
##############################################################

proc runFEPlist { lambdaList nSteps } {
    # Keep track of window number
    global win
    if {![info exists win]} {
      set win 1
    }

    set l1 [lindex $lambdaList 0]
    foreach l2 [lrange $lambdaList 1 end] {
      print [format "Running FEP window %3s: Lambda1 %-6s Lambda2 %-6s \[dLambda %-6s\]"\
        $win $l1 $l2 [expr $l2 - $l1]]
      firsttimestep 0
      alchLambda       $l1
      alchLambda2      $l2
      run $nSteps

      set l1 $l2
      incr win
    }
}


##############################################################
# proc runFEP { start stop dLambda nSteps }
#
# FEP windows of width dLambda between values start and stop
##############################################################

proc runFEP { start stop dLambda nSteps } {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    runFEPlist $ll $nSteps
}


##############################################################
##############################################################

proc runFEPmin { start stop dLambda nSteps nMinSteps temp} {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    if { $nMinSteps > 0 } { 
      alchLambda       $start
      alchLambda2      $start
      minimize $nMinSteps
      reinitvels $temp
    }

    runFEPlist $ll $nSteps
}

##############################################################
##############################################################

proc runTIlist { lambdaList nSteps } {
    # Keep track of window number
    global win
    if {![info exists win]} {
	    set win 1
    }

    foreach l $lambdaList {
	    print [format "Running TI window %3s: Lambda %-6s "	$win $l ]
	    firsttimestep 0
	    alchLambda       $l
	    run $nSteps
	    incr win
    }
}


##############################################################
##############################################################

proc runTI { start stop dLambda nSteps } {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    runTIlist $ll $nSteps
}

##############################################################
# Increment lambda and try to correct truncation errors around
# 0 and 1
##############################################################

proc increment { lambda dLambda } {
    set epsilon 1e-15
    set new [expr { $lambda + $dLambda }]

    if { [expr $new > - $epsilon && $new < $epsilon] } {
      return 0.0
    }
    if { [expr ($new - 1) > - $epsilon && ($new - 1) < $epsilon] } {
      return 1.0
    }
    return $new
}

""")


def em_namd(inpfile):

    VAR = inpfile
    em = f'''#Example NAMD input for Minimization
#Temperature & pressure control using Langevin
#set inp lindex 
#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile        {VAR}.prmtop
ambercoor       {VAR}.inpcrd

seed 35462
temperature     300

######################################################
# NONBOND CUTOFF
######################################################
switching       off
cutoff          998.0 #This_is_for_VDW
switchdist      996.0
pairlistdist    1000.0 #need at least 1.5 between cutoff and pairlist
margin          1 #default=0, once equilibrated margin can be reduced

exclude         scaled1-4
1-4scaling      0.833333
scnb            2.0
outputname      {VAR}-min
restartfreq     250     ;# 500steps = every 1ps
dcdfreq         250
xstFreq         250
binaryoutput    yes

minimize        1000
'''

    if not os.path.exists('em.namd'):
        with open('em.namd', 'w') as w:
            w.write(em)
    cmd = "namd3 +auto-provision em.namd >em-min.log"
    os.system(cmd)

    return

def equil_namd(inpfile, ligname):
    
    VAR = inpfile
    xdim, ydim, zdim, center1, center2, center3 = box_dimensions(inpfile, ligname)


    equil = f'''#############################################################
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
parmfile       {VAR}.prmtop ;#input prmtop
ambercoor      {VAR}.rst7 ;#input inpcrd

set outputname  {VAR}-equil     ;#output name

bincoordinates {VAR}-min.restart.coor
binvelocities  {VAR}-min.restart.vel
extendedSystem {VAR}-min.restart.xsc

#temperature $temperature

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
langevinTemp $temperature
langevinHydrogen off ;# don't couple langevin bath to hydrogens

# Periodic Boundary Conditions
cellBasisVector1    {xdim}  0.   0.
cellBasisVector2     0.    {ydim} 0.
cellBasisVector3     0.      0.  {zdim}
cellOrigin  {center1} {center2} {center3} 
#cellOrigin     0.   0.   0. 

#wrapAll on
wrapWater       on
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
fixedAtomsFile	{VAR}-const-min.pdb
fixedAtomsCol	B

# Output
outputName $outputname

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
reinitvels $temperature

run 50000 ;# 
'''
    if not os.path.exists('equil-amber.namd'):
        with open('equil-amber.namd', 'w') as w:
            w.write(equil)

    cmd = "namd3 +auto-provision equil-amber.namd >equil-amber.log"
    os.system(cmd)

    return 


def prod_simulations(inpfile, ligname, STIME):
    VAR = inpfile
    xdim, ydim, zdim, center1, center2, center3 = box_dimensions(f'{inpfile}-equil', ligname)

    prod_namd = f'''#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       {VAR}.prmtop ;#input prmtop
ambercoor      {VAR}.rst7 ;#input inpcrd

set outputname  {VAR}-prod

bincoordinates {VAR}-equil.restart.coor
binvelocities  {VAR}-equil.restart.vel
extendedSystem {VAR}-equil.restart.xsc

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
#wrapAll         on
wrapWater       on

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
cellBasisVector1    {xdim}  0.   0.
cellBasisVector2     0.    {ydim} 0.
cellBasisVector3     0.      0.  {zdim}
cellOrigin  {center1} {center2} {center3} 
#cellOrigin  0  0  0 
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
outputname     $outputname
restartname     $outputname.restart
binaryrestart   on

run {STIME}

    '''

    if not os.path.exists('prod-amber.namd'):
        with open('prod-amber.namd', 'w') as w:
            w.write(prod_namd)

    cmd = "namd3 +auto-provision prod-amber.namd >prod-amber.log"
    os.system(cmd)

    return

def write_fep_constraint(inpfile, ligname):
    VAR = inpfile
    fep_write = f'''#Box dimensions and constraining heavy atoms for equilibration
mol load parm7 {VAR}.prmtop namdbin {VAR}-prod.coor
set all [atomselect top all]
set lig [atomselect top "resname {ligname}"]

$all set beta 0.0
$lig set beta -1.0

$all writepdb {VAR}-solvate.fep
quit
    '''

    if not os.path.exists('fep_write.tcl'):
        with open('fep_write.tcl', 'w') as w:
            w.write(fep_write)

    cmd = "vmd -dispdev text -e fep_write.tcl"
    os.system(cmd)

    return

def fep_forward(inpfile, ligname):
    VAR = inpfile
    xdim, ydim, zdim, center1, center2, center3 = box_dimensions(f'{inpfile}-prod', ligname)

    forward_write = f'''#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       {VAR}.prmtop ;#input prmtop
ambercoor      {VAR}.rst7 ;#input inpcrd

set outputname  {VAR}-prod

bincoordinates {VAR}-prod.restart.coor
binvelocities  {VAR}-prod.restart.vel
extendedSystem {VAR}-prod.restart.xsc

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
#wrapAll         on
wrapWater       on

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
cellBasisVector1    {xdim}  0.   0.
cellBasisVector2     0.    {ydim} 0.
cellBasisVector3     0.      0.  {zdim}
cellOrigin  {center1} {center2} {center3} 
#cellOrigin  0  0  0 
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
alchFile                {VAR}-solvate.fep
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

while {{$Lambda0 < 0.99}} {{
 alchLambda $Lambda0
 set Lambda0 [expr $Lambda0 + $dLambda]
 alchLambda2 $Lambda0
 run  $nSteps
}}'''

    if not os.path.exists('fep-forward.namd'):
        with open('fep-forward.namd', 'w') as w:
            w.write(forward_write)

    cmd = "namd3 +auto-provision fep-forward.namd >fep-forward.log"
    os.system(cmd)
    return

def fep_backward(inpfile, ligname):
    VAR = inpfile
    xdim, ydim, zdim, center1, center2, center3 = box_dimensions(f'{inpfile}-prod', ligname)

    backward_write = f'''#Example NAMD input for NPT simulation for a protein
#Temperature & pressure control using Langevin
set PME 0

#######################################################
#    PDB, PSF, PAR, VEL & INPUT PARAMETERS
#######################################################
amber           on
parmfile       {VAR}.prmtop ;#input prmtop
ambercoor      {VAR}.rst7 ;#input inpcrd

set outputname  {VAR}-prod

bincoordinates {VAR}-prod.restart.coor
binvelocities  {VAR}-prod.restart.vel
extendedSystem {VAR}-prod.restart.xsc

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
#wrapAll         on
wrapWater       on

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
cellBasisVector1    {xdim}  0.   0.
cellBasisVector2     0.    {ydim} 0.
cellBasisVector3     0.      0.  {zdim}
cellOrigin  {center1} {center2} {center3} 
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
alchFile                {VAR}-solvate.fep
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

while {{$Lambda0 > 0.01}} {{
 alchLambda $Lambda0
 set Lambda0 [expr $Lambda0 - $dLambda]
 alchLambda2 $Lambda0
 run  $nSteps
}}
    '''

    if not os.path.exists('fep-backward.namd'):
        with open('fep-backward.namd', 'w') as w:
            w.write(backward_write)

    cmd = "namd3 +auto-provision fep-backward.namd >fep-backward.log"
    os.system(cmd)

    return

def parsefep(forward, backward):
    parse = f'''package require parsefep
parsefep -forward {forward} -entropy -gc 3 -gauss -backward {backward} -bar
    '''
    if not os.path.exists('parsefep-analysis.tcl'):
        with open('parsefep-analysis.tcl', 'w') as w:
            w.write(parse)
    cmd = "vmd -dispdev text -eofexit <parsefep-analysis.tcl >final-fep-out.txt"
    os.system(cmd)
    return 

def box_dimensions(inpfile, ligname):

    VAR = inpfile 
    dimension = f'''#Box dimensions and constraining heavy atoms for equilibration
mol load parm7 {VAR}.prmtop namdbin {VAR}-min.coor
set all [atomselect top all]
set prot [atomselect top "noh and (protein or nucleic or resname {ligname})"]

$all set beta 0.0
$prot set beta 1.0

$all writepdb {VAR}-const-min.pdb

#finding box dimensions 
set out [open "box-dimensions.txt" w]
set all [atomselect top "all"]
set minmax [measure minmax $all]
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
set center [measure center $all]
puts $out "xdim [lindex $vec 0]
ydim [lindex $vec 1]
zdim [lindex $vec 2]
center1 [lindex $center 0]
center2 [lindex $center 1]
center3 [lindex $center 2]"
close $out 
exit
    '''

    if not os.path.exists('box_dimension.tcl'):
        with open('box_dimension.tcl', 'w') as w:
            w.write(dimension)
    cmd = "vmd -dispdev text -eofexit <box_dimension.tcl"
    os.system(cmd)

    with open("box-dimensions.txt", 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'xdim' in line:
                xdim = str(line.split()[1].strip())
            elif 'ydim' in line:
                ydim = str(line.split()[1].strip())
            elif 'zdim' in line:
                zdim = str(line.split()[1].strip())
            elif 'center1' in line:
                center1 = str(line.split()[1].strip())
            elif 'center2' in line:
                center2 = str(line.split()[1].strip())
            elif 'center3' in line:
                center3 = str(line.split()[1].strip())
            else:
                pass

    return xdim, ydim, zdim, center1, center2, center3

def main(inputname, ligand_name, STIME):

    em_namd(inputname)
    equil_namd(inputname, ligand_name)
    prod_simulations(inputname, ligand_name, STIME)
    write_fep_constraint(inputname, ligand_name)
    fep_forward(inputname, ligand_name)
    fep_backward(inputname, ligand_name)
    forward = 'forward-alchemy.fepout'
    backward = 'backward-alchemy.fepout'
    parsefep(forward, backward)
    return 

if __name__ == "__main__":
    if len(sys.argv) > 3:
        inputname = sys.argv[1]
        ligand_name = sys.argv[2]
        stime = sys.argv[3]
        main(inputname, ligand_name, stime)
    else:
        print("Usauage: python3 run_fep_namd_vmd.py inputname ligand_name time(ps)")
