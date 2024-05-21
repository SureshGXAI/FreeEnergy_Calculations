#############################################################################
# This program calculates absolute binding free energy of a protein-drug complex 
# It uses decoupling approach in GROMACS
# It also uses PMX, PharmED, MDAnalysis libraries to prepare the systems
# Authors: Suresh Gorle
# Date: 17 May 2024
#############################################################################

import os, sys
from pmx import *
from pmx import gmx
import glob, subprocess
import MDAnalysis as mda
import warnings
warnings.filterwarnings('ignore')


em="""
; Run control
integrator               = steep 
nsteps                   = 5000
; EM criteria and other stuff
emtol                    = 100
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10
; Output control
nstlog                   = 1
nstenergy                = 1
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.2
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no
; Free energy control stuff
free_energy              = yes
init_lambda_state        = {1}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = {0}      ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = yes
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
vdw_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; We are not transforming any bonded or restrained interactions
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Not doing simulated temperting here
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10
; No velocities during EM 
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

nvt = """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.2
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 298
; Pressure coupling is off for NVT
Pcoupl                   = No
tau_p                    = 0.5
compressibility          = 4.5e-05
ref_p                    = 1.0
; Free energy control stuff
free_energy              = yes
init_lambda_state        = {1}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = {0}      ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = yes
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
vdw_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; We are not transforming any bonded or restrained interactions
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Not doing simulated temperting here
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10
; Generate velocities to start
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = -1
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Do not constrain the starting configuration
continuation             = no
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

npt = """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000    ; 100 ps
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.2
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 298 
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda_state        = {1}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = {0}      ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = yes
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
vdw_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; We are not transforming any bonded or restrained interactions
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Not doing simulated temperting here
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NVT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""

md = """; Run control
integrator               = sd       ; Langevin dynamics
tinit                    = 0
dt                       = 0.002
nsteps                   = 500000     ; 1 ns
nstcomm                  = 100
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 0
nstlog                   = 500
nstenergy                = 500
nstxout-compressed       = 0
; Neighborsearching and short-range nonbonded interactions
cutoff-scheme            = verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2
; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.2
; van der Waals
vdwtype                  = cutoff
vdw-modifier             = potential-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
; Temperature coupling
; tcoupl is implicitly handled by the sd integrator
tc_grps                  = system
tau_t                    = 1.0
ref_t                    = 298 
; Pressure coupling is on for NPT
Pcoupl                   = Parrinello-Rahman 
tau_p                    = 1.0
compressibility          = 4.5e-05
ref_p                    = 1.0 
; Free energy control stuff
free_energy              = yes
init_lambda_state        = {1}
delta_lambda             = 0
calc_lambda_neighbors    = 1        ; only immediate neighboring windows
couple-moltype           = {0}  ; name of moleculetype to decouple
couple-lambda0           = vdw      ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = yes
; Vectors of lambda specified here
; Each combination is an index that is retrieved from init_lambda_state for each simulation
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
vdw_lambdas              = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; We are not transforming any bonded or restrained interactions 
bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Not doing simulated temperting here
temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3 
nstdhdl                  = 10
; Do not generate velocities
gen_vel                  = no 
; options for bonds
constraints              = h-bonds  ; we only have C-H bonds here
; Type of constraint algorithm
constraint-algorithm     = lincs
; Constrain the starting configuration
; since we are continuing from NPT
continuation             = yes 
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 12
"""



gmx.set_gmxlib()

def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)

def clean_backup_files(path):
    toclean = glob.glob('{0}/*#'.format(path))
    for clean in toclean:
        os.remove(clean)

def split_complex(complex_id, protein_id, ligand_id):
    u = mda.Universe(str(complex_id)+".pdb")
    p = u.select_atoms("protein")
    l = u.select_atoms("not protein")

    mda.Writer(str(protein_id)+".pdb", p.n_atoms).write(p)
    mda.Writer(str(ligand_id)+".pdb", l.n_atoms).write(l)
    return


complex_id = sys.argv[1].split('.')[0]
Receptor = "protein"
Ligand = "ligand"

#Split the complex into protein and ligand
split_complex(complex_id, Receptor, Ligand)


ptop = 'topol.top'
pgro = 'protein.gro'
ltop = 'ligand.top'
lgro = 'ligand.gro'
seed = '001430'


subprocess.call('acpype -i {0}.pdb'.format(Ligand), shell=True)

# Rewriting ligand topology
import parmed as pmd
from parmed import gromacs, amber, unit as u
top = pmd.load_file('{0}.acpype/{1}_GMX.top'.format(Ligand, Ligand), xyz='{0}.acpype/{1}_GMX.gro'.format(Ligand, Ligand))
top.save(ltop)
top.save(lgro)

# get ligand name 
for line in open('{0}.pdb'.format(Ligand), 'r'):
    if line.startswith('ATOM'):
        ligname = line.split()[3]

# Protein preparation
subprocess.call('echo 12 | gmx pdb2gmx -f {0}.pdb -ter -ignh -water tip3p -o {1}'.format(Receptor, pgro), shell=True)

# System Preparation
subprocess.call('pmx abfe -pt {0} -lt {1} -pc {2} -lc {3} --build --seed {4}'.format(ptop, ltop, pgro, lgro, seed), shell=True)


LIGGRO = "ligand/genion.gro"
LIGTOP = "ligand/ligand.top"
LIG = ligname
COMPGRO = "complex/genion.gro"
COMPTOP = "complex/complex.top"

sim_dir = os.getcwd()

#Ligand only calculations
os.mkdir("LigandOnly")
os.chdir("LigandOnly")


for lamda in range(0, 21):
    os.mkdir("lambda"+str(lamda))
    os.chdir("lambda"+str(lamda))
    with open(str(lamda)+'-em.mdp', 'w') as inp: inp.write(em.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-em.mdp', c=str(sim_dir)+'/'+str(LIGGRO), p=str(sim_dir)+'/'+str(LIGTOP), o='enmin.tpr', maxwarn=3)
    gmx.mdrun(s='enmin.tpr', deffnm='enmin', verbose=True)

    with open(str(lamda)+'-nvt.mdp', 'w') as inp: inp.write(nvt.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-nvt.mdp', c='enmin.gro', p=str(sim_dir)+'/'+str(LIGTOP), o='nvt.tpr', maxwarn=3)
    gmx.mdrun(s='nvt.tpr', deffnm='nvt', verbose=True)

    with open(str(lamda)+'-npt.mdp', 'w') as inp: inp.write(npt.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-npt.mdp', c='nvt.gro', p=str(sim_dir)+'/'+str(LIGTOP), o='npt.tpr', maxwarn=3)
    gmx.mdrun(s='npt.tpr', deffnm='npt', verbose=True)

    with open(str(lamda)+'-md.mdp', 'w') as inp: inp.write(md.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-md.mdp', c='npt.gro', p=str(sim_dir)+'/'+str(LIGTOP), o='md.tpr', maxwarn=3)
    gmx.mdrun(s='md.tpr', deffnm='md', verbose=True)
    
    os.system("mv md.xvg md"+str(lamda)+".xvg")

    os.chdir("..")

os.mkdir("analysis")
os.chdir("analysis")
os.system("cp -p ../lambda*/md*xvg .")
os.system("gmx bar -f md*.xvg -o -oi >result.txt")
os.chdir(sim_dir)

# Protein-Ligand complex calculations
os.mkdir("ProteinLigand")
os.chdir("ProteinLigand")

for lamda in range(0, 21):
    os.mkdir("lambda"+str(lamda))
    os.chdir("lambda"+str(lamda))
    with open(str(lamda)+'-em.mdp', 'w') as inp: inp.write(em.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-em.mdp', c=str(sim_dir)+'/'+str(COMPGRO), p=str(sim_dir)+'/'+str(COMPTOP), o='enmin.tpr', maxwarn=3)
    gmx.mdrun(s='enmin.tpr', deffnm='enmin', verbose=True)

    with open(str(lamda)+'-nvt.mdp', 'w') as inp: inp.write(nvt.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-nvt.mdp', c='enmin.gro', p=str(sim_dir)+'/'+str(COMPTOP), o='nvt.tpr', maxwarn=3)
    gmx.mdrun(s='nvt.tpr', deffnm='nvt', verbose=True)

    with open(str(lamda)+'-npt.mdp', 'w') as inp: inp.write(npt.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-npt.mdp', c='nvt.gro', p=str(sim_dir)+'/'+str(COMPTOP), o='npt.tpr', maxwarn=3)
    gmx.mdrun(s='npt.tpr', deffnm='npt', verbose=True)

    with open(str(lamda)+'-md.mdp', 'w') as inp: inp.write(md.format(str(LIG), int(lamda)))
    gmx.grompp(f=str(lamda)+'-md.mdp', c='npt.gro', p=str(sim_dir)+'/'+str(COMPTOP), o='md.tpr', maxwarn=3)
    gmx.mdrun(s='md.tpr', deffnm='md', verbose=True)
    
    os.system("mv md.xvg md"+str(lamda)+".xvg")

    os.chdir("..")

os.mkdir("analysis")
os.chdir("analysis")
os.system("cp -p ../lambda*/md*xvg .")
os.system("gmx bar -f md*.xvg -o -oi >result.txt")
os.chdir(sim_dir)



# calculating final binding free energy: - Eprot + Elig + Erestrain
for line in open(str(sim_dir)+'/restraints.info', 'r'):
    if line.startswith(' '):
        Erestrain = line.split()[1]

for line in open(str(sim_dir)+'/ProteinLigand/analysis/result.txt', 'r'):
    if line.startswith('total'):
        Eprot = line.split()[5]
        recep_std = line.split()[7]

for line in open(str(sim_dir)+'/LigandOnly/analysis/result.txt', 'r'):
    if line.startswith('total'):
        Elig = line.split()[5]
        lig_std = line.split()[7]


final_dg = float(Erestrain) - float(Eprot) + float(Elig)

print("Free energy of binding in kcal/mol: {0}".format(final_dg))

with open('final_dg_binding.txt', 'w') as fe:
    fe.write("Free energy of binding in kcal/mol: {0}".format(final_dg))
    fe.close()
