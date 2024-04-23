import os, sys
import parmed as pmd
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns

warnings.filterwarnings('ignore')

def wrapping(topology, trajectory, trjtype, wrap_traj):
    wrap_protocol = f'''mol new {topology}
mol addfile {trajectory} type {trjtype} first 1 last -1 step 1 filebonds 1 autobonds 1 waitfor all
package require pbctools
pbc wrap -all -compound res -center bb -centersel "protein or nucleic"
animate write dcd {wrap_traj}.dcd
quit'''

    if not os.path.exists('wrap.tcl'):
        with open('wrap.tcl', 'w') as f:
            f.write(wrap_protocol)

    os.system(f'vmd -dispdev text -e wrap.tcl')
    os.system("rm wrap.tcl")
    return


def split_complex(topology, trajectory, ligname):

    u = mda.Universe(topology, trajectory)

    u_complex = u.select_atoms(f'protein or resname {ligname}')
    complex_file = f'complex_only.dcd'

    u_receptor = u.select_atoms('protein')
    receptor_file = f'receptor_only.dcd'

    u_ligand = u.select_atoms(f'resname {ligname}')
    ligand_file = f'ligand_only.dcd'

    u_complex.write(complex_file, frames='all')

    u_receptor.write(receptor_file, frames='all')

    u_ligand.write(ligand_file, frames='all')

    return



def namd_mmgbsa(prmtop, rst7, dcdfile):
    system = prmtop.split('.')[0]
    mmgbsa = f'''amber           on
parmfile       {prmtop}
ambercoor      {rst7} 

outputname gb_bnm

numsteps 0

GBIS on
solventDielectric 74.69
ionConcentration 0.3
alphaCutoff 14
switching on
switchdist 9.
cutoff 10
pairlistdist 11.5

sasa on
surfaceTension 0.0072
timestep 1
nonbondedFreq 1
fullElectFrequency 1
exclude scaled1-4
temperature 310
set ts 0

coorfile open dcd {dcdfile}

while {{ ![coorfile read] }} {{
    firstTimestep $ts
    run 0
    incr ts 1
}}
coorfile close '''

    if not os.path.exists(f'{system}-mmgbsa.conf'):
        with open(f'{system}-mmgbsa.conf', 'w') as f:
            f.write(mmgbsa)

    os.system(f'namd3 +auto-provision {system}-mmgbsa.conf >{system}-mmgbsa.log')
    os.system("rm gb_bnm.*")
     extract_data(f'{system}-mmgbsa.log')

    return 


def extract_data(inpfile):
    out = open('energy_'+inpfile, 'w')
    for line in open(inpfile, 'r').readlines():
        if line.startswith("ENERGY:"):
            time = line.split()[1]
            elec = line.split()[6]
            vdw  = line.split()[7]
            potener = line.split()[13]
            out.write(str(time) +'\t' +str(potener))
            out.write('\n')
        else:
            pass
    out.close()
    return


def binding_free_energy():
    comp = pd.read_csv('complex_only-mmgbsa.log', sep='\s+', names=['Time', 'Elect', 'Vdw', 'Total'])
    rec = pd.read_csv('receptor_only-mmgbsa.log', sep='\s+', names=['Time', 'Elect', 'Vdw', 'Total'])
    lig = pd.read_csv('ligand_only-mmgbsa.log', sep='\s+', names=['Time', 'Elect', 'Vdw', 'Total'])

    bind['Bind_FE'] = comp['Total'] - rec['Total'] - lig['Total']
    bind['Bind_FE_Elect'] = comp['Elect'] - rec['Elect'] - lig['Elect']
    bind['Bind_FE_Vdw'] = comp['Vdw'] - rec['Vdw'] - lig['Vdw']
    final = pd.concate([comp['Time'], bind], keys=['Time', 'Bind_FE', 'Bind_FE_Elect', 'Bind_FE_Vdw'])
    final.to_csv('Final_Bind_FE.txt')

    return 

def main():
    wrapping(topology, trajectory, trjtype, "wrapped_traj")
    split_complex(topology, "wrapped_traj.dcd", ligname)
    namd_mmgbsa("complex_only.prmtop", "complex_only.rst7", "complex_only.dcd")
    namd_mmgbsa("receptor_only.prmtop", "receptor_only.rst7", "receptor_only.dcd")
    namd_mmgbsa("ligand_only.prmtop", "ligand_only.rst7", "ligand_only.dcd")
    binding_free_energy()

    return 

if __name__ == "__main__":
    if len(sys.argv) > 2:
        topology = sys.argv[1]
        trajectory = sys.argv[2]
        ligname = sys.argv[3]
        trjtype = "xtc"
        main()
    else:
        print("Usuage: python3 mm-gbsa-v2.py topology trajectory ligandname")
