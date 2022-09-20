from ase import Atoms, Atom 
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.visualize import view
from build_input_data import *
from math import *
import numpy as np
import shutil
from build_molecules import *
import sys
import os



pref = sys.argv[1]
cell_dofree_ = sys.argv[2]
smearing = sys.argv[3]

if smearing == 'y':
    occupations='smearing'
    smearing='gaussian'
    degauss=0.011
else:
    occupations=None
    smearing=None
    degauss=None


pseudo_dir = "/home/kilian/C3MP/simulations/pseudo/"
out_dir = f"/scratch/kilian/out/conf_B/{pref}/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# create building blocks
c_block_1 = make_carbon_ring(connections=[None,"H",None,"H",None,"H"])  # C6H3
c_block_2 = make_carbon_ring(ring_atoms=["C","N","C","N","C","N"])      # C3N3
c_block_3 = make_carbon_ring(ring_atoms=["B","O","B","O","B","O"])      # B3O3
c_block_4 = Atoms([Atom("C", [0,0,0])])                                 # C 
c_block_5 = Atoms([Atom("N", [0,0,0])])                                 # N

mol11 = combine_blocks_c(c_block_1, c_block_1)
mol12 = combine_blocks_c(c_block_1, c_block_2)
mol13 = combine_blocks_c(c_block_1, c_block_3)
mol14 = combine_blocks_c(c_block_1, c_block_4)
mol15 = combine_blocks_c(c_block_1, c_block_5)
mol22 = combine_blocks_c(c_block_2, c_block_2)
mol23 = combine_blocks_c(c_block_2, c_block_3)
mol24 = combine_blocks_c(c_block_2, c_block_4)
mol25 = combine_blocks_c(c_block_2, c_block_5)
mol33 = combine_blocks_c(c_block_3, c_block_3)
mol34 = combine_blocks_c(c_block_3, c_block_4)
mol35 = combine_blocks_c(c_block_3, c_block_5)
mol44 = combine_blocks_c(c_block_4, c_block_4)
mol45 = combine_blocks_c(c_block_4, c_block_5)
mol55 = combine_blocks_c(c_block_5, c_block_5)
c_block_2.center(vacuum=10, axis=(0,1,2))
c_block_3.center(vacuum=10, axis=(0,1,2))

mol_dict_c = {'C6H3_C6H3':mol11,'C6H3_C3N3':mol12, 'C6H3_B3O3':mol13,'C6H3_C':mol14,'C6H3_N':mol15,
              'C3N3_C3N3':mol22, 'C3N3_B3O3':mol23, 'C3N3_C':mol24, 'C3N3_N':mol25,
              'B3O3_B3O3':mol33, 'B3O3_C':mol34, 'B3O3_N':mol35,'C_C':mol44, 'C_N':mol45,
              'N_N':mol55, 'C3N3':c_block_2, 'B3O3':c_block_3, 'C':c_block_4, 'N':c_block_5
}


mol_to_simulate = mol_dict_c[pref]
ntypes = len(set(mol_to_simulate.get_atomic_numbers()))

etot_conv_thr = 1.0e-6
forc_conv_thr = 1.0e-4
tprnfor = True
natoms = len(mol_to_simulate)
control = set_control_parameters(calculation='vc-relax', prefix=pref, pseudo_dir=pseudo_dir, outdir=out_dir, etot_conv_thr=etot_conv_thr, forc_conv_thr=forc_conv_thr, tprnfor=tprnfor)
system = set_system_parameters(ibrav=0,  nat=natoms, ntyp=ntypes, ecutwfc=42.0, occupations=occupations, smearing=smearing, degauss=degauss)
electrons = set_electrons_parameters(conv_thr=1.0e-8)
cell = set_cell_parameters(cell_dofree=cell_dofree_)
input_data = {'control':control,
              'system':system,
              'electrons':electrons,
              'cell':cell}

pseudopotentials = {}
for e in set(mol_to_simulate.get_chemical_symbols()):
    pseudopotentials[e] = e+".upf"

calc = Espresso(input_data=input_data, pseudopotentials=pseudopotentials)
mol_to_simulate.calc = calc

pe = mol_to_simulate.get_potential_energy() 
