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
# cell_dofree_ = '2Dxy'

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
b_block_1 = make_carbon_ring(connections=[None,"H",None,"H","H","H"])                                         # C6H4
b_block_2 = make_carbon_ring(ring_atoms=["C","N","C","C","C","C"],connections=[None,None,None,"H","H","H"])  # C5NH3
b_block_3 = make_carbon_ring(ring_atoms=["C","C","C","C","N","C"],connections=[None,"H",None,"H",None,"H"])  # C5NH3a
b_block_4 = Atoms([Atom("C", [0,0,0]), Atom("H", [0,l_ch,0])])                                                # CH
b_block_5 = Atoms([Atom("N", [0,0,0])])                                                                      # N
b_block_6 = Atoms([Atom("O", [0,0,0])])                                                                      # O


benz1 = make_carbon_ring(connections=["H","H",None,"H","H","H"])
benz2 = make_carbon_ring(connections=["H","H",None,"H","H","H"])
benzz = combine_blocks_b(benz1, benz2)
benzz.center(vacuum=10, axis=(0,1,2))

mol1 = make_carbon_ring(connections=["H","H","H","H","H","H"])
mol1.center(vacuum=10, axis=(0,1,2))
mol2 = make_carbon_ring(ring_atoms=["C","N","C","C","C","C"],connections=["H",None,"H","H","H","H"])
mol2.center(vacuum=10, axis=(0,1,2))
mol3 = make_carbon_ring(ring_atoms=["C","C","C","C","N","C"],connections=["H","H","H","H",None,"H"])
mol3.center(vacuum=10, axis=(0,1,2))
mol11 = combine_blocks_b(b_block_1, b_block_1)
mol12 = combine_blocks_b(b_block_1, b_block_2)
mol13 = combine_blocks_b(b_block_1, b_block_3)
mol14 = combine_blocks_b(b_block_1, b_block_4)
mol15 = combine_blocks_b(b_block_1, b_block_5)
mol16 = combine_blocks_b(b_block_1, b_block_6)
mol22 = combine_blocks_b(b_block_2, b_block_2)
mol23 = combine_blocks_b(b_block_2, b_block_3)
mol24 = combine_blocks_b(b_block_2, b_block_4)
mol25 = combine_blocks_b(b_block_2, b_block_5)
mol26 = combine_blocks_b(b_block_2, b_block_6)
mol33 = combine_blocks_b(b_block_3, b_block_3)
mol34 = combine_blocks_b(b_block_3, b_block_4)
mol35 = combine_blocks_b(b_block_3, b_block_5)
mol36 = combine_blocks_b(b_block_3, b_block_6)
mol44 = combine_blocks_b(b_block_4, b_block_4)
mol45 = combine_blocks_b(b_block_4, b_block_5)
mol46 = combine_blocks_b(b_block_4, b_block_6)
mol55 = combine_blocks_b(b_block_5, b_block_5)
mol56 = combine_blocks_b(b_block_5, b_block_6)
mol66 = combine_blocks_b(b_block_6, b_block_6)
b_block_4.center(vacuum=10, axis=(0,1,2))
b_block_5.center(vacuum=10, axis=(0,1,2))
b_block_6.center(vacuum=10, axis=(0,1,2))

mol_dict_b = {'C6H4_C6H4':mol11,'C6H4_C5NH3':mol12, 'C6H4_C5NH3a':mol13,'C6H4_CH':mol14,'C6H4_N':mol15,'C6H4_O':mol16,
              'C5NH3_C5NH3':mol22, 'C5NH3_C5NH3a':mol23, 'C5NH3_CH':mol24, 'C5NH3_N':mol25, 'C5NH3_O':mol26,
              'C5NH3a_C5NH3a':mol33, 'C5NH3a_CH':mol34, 'C5NH3a_N':mol35, 'C5NH3a_O':mol36, 'CH_CH':mol44, 'CH_N':mol45,
              'CH_O':mol46, 'N_N':mol55, 'N_O':mol56, 'O_O':mol66, 'C6H4':mol1, 'C5NH3':mol2, 'C5NH3a':mol3, 'CH':b_block_4, 'N':b_block_5, 'O':b_block_6, "benz":benzz
}


mol_to_simulate = mol_dict_b[pref]
ntypes = len(set(mol_to_simulate.get_atomic_numbers()))

etot_conv_thr = 1.0e-6
forc_conv_thr = 1.0e-4
tprnfor = True
natoms = len(mol_to_simulate)
control = set_control_parameters(calculation='relax', prefix=pref, pseudo_dir=pseudo_dir, outdir=out_dir, etot_conv_thr=etot_conv_thr, forc_conv_thr=forc_conv_thr, tprnfor=tprnfor)
system = set_system_parameters(ibrav=0,  nat=natoms, ntyp=ntypes, ecutwfc=60.0, occupations=occupations, smearing=smearing, degauss=degauss)
electrons = set_electrons_parameters(conv_thr=1.0e-8, mixing_beta=.2)
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


#calc.calculate(mol_to_simulate)
te = mol_to_simulate.get_total_energy()
