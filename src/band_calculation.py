from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.visualize import view
from build_input_data import *
from ase.io import read
import shutil
from build_molecules import *
import sys
from k_path import *

pref = sys.argv[1]
shape = sys.argv[2]
smearing = sys.argv[3]

if smearing == 'y':
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.011
else:
    occupations = 'fixed'
    smearing = None
    degauss = None

pseudo_dir = "/home/kilian/C3MP/simulations/pseudo/"
out_dir = f"/scratch/kilian/out/conf_B/{pref}_test/"

path_to_relax_output = f"/scratch/kilian/out/conf_B/{pref}/{pref}_relax.pwo"
atoms = read(path_to_relax_output, format='espresso-out')

nbands = int((sum(atoms.get_atomic_numbers())/2) + 10)
natoms = len(atoms)
ntypes = len(set(atoms.get_atomic_numbers()))

control = set_control_parameters(calculation='bands', prefix=pref, pseudo_dir=pseudo_dir, outdir=out_dir, verbosity='high')
system = set_system_parameters(ibrav=0, nat=natoms, ntyp=ntypes, nbnd=nbands, ecutwfc=42.0,  occupations=occupations, smearing=smearing, degauss=degauss)
electrons = set_electrons_parameters(conv_thr=1.0e-7)
input_data = {'control':control,
              'system':system,
              'electrons':electrons}

pseudopotentials = {}
for e in set(atoms.get_chemical_symbols()):
    pseudopotentials[e] = e+".upf"

calc = Espresso(input_data=input_data, pseudopotentials=pseudopotentials)
atoms.calc = calc

lat = atoms.cell.get_bravais_lattice()
pts_list = lat.get_special_points()
print(pts_list)
if shape == 'crr':
    pts = 'GYHG'
elif shape == 'brr':
    pts = "GX"
else:
    pts = ''.join([str(item) for item in pts_list])
print('pints', pts)
path = atoms.cell.bandpath(pts, npoints=100)
print(path)
calc.set(kpts=path)

calc.calculate(atoms)
