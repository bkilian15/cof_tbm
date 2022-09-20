from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.visualize import view
from build_input_data import *
from ase.io import read
import shutil
from build_molecules import *
import sys
from ase import Atoms, Atom 
from ase.data.pubchem import pubchem_atoms_search, pubchem_atoms_conformer_search

pref = sys.argv[1]
smearing = sys.argv[2]

if smearing == 'y':
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.011
else:
    occupations = 'fixed'
    smearing = None
    degauss = None


pseudo_dir = "/home/kilian/C3MP/simulations/pseudo/"
out_dir = f"/scratch/kilian/out/conf_B/{pref}/"

if pref == "C3N3":
    atoms = make_carbon_ring(ring_atoms=["C","N","C","N","C","N"])      # C3N3
    atoms.center(vacuum=10, axis=(0,1,2))
elif pref == "B3O3":
    atoms = make_carbon_ring(ring_atoms=["B","O","B","O","B","O"])      # B3O3
    atoms.center(vacuum=10, axis=(0,1,2))
elif pref == "C6H4":
    atoms = pubchem_atoms_search(cid=241)
    atoms.center(vacuum=10, axis=(0,1,2))
    print('mol centered')
else:
    path_to_relax_output = f"/scratch/kilian/out/conf_B/{pref}/{pref}_relax.pwo"
    atoms = read(path_to_relax_output, format='espresso-out')

if pref == "benz":
    path_to_relax_output = f"/scratch/kilian/out/conf_B/{pref}/{pref}_relax.pwo"
    atoms = read(path_to_relax_output, format='espresso-out')

    h1 = Atom("H", [6.13492959+(3.409-2.504), 9.31035986+(12.054-11.495),10])
    h2 = Atom("H", [0.06574371+(2.791-3.706),10.11291669+(7.368-7.928),10])
    atoms.append(h1)
    atoms.append(h2)
    atoms.center(vacuum=10, axis=(0,1,2))

    
nbands = int((sum(atoms.get_atomic_numbers())/2) + 10)
natoms = len(atoms)
ntypes = len(set(atoms.get_atomic_numbers()))
control = set_control_parameters(calculation='scf', prefix=pref, pseudo_dir=pseudo_dir, outdir=out_dir)
system = set_system_parameters(ibrav=0, nat=natoms, ntyp=ntypes, nbnd=nbands, ecutwfc=84.0, occupations=occupations, smearing=smearing, degauss=degauss)
electrons = set_electrons_parameters(conv_thr=1.0e-9)
input_data = {'control':control,
              'system':system,
              'electrons':electrons}

pseudopotentials = {}
for e in set(atoms.get_chemical_symbols()):
    pseudopotentials[e] = e+".upf"

calc = Espresso(input_data=input_data, pseudopotentials=pseudopotentials, kpts=(1,1,1))
atoms.calc = calc

pe = atoms.get_potential_energy()
print(pe)

shutil.move(out_dir+"espresso.pwo", out_dir+pref+"_scf.pwo")
shutil.move(out_dir+"espresso.pwi", out_dir+pref+"_scf.pwi")
#shutil.copy(out_dir+pref+".save", out_dir+pref+"_scf.save")
shutil.copy(out_dir+pref+".xml", out_dir+pref+"_scf.xml")
