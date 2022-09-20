from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.visualize import view
from ase.io.espresso import write_espresso_in
from build_input_data import *
from ase.io import read
import shutil
import os
import sys

pref = sys.argv[1]
nbands = int(sys.argv[2])
smearing = sys.argv[3]
conf = "B" # !!!! change

if not os.path.exists(f"conf_{conf}/{pref}/"):
    os.makedirs(f"conf_{conf}/{pref}/")

if smearing == 'y':
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.011
else:
    occupations = 'fixed'
    smearing = None
    degauss = None


pseudo_dir = "/home/kilian/C3MP/simulations/pseudo/"
out_dir = f"/scratch/kilian/out/conf_{conf}/{pref}/"

path_to_relax_output = f"/scratch/kilian/out/conf_{conf}/{pref}/{pref}_relax.pwo"
atoms = read(path_to_relax_output, format='espresso-out')

natoms = len(atoms)
ntypes = len(set(atoms.get_atomic_numbers()))
control = set_control_parameters(calculation='nscf', prefix=pref, pseudo_dir=pseudo_dir, outdir=out_dir, restart_mode='from_scratch')
system = set_system_parameters(ibrav=0, nat=natoms, ntyp=ntypes, nbnd=nbands, ecutwfc=42.0, occupations=occupations, smearing=smearing, degauss=degauss)
electrons = set_electrons_parameters(conv_thr=1.0e-10, diago_full_acc=True)
input_data = {'control':control,
              'system':system,
              'electrons':electrons}

pseudopotentials = {}
for e in set(atoms.get_chemical_symbols()):
    pseudopotentials[e] = e+".upf"


kpts=(7,7,7)
calc = Espresso(input_data=input_data, pseudopotentials=pseudopotentials, kpts=kpts)
atoms.calc = calc

infile = f"conf_{conf}/{pref}/{pref}_nscf.pwi"
with open(infile, 'w') as myfile:
    write_espresso_in(myfile, atoms, input_data=input_data, pseudopotentials=pseudopotentials, kpts=kpts)

## fix kpoints
with open(infile) as file:
        filedata = file.read()

real_kpts = os.popen('/work/c3mp/kilian/kmesh.pl 5 5 5').read()
filedata = filedata.replace("K_POINTS automatic",'')
filedata = filedata.replace("7 7 7  0 0 0",real_kpts)

with open(infile, 'w') as file:
    file.write(filedata)

#shutil.move(out_dir+"espresso.pwo", out_dir+pref+"_scf.pwo")
#shutil.move(out_dir+"espresso.pwi", out_dir+pref+"_scf.pwi")
#shutil.copy(out_dir+pref+".save", out_dir+pref+"_scf.save")
#shutil.copy(out_dir+pref+".xml", out_dir+pref+"_scf.xml")
