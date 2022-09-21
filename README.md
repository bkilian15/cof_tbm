# cof_tbm
Tight-binding model for covalent organic frameworks. Here there is the work completed by Benjamin Kilian during the masters thesis at EPFL, Spring 2022. For any questions, contact Benjamin Kilian (benjamin.kilian@epfl.ch, note this email will expire after Jan. 2023), Marta Brzezinska (marta.brzezinska@epfl.ch) or Prof. Oleg Yazyev (oleg.yazyev@epfl.ch).

## existing data
There is too much data to put it all here on github, but it should be found in the work directory of the fidis machine. it should be in /work/c3mp/kilian. From there, there are several sub-directories which I saved at different times. The relevant one should be "out_27jun" as this is the last checkpoint I made and this is the data I ran further calculations on.

In this /work/c3mp/kilian/out_27jun directory there are 2 relevant directories (plus one 'test' directory), conf_B and conf_C. The conf_B directory corresponds to the case of 1-Dimensional periodicity and the conf_C is the 2-Dimensional periodicity. The subdirectories inside conf_B and conf_C correspond to the materials themselves (e.g. C6H4_C5NH3 corresponds to one building block of a benzene ring and the other a benzene ring with a nitrogen replacing one CH). These directories will have the output from all of the calculations that were run on the material. Note that some of these are incomplete. If I remember correctly, we had some problems especially with the simple-simple cases, i.e. the case where both building blocks are single atoms (or CH). All of the ring-ring configurations should have band structure calculations completed

## pipeline: config files and submission scripts.
The general pipeline we followed on the materials was as follows: relaxation, SCF, band structure calculation, PDOS, vacuum potential calculation, wannierization (though no formal conclusion was reached for the latter and we therefore don't have any real data).

Note that for all of the .run submission files you have to change the tag corresponding to your material manually in the header or come up with some strategy to do this automatically. Also, you'll have to change my absolute paths in my files to whichever paths you are using on your machine.
### Relaxation
I used the src/relax_molecule.py script to create my structures and then run the calculation using ase, and the submission script src/submit_relax.run to submit this to fidis.

### SCF
As before, I use src/scf_calculation.py to run calculation with ase, and src/submit_scf.run to submit to fidis.

### band structure
same as before, with src/band_calculation.py and src/submit_bands.run.

Once this is completed, see src/in/submit_band_plots.run to get your output files to be able to plot the band structure

### PDOS
see /src/pdos/submit_pdos.run. For this you'll need to create config files for your materials, i'll upload mine as examples shortly.

### vacuum potential
these submission files can be found in the src/avg folder. Example config files coming shortly (email me if they're still not there when you need them).

## jupyter notebooks for analysis
In the folder src/ can be found most of the jupyter notebooks that I used for analysis throughout the course of the project. src/plot_bands.ipynb and src/plot_bands_pdos.ipynb are two notebooks used for plotting the band structure of a material after DFT calculations have been run, and plotting the band structure with partial density of states, respectively. In each case I added several helper functions to make the plotting process easier and to give the user more control. In certain cases there may be some things that are hardcoded and can be adjusted to be passed as arguments to functions.

The notebook src/plot_wannier.ipynb is used to confirm that wannierization of a material's Bloch states yields the same bands as the DFT calculations did.

The rest of the jupyter notebooks in src are mainly used for scratch work, but I will leave them here as some of these ideas may be useful to you in the future.


If there seem to be any missing files, please contact me (Benjamin Kilian) and I will try my best to help you out. Good luck
