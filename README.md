# KMC-code-nanoporous-Au
Code for simulating the selective dissolution of Ag from AuAg alloys

This repository implements a kinetic Monte Carlo (KMC) + molecular dynamics (MD) framework to simulate dealloying of bimetallic nanoparticle.
The model incorporates diffusion, dissolution, defects, pore formation, ligament formation and structural collapse as seen during dealloying process.
The code is intended for researchers working on simulating dealloying and nanoporous materials.

############# Features #####################
- During KMC move diffusion and dissolution moves are incorporated.
- MD takes care of structural collapse and ligament formation.


############# Operations #####################
The KMC-MD is implemented in 5 steps:
Step 1: The 18 nm NP is copied into the folder and KMC is performed on the particle. Here Ag is allowed to surface diffuse and dissolve. While Au can surface diffuse.
Step 2: post KMC structure needs to be relaxed. This performed using a short MD.
Step 3: Solvent is added during this step. This is perfromed using MATLAB. Solvent is added all around the relaxed nanoaprticle.
Step 4: After the plain addition of solvent in previous step, the solvent needs to be brought in equilibrium along with the nanoaprticle. So during this step another short MD step where temperture is slowly relaxed and slovent-nanoaprticle system is brought to equilibrium.
Step 5: The MD step of 1 ns is performed. This time is sufficient for the nanoaprticle to relax which causes a collapse of structure. Collapse leads to defect formation and we see ligaments appearing. 

######### Files and Folders ####
KMC_run (folder): contains all the files necessary to run the KMC step
all_scripts (folder): contains all the files to perform the complete operation of KMC-MD. Necessary files to perform stepwise operations one after another.
make-directories.sh (file): creates the necessary stepwise folders and copies required files to that folder and runs the operations from that folder.
18nm_trunc_octa (file): contains the coordinates of truncated octahedron nanoparticle of size 18 nm.
