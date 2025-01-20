# slc33a1_md

This project is focused on preparing and simulating an atomic model of SLC33A1 with gluthatione disulfide (GSSG) bound. This project is associated with "SLC33A1 exports oxidized glutathione to maintain endoplasmic reticulum redox homeostasis" (citation here).

# slc33a1

SLC33A1 is a GSSG transporter which belongs to the major facilitator superfamily (MFS) protein superfamily. SLC33A1 transports GSSG out of the endoplasmic reticulum (ER) and into the cell cytoplasm. This transport activity is critical for maintaining redox homeostasis in the cell. 

The structure of SLC33A1 in an inward-open facing state was determined via cryo-EM. A non-protein density corresponding to GSSG was found in the transport cavity of SLC33A1. The purpose of this project is to use molecular dynamics (MD) simulations to predict the stability of GSSG in the SLC33A1 transport cavity as determined by cryo-EM and to predict the interaction network which may facilitate GSSG transport by SLC33A1. 

The general simulation and analysis scheme for this project is as follows. The structure of SLC33A1 bound with GSSG was placed in a lipid bilayer membrane and water box with neutralizing ions using CHARMMGUI (https://charmm-gui.org/). This system was simulated using OpenMM 7 (https://openmm.org/) using the CHARMM36m forcefield (https://academiccharmm.org/showcase/natmeth_2016_14_71). The conformational flexibility of GSSG was assesed by computing the radius of gyration of the peptide throughout the simulation. The stability of GSSG in the SLC33A1 binding pocket was assesed by computing the distance of the peptide from a nearby protein resiude (leucine 335 in `slc33a1map_gssg_1.pdb`). The interactivity of GSSG with protein residues was assesed by computing the frequency that GSSG atoms came within a cutoff distance of protein residue atoms. 

# Simulation system setup

1. `slc33a1map_gssg_1.pdb` was prepped for simulation using the Schrödinger Maestro Protein Preparation Wizard at pH 7.5 (https://www.schrodinger.com/platform/products/maestro/). The position of the resulting protein structure in a membrane was then calulcated using the Orientations of Proteins in Membranes (OPM) PPM 2.0 webserver (https://opm.phar.umich.edu/ppm_server), resulting in `slc33a1map_gssg_2.pdb`.
2. GSSG in partially protonated state (net charge -2) was parameterized for simulation in the CHARMM36m forcefield using the CHARMMGUI Ligand Reader and Modeler (https://charmm-gui.org/?doc=input/ligandrm). The atom naming for GSSG in the output PDB file was used to replace the atom naming for GSSG in `slc33a1map_gssg_2.pdb`, resulting in `slc33a1map_gssg_3.pdb`.
   

# Directory structure
