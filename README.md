# slc33a1_md

This project is focused on preparing and simulating an atomic model of SLC33A1 with gluthatione disulfide (GSSG) bound. 

# slc33a1

SLC33A1 is a GSSG transporter which belongs to the major facilitator superfamily (MFS) protein superfamily. SLC33A1 transports GSSG out of the endoplasmic reticulum (ER) and into the cell cytoplasm. This transport activity is critical for maintaining redox homeostasis in the cell. 

The structure of SLC33A1 in an inward-open facing state was determined via cryo-EM. A non-protein density corresponding to GSSG was found in the transport cavity of SLC33A1. The purpose of this project is to use molecular dynamics (MD) simulations to predict the stability of GSSG in the SLC33A1 transport cavity as determined by cryo-EM and to predict the interaction network which may facilitate GSSG transport by SLC33A1. 

The general simulation and analysis scheme for this project is as follows. The structure of SLC33A1 bound with GSSG was placed in a lipid bilayer membrane and water box with neutralizing ions using CHARMMGUI (https://charmm-gui.org/). This system was simulated using OpenMM 7 (https://openmm.org/) using the CHARMM36m forcefield (https://academiccharmm.org/showcase/natmeth_2016_14_71). The conformational flexibility of GSSG was assesed by computing the radius of gyration of the peptide throughout the simulation. The stability of GSSG in the SLC33A1 binding pocket was assesed by computing the distance of the peptide from a nearby protein resiude (leucine 335 in `slc33a1map_gssg_1.pdb`).

# Simulation system setup

# Directory structure
