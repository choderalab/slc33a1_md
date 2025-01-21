# slc33a1_md

This project is focused on preparing and simulating an atomic model of SLC33A1 with gluthatione disulfide (GSSG) bound. This project is associated with "SLC33A1 exports oxidized glutathione to maintain endoplasmic reticulum redox homeostasis" (citation here).

# slc33a1

SLC33A1 is a GSSG transporter which belongs to the major facilitator superfamily (MFS) protein superfamily. SLC33A1 transports GSSG out of the endoplasmic reticulum (ER) and into the cell cytoplasm. This transport activity is critical for maintaining redox homeostasis in the cell. 

The structure of SLC33A1 in an inward-open facing state was determined via cryo-EM. A non-protein density corresponding to GSSG was found in the transport cavity of SLC33A1. The purpose of this project is to use molecular dynamics (MD) simulations to predict the stability of GSSG in the SLC33A1 transport cavity as determined by cryo-EM and to predict the interaction network which may facilitate GSSG transport by SLC33A1. 

The general simulation and analysis scheme for this project is as follows. The structure of SLC33A1 bound with GSSG was placed in a lipid bilayer membrane and water box with neutralizing ions using CHARMMGUI (https://charmm-gui.org/). This system was simulated using OpenMM 7 (https://openmm.org/) using the CHARMM36m forcefield (https://academiccharmm.org/showcase/natmeth_2016_14_71). The conformational flexibility of GSSG was assesed by computing the radius of gyration of the peptide throughout the simulation. The stability of GSSG in the SLC33A1 binding pocket was assesed by computing the distance of the peptide from a nearby protein resiude (leucine 335 in `slc33a1map_gssg_1.pdb`). The interactivity of GSSG with protein residues was assesed by computing the frequency that GSSG atoms came within a cutoff distance of protein residue atoms. 

# Simulation system setup

1. `/1starting_models/slc33a1map_gssg_1.pdb` was prepped for simulation using the Schrödinger Maestro Protein Preparation Wizard at pH 7.5 (https://www.schrodinger.com/platform/products/maestro/). The position of the resulting protein structure in a membrane was then calulcated using the Orientations of Proteins in Membranes (OPM) PPM 2.0 webserver (https://opm.phar.umich.edu/ppm_server), resulting in `/2model_preparation/slc33a1map_gssg_2.pdb`.
2. GSSG in partially protonated state (net charge -2) was parameterized for simulation in the CHARMM36m forcefield using the CHARMMGUI Ligand Reader and Modeler (https://charmm-gui.org/?doc=input/ligandrm). The atom naming for GSSG in the output PDB file was used to replace the atom naming for GSSG in `/2model_preparation/slc33a1map_gssg_2.pdb`, resulting in `/2model_preparation/slc33a1map_gssg_3.pdb`.
3. `/2model_preparation/slc33a1map_gssg_3.pdb` was used the input for the CHARMMGUI Bilayer Builder (https://charmm-gui.org/?doc=input/membrane.bilayer). The protein-ligand system was placed in a POPC:POPE = 3:1 membrane bilayer. The system was then solvated in a water box with 0.15 M KCl, resulting in `/3equilibration/openmm/step5_input.pdb`. 
   

# Directory structure

* `/1starting_models` contains the PDB file of the atomic model of SLC33A1 bound with GSSG as determined via cryo-EM.
* `/2model_preparation` contains PDB files that are outputs from Maestro, OPM PPM 2.0, and CHARMMGUI Ligand Reader and Modeler.
* `/3equilibration` contains outputs from CHARMMGUI Bilayer Builder and scripts to energy minimize and equilibrate the system using OpenMM.
* `/4production` contains scripts to simulate the system using OpenMM.
* `/5analysis` contains all scripts used to analyze simulation trajectories and to make figures for "SLC33A1 exports oxidized glutathione to maintain endoplasmic reticulum redox homeostasis" (citation here).

# Citations

* This work: Please cite (citation here)
* CHARMMGUI:
  1. Jo, S., Kim, T., Iyer, V. G. & Im, W. CHARMM-GUI: A web-based graphical user interface for CHARMM. J. Comput. Chem. 29, 1859–1865 (2008).
  2. Kim, S. et al. CHARMM-GUI ligand reader and modeler for CHARMM force field generation of small molecules. J. Comput. Chem. 38, 1879–1886 (2017).
  3. Jo, S., Kim, T. & Im, W. Automated Builder and Database of Protein/Membrane Complexes for Molecular Dynamics Simulations. PLOS ONE 2, e880 (2007).
  4. Lee, J. et al. CHARMM-GUI Input Generator for NAMD, GROMACS, AMBER, OpenMM, and CHARMM/OpenMM Simulations Using the CHARMM36 Additive Force Field. J. Chem. Theory Comput. 12, 405–413 (2016).
* OPM:
  1. Lomize, M. A., Pogozheva, I. D., Joo, H., Mosberg, H. I. & Lomize, A. L. OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Res. 40, D370–D376 (2012).
* CHARMM forcefield:
  1. Huang, J. et al. CHARMM36m: an improved force field for folded and intrinsically disordered proteins. Nat. Methods 14, 71–73 (2017).
* OpenMM:
  1. Eastman, P. et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLOS Comput. Biol. 13, e1005659 (2017).
* Schrödinger Maestro Protein Preparation Wizard:
  1. Madhavi Sastry, G., Adzhigirey, M., Day, T., Annabhimoju, R. & Sherman, W. Protein and ligand preparation: parameters, protocols, and influence on virtual screening enrichments. J. Comput. Aided Mol. Des. 27, 221–234 (2013).
