# Manifest

To continue simulations following minimization and equilibration:

1. Download CHARMM topology and parameter files from: https://charmm-gui.org/?doc=toppar and ensure that simulation scripts (`*.py` in this directory) can properly point to a directory containing these files. Topology and parameter files for GSSG exist in `/3equilibration/gds`.
2. Run `slc33a1_gssg_sim.py` with the proper dependencies. Note that input and output filenames, and simulation parameters are contained in `sim_params.yaml`. 
