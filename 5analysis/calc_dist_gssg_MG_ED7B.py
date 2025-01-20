import mdtraj as md 
from matplotlib import pyplot as plt
import numpy as np
import pymbar as pr
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import sys
import random
import seaborn as sb
import scipy.stats as st
from collections import Counter
import matplotlib
matplotlib.use('pdf')  # Set the backend to PDF for vector graphics

# Set up the font configuration
plt.rcParams['pdf.fonttype'] = 42  # Ensure text is stored as true text in PDF
plt.rcParams['ps.fonttype'] = 42   # Ensure text is stored as true text in PostScript
plt.rcParams['font.size'] = 12     # Base font size
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']  # Use Helvetica font

# Define individual trajectories for 2
traj2_1=md.load('trajectories/1.dcd',top='topology.psf')
traj2_2=md.load('trajectories/2.dcd',top='topology.psf')
traj2_3=md.load('trajectories/3.dcd',top='topology.psf')
traj2_4=md.load('trajectories/4.dcd',top='topology.psf')
traj2_5=md.load('trajectories/5.dcd',top='topology.psf')
traj2_6=md.load('trajectories/6.dcd',top='topology.psf')
traj2_7=md.load('trajectories/7.dcd',top='topology.psf')
traj2_8=md.load('trajectories/8.dcd',top='topology.psf')
traj2_9=md.load('trajectories/9.dcd',top='topology.psf')
traj2_12=md.load('trajectories/12.dcd',top='topology.psf')
traj2_13=md.load('trajectories/13.dcd',top='topology.psf')

trajs2=[traj2_1,traj2_2,traj2_3,traj2_4,traj2_5,traj2_6,traj2_7,traj2_8,traj2_9,traj2_1,traj2_13]

# Join trajs 2
trajs_traj_2=md.join(trajs2)

# Calculate center of mass of GSSG 
gssg1_com2=md.compute_center_of_mass(trajs_traj_2,select='resname GDS')
leu335_2 = trajs_traj_2.topology.select('residue 335 and name CG')[0]

# Define the atomic coordinates of the CG of Leu335 (centrally located carbon)
leu335_2 = trajs_traj_2.topology.select('residue 335 and name CG')[0]
protein_atom_coords2 = trajs_traj_2.xyz[:, leu335_2, :]

# Calculate distances between GSSG COM and Leu335 CG for all frames
distances2 = np.linalg.norm(gssg1_com2 - protein_atom_coords2, axis=1)*10  # Euclidean distance

# Create figure with high-quality settings
plt.figure(figsize=(8, 6), dpi=300)

# Create the plot with seaborn
sb.set_style("ticks")  # Use ticks style without grid
plt.grid(False)  # Ensure no grid is shown
sb.kdeplot(distances2, linewidth=2, color='cornflowerblue', label='GSSG')

# Customize the plot
plt.xlabel('GSSG - Leu335 COM distance (Å)', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.legend(fontsize=10)

# Adjust layout to prevent label clipping
plt.tight_layout()

# Save as PDF with vector graphics
plt.savefig('dist_comparison.pdf', format='pdf', bbox_inches='tight')

# Save distances to CSV
np.savetxt('distances.csv', distances2, delimiter=',', fmt='%.3f')

# Close the figure to free memory
plt.close()