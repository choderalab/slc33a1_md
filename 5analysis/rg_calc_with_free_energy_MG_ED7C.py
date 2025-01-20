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

# Select GSSG atoms2
gssg_atoms2 = trajs_traj_2.topology.select('resname GDS')

# Calculate radius of gyration2
rg2 = md.compute_rg(trajs_traj_2.atom_slice(gssg_atoms2))*10

np.savetxt('rg.csv', rg2, delimiter=',', fmt='%.3f')

# Create first figure - RG KDE plot
plt.figure(figsize=(8, 6), dpi=300)
sb.set_style("ticks")  # Use ticks style without grid
plt.grid(False)  # Ensure no grid is shown
sb.kdeplot(rg2, label='GSSG 2', linewidth=2, color='cornflowerblue')
plt.xlabel('Radius of gyration (Å)', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('rg_comparison_e7C.pdf', format='pdf', bbox_inches='tight')
plt.close()

# Calculate the free energy associated with the radius of gyration
kB = 0.001987  # kcal/(mol*K)
T = 303.15  # temperature in Kelvin (as used in the simulation)

rg_values = rg2
hist, bin_edges = np.histogram(rg_values, bins=35, density=True)
P_x = hist / np.sum(hist)
free_energy = -kB * T * np.log(P_x + 1e-10)  # avoid log(0)

# Create second figure - Free energy plot
plt.figure(figsize=(8, 6), dpi=300)
sb.set_style("ticks")  # Use ticks style without grid
plt.grid(False)  # Ensure no grid is shown
plt.plot(bin_edges[:-1], free_energy, linewidth=2, color='gray')
plt.xlabel("Radius of Gyration (Å)", fontsize=12)
plt.ylabel("Free Energy (kcal/mol)", fontsize=12)
plt.tight_layout()
plt.savefig('rg_free_energy.pdf', format='pdf', bbox_inches='tight')
plt.close()