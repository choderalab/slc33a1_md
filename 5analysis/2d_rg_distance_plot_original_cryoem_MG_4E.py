import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib import rcParams
import mdtraj as md
import numpy as np

# Load the PDB file 
pdb_file = "SLC33A1-GSSG-pose2.pdb"
traj = md.load_pdb(pdb_file)

# Calculate GSSG COM and distance to Leu335 CG
# 1. Select GSSG atoms and calculate COM
gssg_com = md.compute_center_of_mass(traj, select='resname GDS')

# 2. Get Leu335 CG atom coordinates
leu335_cg = traj.topology.select('residue 335 and name CG')[0]
leu335_coords = traj.xyz[:, leu335_cg, :]

# 3. Calculate distance between GSSG COM and Leu335 CG
initial_distance = np.linalg.norm(gssg_com - leu335_coords, axis=1)*10 #convert to Angstroms

# Calculate Radius of Gyration for GSSG
# 1. Select GSSG atoms
gssg_atoms = traj.topology.select('resname GDS')

# 2. Calculate radius of gyration
initial_rg = md.compute_rg(traj.atom_slice(gssg_atoms)) * 10  # Convert to Angstroms

print(f"Initial GSSG-Leu335 distance: {initial_distance[0]:.2f} Å")
print(f"Initial GSSG radius of gyration: {initial_rg[0]:.2f} Å")

# Configure text to be editable in vector formats
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['svg.fonttype'] = 'none'

# Set font family globally
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Reading files with no headers
distances1 = pd.read_csv('distances.csv', header=None)[0]
rg1 = pd.read_csv('rg.csv', header=None)[0]

# Create figure with specific size
width_pts = 500
height_pts = 300
plt.figure(figsize=(width_pts/72, height_pts/72))

# Create KDE plot with normalized density using common_norm=True
kde = sb.kdeplot(x=distances1, y=rg1, fill=True, cmap="Blues", thresh=0, levels=100,
                 cbar=True, 
                 cbar_kws={'label': '2D Probability Density'},
                 common_norm=False)  # This ensures proper normalization

plt.scatter(initial_distance[0], initial_rg[0], color='#ee8353', s=50, label='cryo-EM modeled GSSG')

print(distances1[0])

plt.xlabel('GSSG-Leu335 distance (Å)', fontname='Helvetica')
plt.ylabel('Radius of Gyration (Rg)', fontname='Helvetica')
plt.title('Correlation between Rg and GSSG-Leu335 distance (GSSG pose 2)', fontname='Helvetica')
plt.xlim(distances1.min(), distances1.max())
plt.ylim(rg1.min(), rg1.max())
plt.legend()

# Adjust layout to prevent text cutoff
plt.tight_layout()

# Save in multiple formats
plt.savefig('corr_1.png', dpi=300)  # Regular PNG
plt.savefig('corr_1_4E.pdf',           # Vector PDF
            format='pdf', 
            bbox_inches='tight',
            pad_inches=0.1)
plt.savefig('corr_1.svg',           # Vector SVG
            format='svg',
            bbox_inches='tight',
            pad_inches=0.1)

plt.show()
plt.close()