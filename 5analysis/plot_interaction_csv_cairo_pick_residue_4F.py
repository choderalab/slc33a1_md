import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
from matplotlib import rcParams
import re

# Configure text to be editable in vector formats
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['svg.fonttype'] = 'none'

# Set font family globally
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Load the saved data
df = pd.read_csv('interaction_data_final.csv', index_col=0)

# Function to extract number from residue name
def get_residue_number(residue_name):
    return int(re.search(r'\d+', residue_name).group())

# Select specific residues
selected_residues = ['TYR225', 'ASN229', 'GLN356', 'TYR390', 'TYR421']

# Sort the selected residues by their numbers in descending order
sorted_residues = sorted(selected_residues, key=get_residue_number, reverse=True)

# Get the filtered DataFrame with sorted residues
filtered_df = df.loc[sorted_residues]

# Create plot with narrow width and tall height
width_pts = 300
height_pts = 300
plt.figure(figsize=(width_pts/72, height_pts/72))

# Create meshgrid for pcolormesh
x = np.arange(filtered_df.shape[1] + 1)
y = np.arange(filtered_df.shape[0] + 1)
X, Y = np.meshgrid(x, y)

# Use pcolormesh instead of imshow
plt.pcolormesh(X, Y, filtered_df, cmap='Blues', shading='auto')
colorbar = plt.colorbar(label='Interaction Frequency (%)')
colorbar.ax.set_ylabel('Interaction Frequency (%)', fontname='Helvetica')

# Adjust tick positions for pcolormesh
plt.xticks(ticks=np.arange(len(filtered_df.columns)) + 0.5, 
          labels=filtered_df.columns, 
          rotation=90, 
          fontname='Helvetica')
plt.yticks(ticks=np.arange(len(filtered_df.index)) + 0.5, 
          labels=filtered_df.index, 
          fontname='Helvetica')

plt.title('Residue Interaction Likelihood with GSSG Groups', fontname='Helvetica')
plt.xlabel('GSSG Groups', fontname='Helvetica')
plt.ylabel('Protein Residues (Residue Name + Index)', fontname='Helvetica')

# Adjust layout to prevent text cutoff
plt.tight_layout()

# Regular PNG export
plt.savefig('gssg1_res_interaction.png', dpi=300)

# Vector graphics export with editable text
plt.savefig('gssg1_res_interaction_specific_residue.pdf', 
            format='pdf', 
            bbox_inches='tight',
            pad_inches=0.1)

plt.savefig('gssg1_res_interaction_specific_residue.svg', 
            format='svg',
            bbox_inches='tight',
            pad_inches=0.1)

plt.close()