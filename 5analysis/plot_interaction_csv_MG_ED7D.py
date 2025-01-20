import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
import matplotlib

# Set font family globally
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['font.size'] = 6  # Set global font size to 6

# Load the saved data
df = pd.read_csv('interaction_data_final.csv', index_col=0)

# Function to adjust residue numbers
def adjust_residue_label(label):
    # Split the label into residue name and number
    residue_name = ''.join(filter(str.isalpha, label))
    residue_num = int(''.join(filter(str.isdigit, label)))
    
    # Adjust numbers after 404
    if residue_num > 404:
        residue_num -= 3
    
    # Return the adjusted label
    return f"{residue_name}{residue_num}"

# Function to extract number from residue label
def get_residue_number(label):
    return int(''.join(filter(str.isdigit, label)))

# Function to format column labels (remove underscores)
def format_column_label(label):
    return label.replace('_', '')

# Adjust the index labels
df.index = [adjust_residue_label(idx) for idx in df.index]

# Create a list of tuples with (index, number) for sorting
index_with_numbers = [(idx, get_residue_number(idx)) for idx in df.index]
sorted_indices = [idx for idx, _ in sorted(index_with_numbers, key=lambda x: x[1])]

# Reorder the DataFrame rows according to the sorted indices
df = df.loc[sorted_indices]

# Convert pts to inches
width_pts = 300
height_pts = 750
width_inches = width_pts / 72
height_inches = height_pts / 72

# Filter out residues that never reach the cutoff
interaction_cutoff = 0
filtered_df = df[df.ge(interaction_cutoff).any(axis=1)]
filtered_df = filtered_df.where(filtered_df >= interaction_cutoff)

# Create plot with dimensions in pts
matplotlib.use('cairo')

plt.figure(figsize=(width_inches, height_inches))
plt.imshow(filtered_df, cmap='Blues', aspect='auto')
colorbar = plt.colorbar(label='Interaction Frequency (%)')

# Set font properties for colorbar label
colorbar.ax.set_ylabel('Interaction Frequency (%)', 
                      fontname='Helvetica', 
                      fontsize=6)
# Set font size for colorbar tick labels
colorbar.ax.tick_params(labelsize=6)

# Format column labels by removing underscores
formatted_column_labels = [format_column_label(col) for col in filtered_df.columns]

# Set font properties for x and y axis labels and ticks
plt.xticks(ticks=np.arange(len(filtered_df.columns)), 
           labels=formatted_column_labels, 
           rotation=90, 
           fontname='Helvetica', 
           fontsize=6)
plt.yticks(ticks=np.arange(len(filtered_df.index)), 
           labels=filtered_df.index, 
           fontname='Helvetica', 
           fontsize=6)

# Set font properties for title and axis labels
plt.title('Residue Interaction Likelihood with GSSG Groups', 
          fontname='Helvetica', 
          fontsize=6)
plt.xlabel('GSSG Groups', 
          fontname='Helvetica', 
          fontsize=6)
plt.ylabel('Protein Residues (Residue Name + Index)', 
          fontname='Helvetica', 
          fontsize=6)

# Adjust layout to prevent text cutoff
plt.tight_layout()

# Save the figure in different formats
plt.savefig('gssg1_res_interaction.png', dpi=300)
plt.savefig('gssg1_res_interaction.pdf', format='pdf')
plt.savefig('gssg1_res_interaction.svg', format='svg')

plt.close()