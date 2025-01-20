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
import matplotlib as mpl

# Set global font to Helvetica (font family 42)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']


# Define individual trajectories for pose 2
def load_trajectories():
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
    
    return [traj2_1,traj2_2,traj2_3,traj2_4,traj2_5,traj2_6,traj2_7,traj2_8,traj2_9,traj2_1,traj2_13]

# Define RMSD computation function
def compute_rmsd(traj, reference=None):
    """
    Compute RMSD for a trajectory
    Args:
        traj: MDTraj trajectory object
        reference: Reference structure (defaults to first frame)
    Returns:
        rmsd: Array of RMSD values in Angstroms
    """
    if reference is None:
        reference = traj[0]
    backbone_indices = traj.topology.select('protein and name CA')
    rmsd = md.rmsd(traj, reference, atom_indices=backbone_indices) * 10  # Convert to Angstroms
    return rmsd

def process_trajectories(trajs):
    """
    Process trajectories and compute RMSD
    Args:
        trajs: List of trajectory objects
    Returns:
        all_rmsd: List of RMSD arrays truncated to same length
    """
    # Print trajectory lengths
    for i, traj in enumerate(trajs):
        print(f"Length of trajectory {i+1}: {len(traj)} frames")
    
    # Find minimum length across all trajectories
    min_length = min(len(traj) for traj in trajs)
    print(f"Truncating all trajectories to {min_length} frames")
    
    # Calculate RMSD for each trajectory and truncate to minimum length
    all_rmsd = []
    for traj in trajs:
        rmsd = compute_rmsd(traj)
        all_rmsd.append(rmsd[:min_length])
    
    return np.array(all_rmsd)

def plot_rmsd(all_rmsd, output_file='gssg1_rmsd'):
    """
    Plot RMSD with confidence intervals
    Args:
        all_rmsd: numpy array of RMSD values
        output_file: filename for saving the plot (without extension)
    """
    # Calculate statistics
    mean_rmsd = np.mean(all_rmsd, axis=0)
    std_error = st.sem(all_rmsd, axis=0)
    ci_95 = 1.96 * std_error
    
    # Create time axis (assuming 1 frame = 1 ns, adjust if different)
    time = np.linspace(0, len(mean_rmsd), len(mean_rmsd))
    
    # Create plot with specific font settings
    plt.figure(figsize=(10, 6))
    
    # Plot data
    plt.plot(time, mean_rmsd, label='Mean RMSD', color='cornflowerblue')
    plt.fill_between(time, 
                     mean_rmsd - ci_95, 
                     mean_rmsd + ci_95, 
                     color='cornflowerblue', 
                     alpha=0.2, 
                     label='95% CI')
    
    # Set labels with specific font
    plt.xlabel('Time (ns)', fontname='Helvetica', fontsize=12)
    plt.ylabel('RMSD (Å)', fontname='Helvetica', fontsize=12)
    plt.title('Average RMSD with 95% Confidence Interval', 
              fontname='Helvetica', 
              fontsize=14)
    
    # Customize tick labels
    plt.tick_params(axis='both', which='major', labelsize=10)
    for tick in plt.gca().get_xticklabels():
        tick.set_fontname('Helvetica')
    for tick in plt.gca().get_yticklabels():
        tick.set_fontname('Helvetica')
    
    # Customize legend
    legend = plt.legend()
    for text in legend.get_texts():
        text.set_fontname('Helvetica')
    
    # Save as PDF (vector format) and PNG
    plt.savefig(f'{output_file}.pdf', 
                dpi=300, 
                bbox_inches='tight', 
                format='pdf')
    plt.savefig(f'{output_file}.png', 
                dpi=300, 
                bbox_inches='tight', 
                format='png')
    
    plt.show()
 

def main():
    # Load trajectories
    print("Loading trajectories...")
    trajs = load_trajectories()
    
    # Join trajectories (if needed)
    trajs_joined = md.join(trajs)
    
    # Process trajectories and compute RMSD
    print("Computing RMSD...")
    all_rmsd = process_trajectories(trajs)
    
    # Plot results
    print("Generating plot...")
    plot_rmsd(all_rmsd)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()