import mdtraj as md
import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm
import time

# Load trajectories with progress bar
print("Loading trajectories...")
trajs2 = []
for i in tqdm(range(1,13), desc="Loading DCD files"):
    try:
        traj = md.load_dcd(f'trajectories/{i}.dcd', top='topology.psf')
        trajs2.append(traj)
    except (IOError, FileNotFoundError) as e:
        print(f"Warning: Could not load trajectories/{i}.dcd - File not found")
        continue

print("Joining trajectories...")
trajs_traj_2 = md.join(trajs2)
traj = trajs_traj_2
print(f"Total frames loaded: {len(traj)}")

# GSSG residue atom groups definition
gssg_residues_atoms = {
    'Glu_1': ['O1', 'O2', 'C1', 'C13', "C5", 'C9', 'N1', "C11", "O7"],
    'Cys_1': ['N2', 'O8', 'C2', 'C6', 'S1', 'C14'],
    'Gly_1': ['OXT1', 'O3', 'C15', 'CA1', 'N3'],
    'Glu_2': ['O5', 'O6', 'C18', 'N6', 'C4', 'C8', 'C10', 'C12', 'O10'],
    'Cys_2': ['N5', 'O9', 'C3', 'C17', 'C7', 'S2'],
    'Gly_2': ['O4', 'OXT2', 'C16', 'CA2', 'N4']
}

print("Selecting atoms...")
# Select GSSG atoms and protein atoms
gssg_groups_indices = {
    group: traj.topology.select(f'resname GDS and name {" ".join(atoms)}')
    for group, atoms in gssg_residues_atoms.items()
}
protein_indices = traj.topology.select('protein')

# Analysis parameters
distance_cutoff = 0.4  # 4 Å
interaction_matrix = {group: Counter() for group in gssg_residues_atoms}

# Run analysis with progress bars
print("Running interaction analysis...")
start_time = time.time()
for group_name, gssg_indices in tqdm(gssg_groups_indices.items(), desc="Processing GSSG groups"):
    for frame_idx, frame in enumerate(tqdm(traj, desc=f"Processing frames for {group_name}", leave=False)):
        pairs = md.compute_neighbors(frame, distance_cutoff, gssg_indices, haystack_indices=protein_indices)[0]
        interacting_residues = set(frame.topology.atom(pair).residue.index for pair in pairs)
        interaction_matrix[group_name].update(interacting_residues)
        
        # Print estimated time remaining every 100 frames
        if frame_idx % 100 == 0 and frame_idx > 0:
            elapsed_time = time.time() - start_time
            frames_per_second = frame_idx / elapsed_time
            remaining_frames = len(traj) - frame_idx
            estimated_remaining_time = remaining_frames / frames_per_second
            print(f"\nEstimated time remaining for {group_name}: {estimated_remaining_time/60:.1f} minutes")

# Process results
print("Processing results...")
residue_names = {res.index: str(res) for res in traj.topology.residues}
total_frames = len(traj)

heatmap_data = {}
for group_name, residue_count in interaction_matrix.items():
    heatmap_data[group_name] = {residue_names[res_index]: (count / total_frames) * 100 
                               for res_index, count in residue_count.items()}

df = pd.DataFrame.from_dict(heatmap_data, orient='index').fillna(0)
df = df.T

# Save the results
print("Saving results to CSV...")
df.to_csv('interaction_data_final.csv')
print("Analysis complete!")

total_time = (time.time() - start_time) / 60
print(f"Total analysis time: {total_time:.1f} minutes")