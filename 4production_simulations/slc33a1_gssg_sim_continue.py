#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:26:18 2025

@author: belayv
"""

from sys import stdout
import mdtraj as md 
import openmm as mm
import openmm.app as app
from openmm import LangevinMiddleIntegrator
from openmm.app import CharmmPsfFile, PDBFile, PDBxFile, CharmmParameterSet
import simtk.unit as unit 
import os, time, yaml, bz2
import json

### Open YAML file containing simulation parameters
with open('sim_params.yaml') as file:  # Load YAML containing location of projects
    params = yaml.load(file, Loader=yaml.FullLoader)

### Load PSF Topology:
input_filepath = params['input_paths']['input_filepath']
psf = CharmmPsfFile(input_filepath + 'step5_input.psf')
pdb = PDBFile(input_filepath + 'step5_input.pdb')

### Load parameter files
param_filepath = params['input_paths']['param_filepath']
param_filenames = params['input_paths']['param_filenames']
chk_path=params['input_paths']['chk_path']

param_paths = [os.path.join(param_filepath, file) for file in param_filenames]
params = CharmmParameterSet(*param_paths)

### Set system parameters (load from sysinfo.dat)
sys_info_path = '../input/openmm/sysinfo.dat'
with open(sys_info_path) as sysdata:
    data = json.load(sysdata)
    x, y, z = map(float, data['dimensions'][:3]) * unit.angstroms

psf.setBox(x, y, z)  # Matches box size set by CHARMM

nonbonded_method = app.PME
constraints = app.HBonds
hydrogen_mass = 4.0 * unit.amu
temperature = 303.15 * unit.kelvin
friction = 1 / unit.picosecond
time_step = 0.004 * unit.picoseconds
pressure = 1 * unit.bar
surface_tension = 0

### Set up the system
system = psf.createSystem(
    params,
    nonbondedMethod=nonbonded_method,
    constraints=constraints,
    removeCMMotion=False,
    hydrogenMass=hydrogen_mass,
)

integrator = LangevinMiddleIntegrator(
    temperature,
    friction,
    time_step
)

barostat = mm.MonteCarloMembraneBarostat(
    pressure,
    surface_tension,
    temperature,
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree
)
barostat.setFrequency(50)
system.addForce(barostat)

simulation = app.Simulation(psf.topology, system, integrator)

# Load the checkpoint file
checkpoint_filename = chk_path
if os.path.exists(checkpoint_filename):
    print(f"Loading simulation state from {checkpoint_filename}")
    simulation.loadCheckpoint(checkpoint_filename)
else:
    raise FileNotFoundError(f"Checkpoint file {checkpoint_filename} not found. Ensure the file exists and try again.")

### Simulation parameters
nsteps = 125000000  # 500 ns, 100*10^6 fs
report_freq = 25000  # every 0.1 ns
chk_freq = 12500000  # every 50 ns
traj_freq = 250000  # every 1ns at 4fs / step

current_time = simulation.context.getState().getTime() / unit.nanoseconds
total_simulation_time = nsteps * time_step / unit.nanoseconds
simulation_time = total_simulation_time - current_time
steps_left = round(simulation_time * unit.nanoseconds / time_step)

traj_freq_time = traj_freq * time_step / unit.nanoseconds
report_freq_time = report_freq * time_step / unit.nanoseconds
chk_freq_time = chk_freq * time_step / unit.nanoseconds

### Set output stuff
# Set file names
integrator_xml_filename = "integrator.xml.bz2"
state_xml_filename = "state.xml.bz2"
state_pdb_filename = "final_state.pdb"
state_pdbx_filename = "final_state.cif"
system_xml_filename = "system.xml.bz2"
traj_output_filename = "simulation_output.dcd"

# Write limited state information to standard out:
simulation.reporters.append(
    app.StateDataReporter(
        stdout,
        reportInterval=report_freq,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=steps_left,
        separator="\t",
    )
)

# Write to checkpoint files regularly:
simulation.reporters.append(
    app.CheckpointReporter(
        file=checkpoint_filename,
        reportInterval=chk_freq
    )
)

# Write out the trajectory
simulation.reporters.append(
    md.reporters.DCDReporter(
        file=traj_output_filename,
        reportInterval=traj_freq
    )
)

# Run NPT dynamics
print("Running dynamics in the NPT ensemble...")
initial_time = time.time()
simulation.step(steps_left)
elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = nsteps * time_step
print('    Simulation took %.3f s for %.3f ns (%8.3f ns/day)' % (
    elapsed_time / unit.seconds, simulation_time / unit.nanoseconds, simulation_time / elapsed_time * unit.day / unit.nanoseconds
))

# Save and serialize the final state
print("Serializing state to %s" % state_xml_filename)
state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=True
)

with bz2.open(state_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(state)
    outfile.write(xml)

# Save the final state as a PDB
print("Saving final state as %s" % state_pdb_filename)
with open(state_pdb_filename, "wt") as outfile:
    PDBxFile.writeFile(
        simulation.topology,
        simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True
    )

# Save the final state as a PDBx file (.cif)
print("Saving final state as %s" % state_pdbx_filename)
with open(state_pdbx_filename, 'wt') as outfile:
    PDBxFile.writeFile(
        simulation.topology,
        simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True
    )

# Save and serialize system
print("Serializing system to %s" % system_xml_filename)
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with bz2.open(system_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(simulation.system)
    outfile.write(xml)

# Save and serialize integrator
print("Serializing integrator to %s" % integrator_xml_filename)
with bz2.open(integrator_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(integrator)
    outfile.write(xml)
