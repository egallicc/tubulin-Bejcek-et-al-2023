from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "B-aS"

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds, hydrogenMass= 1.5*amu)
atm_utils = ATMMetaForceUtils(system)

#restrain the positions of C-alpha atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
posrestr_atoms = [ ]
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

temperature = 300 * kelvin

#add barostat
barostat = MonteCarloBarostat(1*bar, temperature)
system.addForce(barostat)
barostat.setFrequency(0)

#set up integrator
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.004 * picosecond
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)
integrator = MTSLangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, [(0,4), (nonbonded_force_group,1)])
integrator.setConstraintTolerance(0.00001)

#platform_name = 'OpenCL'
platform_name = 'HIP'
platform = Platform.getPlatformByName(platform_name)
properties = {}
properties["Precision"] = "mixed"
properties["DeviceIndex"] = "1"
simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

pote = simulation.context.getState(getEnergy = True).getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_npt.xml')

print("Potential Energy =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("MD ...")

nprint = 100000
totalSteps = 25000000
simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, potentialEnergy = True, temperature=True, volume=True))
simulation.reporters.append(DCDReporter(jobname + "_md.dcd", nprint))
simulation.step(totalSteps)

print( "SaveState ...")
simulation.saveState(jobname + '_md.xml')

#save a pdb file that can be used as a topology to load .dcd files in vmd
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
positions = simulation.context.getState(getPositions=True).getPositions()
with open(jobname + '_md.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
