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

jobname = "5itz-p-C-aS"

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds, hydrogenMass= 1.5*amu)
atm_utils = ATMMetaForceUtils(system)

#restrain the positions of C-alpha atoms
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
posrestr_atoms = [8, 25, 49, 64, 75, 94, 105, 124, 141, 157, 164, 181, 191, 198, 214, 231, 250, 257, 271, 281, 292, 316, 331, 350, 371, 382, 401, 416, 433, 440, 459, 484, 490, 502, 509, 526, 551, 557, 568, 580, 602, 616, 635, 642, 649, 656, 668, 680, 691, 711, 725, 739, 759, 779, 790, 806, 820, 827, 837, 844, 866, 883, 907, 913, 937, 947, 963, 983, 999, 1011, 1030, 1053, 1059, 1073, 1089, 1108, 1120, 1135, 1151, 1175, 1189, 1196, 1210, 1231, 1255, 1272, 1291, 1311, 1336, 1342, 1357, 1374, 1393, 1412, 1426, 1433, 1455, 1470, 1482, 1492, 1502, 1516, 1530, 1551, 1561, 1585, 1592, 1609, 1630, 1644, 1663, 1670, 1692, 1707, 1726, 1745, 1757, 1776, 1792, 1811, 1823, 1847, 1866, 1890, 1912, 1931, 1941, 1953, 1970, 1981, 1995, 2002, 2021, 2038, 2045, 2065, 2084, 2100, 2120, 2138, 2149, 2169, 2176, 2183, 2190, 2204, 2211, 2222, 2229, 2249, 2263, 2274, 2293, 2312, 2329, 2344, 2368, 2387, 2398, 2414, 2426, 2447, 2454, 2476, 2498, 2509, 2531, 2550, 2565, 2585, 2596, 2615, 2644, 2650, 2668, 2674, 2691, 2707, 2718, 2732, 2742, 2758, 2774, 2798, 2804, 2825, 2839, 2850, 2869, 2888, 2902, 2916, 2934, 2948, 2962, 2981, 2996, 3013, 3024, 3036, 3047, 3057, 3077, 3094, 3110, 3122, 3136, 3151, 3161, 3180, 3201, 3213, 3232, 3243, 3267, 3291, 3305, 3324, 3336, 3355, 3370, 3402, 3408, 3422, 3443, 3457, 3471, 3490, 3504, 3528, 3547, 3566, 3577, 3594, 3613, 3629, 3640, 3651, 3670, 3684, 3694, 3705, 3724, 3748, 3768, 3780, 3787, 3797, 3816, 3830, 3846, 3858, 3877, 3891, 3906, 3926, 3943, 3957, 3971, 3990, 4014, 4020, 4049, 4055, 4079, 4098, 4115, 4143, 4149, 4168, 4178, 4192, 4213, 4231, 4237, 4253, 4272, 4283, 4293, 4308, 4330, 4340, 4361, 4378, 4393, 4410, 4429, 4440, 4456, 4466, 4481, 4500, 4514, 4528, 4538, 4549, 4569, 4592, 4598, 4608, 4622, 4639, 4656, 4672, 4694, 4705, 4725, 4731, 4755, 4772, 4779, 4801, 4822, 4839, 4849, 4860, 4871, 4890, 4909, 4930, 4954, 4961, 4973, 4989, 5013, 5019, 5041, 5053, 5069, 5083, 5093, 5103, 5122, 5132, 5146, 5165, 5187, 5201, 5223, 5247, 5258, 5277, 5294, 5314, 5330, 5342, 5366, 5385, 5391, 5405, 5412, 5432, 5454, 5470, 5477, 5496, 5510, 5531, 5556, 5570, 5576, 5590, 5606, 5630, 5636, 5643, 5650, 5662, 5681, 5691, 5713, 5729, 5746, 5770, 5780, 5796, 5807, 5824, 5843, 5854, 5868, 5882, 5896, 5906, 5925, 5935, 5950, 5960, 5984, 5994, 6018, 6037, 6050, 6067, 6089, 6109, 6121, 6140, 6157, 6178, 6188, 6210, 6234, 6244, 6264, 6280, 6297, 6321, 6342, 6358, 6365, 6380, 6387, 6404, 6419, 6434, 6441, 6456, 6476, 6487, 6502, 6512, 6536, 6551, 6563, 6580, 6590, 6600, 6619, 6634, 6656, 6668, 6689, 6704, 6719, 6735, 6742, 6758, 6782, 6799, 6823, 6839, 6858, 6874, 6892, 6911, 6928, 6938, 6945, 6962, 6973, 6980, 6994, 7011, 7030, 7037, 7047, 7069, 7089, 7113, 7128, 7144, 7163, 7174, 7186, 7201, 7218, 7225, 7244, 7264, 7270, 7284, 7291, 7302, 7323, 7340, 7347, 7359, 7370, 7382, 7401, 7418, 7437, 7452, 7476, 7495, 7509, 7525, 7546, 7567, 7581, 7596, 7606, 7620, 7627, 7641, 7663, 7684, 7708, 7714, 7738, 7748, 7767, 7786, 7802, 7814, 7833, 7856, 7862, 7869, 7883, 7900, 7913, 7924, 7940, 7964, 7975, 7990, 7996, 8016, 8023, 8040, 8059, 8079, 8111, 8117, 8129, 8143, 8163, 8179, 8199, 8206, 8223, 8234, 8241, 8251, 8258, 8272, 8286, 8310, 8320, 8342, 8349, 8366, 8387, 8401, 8416, 8423, 8433, 8448, 8467, 8483, 8495, 8506, 8522, 8541, 8553, 8569, 8585, 8609, 8631, 8646, 8657, 8672, 8683, 8694, 8706, 8717, 8736, 8753, 8760, 8780, 8797, 8816, 8830, 8847, 8858, 8877, 8884, 8891, 8898, 8912, 8919, 8930, 8937, 8954, 8961, 8975, 8994, 9013, 9032, 9043, 9065, 9084, 9108, 9123, 9138, 9167, 9173, 9185, 9209, 9228, 9245, 9259, 9273, 9293, 9304, 9320, 9345, 9351, 9370, 9376, 9398, 9414, 9425, 9437, 9451, 9467, 9483, 9507, 9513, 9534, 9548, 9558, 9572, 9591, 9602, 9618, 9635, 9652, 9671, 9687, 9702, 9716, 9730, 9742, 9758, 9772, 9793, 9804, 9823, 9835, 9849, 9864, 9874, 9893, 9914, 9926, 9945, 9956, 9976, 10000, 10014, 10033, 10055, 10074, 10088, 10110, 10116, 10130, 10151, 10158, 10170, 10189, 10203, 10220, 10239, 10255, 10266, 10276, 10290, 10307, 10318, 10325, 10341, 10355, 10369, 10380, 10399, 10423, 10451, 10457, 10464, 10481, 10500, 10514, 10524, 10537, 10556, 10580, 10602, 10621, 10631, 10647, 10661, 10678, 10702, 10708, 10736, 10742, 10766, 10785, 10802, 10822, 10842, 10867, 10873, 10880, 10900, 10918, 10924, 10943, 10957, 10968, 10992, 10999, 11010, 11027, 11044, 11065, 11089, 11099, 11118, 11132, 11156, 11162, 11177, 11196, 11210, 11227, 11244, 11261, 11281, 11293, 11304, 11326, 11340, 11357, 11374, 11384, 11394, 11405, 11425, 11431, 11455, 11473, 11480, 11504, 11525, 11544, 11558, 11574, 11584, 11594, 11613, 11633, 11657, 11664, 11688, 11705, 11716, 11733, 11755, 11770, 11786, 11798, 11813, 11830, 11847, 11866, 11880, 11896, 11913, 11927, 11949, 11963, 11974, 11985, 12006, 12026, 12042, 12057, 12081, 12108, 12114, 12128, 12142, 12158, 12180, 12194, 12204, 12220, 12231, 12243, 12270, 12284, 12290, 12314, 12321, 12340, 12362, 12379, 12390, 12400, 12414, 12434, 12453, 12460, 12474, 12485, 12499, 12509, 12528, 12545, 12560, 12579, 12599, 12621, 12645, 12664, 12675, 12690, 12707, 12727, 12741, 12751, 12768, 12788, 12812, 12836, 12858, 12868, 12888, 12907, 12924, 12948, 12969, 12983, 12990, 13005, 13012, 13029, 13041, 13056, 13073, 13089, 13109, 13123, 13138, 13148, 13163, 13174, 13188, 13205, 13219, 13231, 13250, 13266, 13277, 13292, 13313, 13330, 13347, 13368, 13385, 13397, 13407, 13421, 13443, 13468, 13474, 13494, 13513, 13535, 13559, 13566, 13581, 13588, 13607, 13617, 13641]
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
properties["DeviceIndex"] = "0"

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
