
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil
from dmsreader import *
from datetime import datetime
from BEDAMplugin import *

def binding_energy(dms,simulation,positions):
    context = simulation.context
    #energy of bound system
    simulation.context.setPositions(positions)
    state=context.getState(getEnergy=True)
    U1=state.getPotentialEnergy()
    print("U1=%f" % U1.value_in_unit(kilocalorie_per_mole))
    #displace ligand
    new_positions = dms.displaceLigand(positions)
    simulation.context.setPositions(new_positions)
    state=context.getState(getEnergy=True)
    #energy of unbound energy
    U0=state.getPotentialEnergy()
    print("U0=%f" % U0.value_in_unit(kilocalorie_per_mole))
    #simulation.context.setPositions(positions)
    return U1-U0

print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = 'test.out'
f = open(binding_file, 'w')

lmbd = 0.1
print("lambda = %f" %(lmbd))
#rcptfile_input  = 'bcd_recpt_c.dms'
#ligfile_input   = 'benzene_lig_c.dms'
#rcptfile_output = 'bcd_recpt_tmp.dms'
#ligfile_output  = 'benzene_lig_tmp.dms'

rcptfile_input  = 'mcl1-23_rcpt_2.dms'
ligfile_input   = 'mcl1-23_lig_2.dms'
rcptfile_output = 'mcl1-23_rcpt_tmp.dms'
ligfile_output  = 'mcl1-23_lig_tmp.dms'


shutil.copyfile(rcptfile_input, rcptfile_output)
shutil.copyfile(ligfile_input, ligfile_output)
testDes = DesmondDMSFile(ligfile_output, rcptfile_output, BEDAM=True) 

system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.2*nanometer, OPLS = True, implicitSolvent='AGBNP', AGBNPVersion=1)

temperature = 300.0 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

#mcl1-23
natoms_ligand = 43
lig_atom_restr = 12   #index of ligand atom for Vsite (267CA -- central C of ligand)
rcpt_atom_restr = 1535   #index of rcpt atom center of Vsite#267CA -- central C of ligand
kf = 1225.0 #force constant for Vsite in (kj/mol)/nm^2
r0 = 0.7 #radius of Vsite sphere

#benzene
#natoms_ligand = 12
#rcpt_atom_restr = 39  #index of rcpt atom center of Vsite#267CA -- central C of ligand
#lig_atom_restr = 147  #index of ligand atom for Vsite (267CA -- central C of ligand)
#kf = 1000 #1225.0 #force constant for Vsite in (kj/mol)/nm^2
#r0 = 0.8 #0.5 #radius of Vsite sphere


integrator = LangevinIntegratorBEDAM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond,natoms_ligand,lmbd,lig_atom_restr, rcpt_atom_restr,kf, r0)
#integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

#platform = Platform.getPlatformByName('Reference')
platform = Platform.getPlatformByName('OpenCL')

properties = {}
device = "1:0"
m = re.match("(\d+):(\d+)", device)
if m:
    platformid = m.group(1)
    deviceid = m.group(2)
    properties["OpenCLPlatformIndex"] = platformid
    properties["DeviceIndex"] = deviceid
    print("Using platform id: %s, device id: %s" % ( platformid ,  deviceid) )

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
print "Using platform %s" % simulation.context.getPlatform().getName()
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)
context=simulation.context

#state=context.getState(getEnergy=True,getForces=True)
#simulation.reporters.append(StateDataReporter(stdout, 10, step=True, temperature=True, potentialEnergy = True, totalEnergy = True))
#simulation.step(3000)
#U1=state.getPotentialEnergy()
#print("U1=%f" % U1.value_in_unit(kilocalorie_per_mole))
#forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
#at = 0;
#for f in forces:
#    print("F: %d %f %f %f" % (at,f[0],f[1],f[2]))
#    at += 1

#simulation context to compute binding energies without image and without a cutoff (large cutoff)

dms_bindingE = DesmondDMSFile(ligfile_output, rcptfile_output)
system_bindingE = dms_bindingE.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=20*nanometer, OPLS = True, implicitSolvent='AGBNP', AGBNPVersion=1)
#platform_bindingE = Platform.getPlatformByName('Reference')
#platform_properties = {}
integrator_bindingE = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
simulation_bindingE = Simulation(dms_bindingE.topology, system_bindingE, integrator_bindingE, platform, properties)
natoms_bindingE = len(dms_bindingE.positions)
print("Number of receptor atoms:"+str(natoms_bindingE - natoms_ligand))
print("Number of ligand atoms:"+str(natoms_ligand))
print("Number of atoms in complex:"+str(natoms_bindingE))
print("Number of atoms in BEDAM system:"+str(2*natoms_bindingE))

#compute initial binding energy (optional)
positions = simulation.context.getState(getPositions=True).getPositions()
be = binding_energy(dms_bindingE,simulation_bindingE,positions[0:natoms_bindingE])
print("initial binding energy="+str(be))

stepId = 1000
totalSteps = 50000
loopStep = totalSteps/stepId
simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, temperature=True))

for i in range(loopStep):
    start_sim = datetime.now()
    simulation.step(stepId)
    end_sim = datetime.now()
    elapsed=end_sim - start_sim
    print("elapsed simulation time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
    if simulation.currentStep%stepId==0:
        positions = simulation.context.getState(getPositions=True).getPositions()
    	be = binding_energy(dms_bindingE,simulation_bindingE,positions[0:natoms_bindingE])
        Ut = 0;#total energy, not used for now
        print("binding energy="+str(be))
        f.write("%i %f %f %f %f\n" % (simulation.currentStep, temperature.value_in_unit(kelvin),lmbd, be.value_in_unit(kilocalorie_per_mole), Ut))	

f.close()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
print "Updating positions and velocities"
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")

