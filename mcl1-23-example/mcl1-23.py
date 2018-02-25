
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
    #print(U1)
    #displace ligand
    new_positions = dms.displaceLigand(positions)
    simulation.context.setPositions(new_positions)
    state=context.getState(getEnergy=True)
    #energy of unbound energy
    U0=state.getPotentialEnergy()
    #print(U0)
    #simulation.context.setPositions(positions)
    return U1-U0


print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = '@jobname@_@n@.out'
f = open(binding_file, 'w')

lmbd = @lambda@
print("lambda = %f" %(lmbd))
shutil.copyfile('@jobname@_rcpt_@nm1@.dms','@jobname@_rcpt_@n@.dms') 
shutil.copyfile('@jobname@_lig_@nm1@.dms','@jobname@_lig_@n@.dms') 

testDes = DesmondDMSFile('@jobname@_rcpt_@n@.dms','@jobname@_lig_@n@.dms',BEDAM=True) 

#system = testDes.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent='AGBNP')
system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.2*nanometer, OPLS = True,implicitSolvent='AGBNP', AGBNPVersion=1)

temperature = @temperature@ * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
natoms_ligand = 43
rcpt_atom_restr = 1535  #index of rcpt atom center of Vsite#267CA -- central C of ligand
lig_atom_restr = 2435  #index of ligand atom for Vsite (267CA -- central C of ligand)
kf = 1255.8 #force constant for Vsite in (kj/mol)/nm^2
r0 = 0.7 #radius of Vsite sphere
    
integrator = LangevinIntegratorBEDAM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond,natoms_ligand,lmbd,rcpt_atom_restr,lig_atom_restr,kf, r0)

platform = Platform.getPlatformByName('@platform@')

properties = {}

if '@platform@'=='OpenCL':
    #expected "platformid:deviceid" or empty
    device = "@pn@"
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

#simulation context to compute binding energies without image and without a cutoff (large cutoff)
dms_bindingE = DesmondDMSFile('@jobname@_rcpt_@n@.dms','@jobname@_lig_@n@.dms')
system_bindingE = dms_bindingE.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=20*nanometer, OPLS = True,implicitSolvent='AGBNP', AGBNPVersion=1)
integrator_bindingE = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
simulation_bindingE = Simulation(dms_bindingE.topology, system_bindingE, integrator_bindingE,platform, properties)
natoms_bindingE = len(dms_bindingE.positions)
print("Number of receptor atoms:"+str(natoms_bindingE - natoms_ligand))
print("Number of ligand atoms:"+str(natoms_ligand))
print("Number of atoms in complex:"+str(natoms_bindingE))
print("Number of atoms in BEDAM system:"+str(2*natoms_bindingE))

#compute initial binding energy (optional)
positions = simulation.context.getState(getPositions=True).getPositions()
be = binding_energy(dms_bindingE,simulation_bindingE,positions[0:natoms_bindingE])
print("initial binding energy="+str(be))

stepId = @stepgap@
totalSteps = @totalsteps@
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
        f.write("%i %f %f %f %f\n" % (simulation.currentStep, @temperature@,lmbd, be.value_in_unit(kilocalorie_per_mole), Ut))	

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
