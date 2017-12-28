
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil
import multiprocessing as mp
from dmsreader import *
from datetime import datetime
from BEDAMplugin import *

def binding_energy(dms,simulation,positions):
    context = simulation.context
    #energy of bound system
    simulation.context.setPositions(positions)
    state=context.getState(getEnergy=True)
    U1=state.getPotentialEnergy()
    #print("U1=%f" % U1.value_in_unit(kilocalorie_per_mole))
    #displace ligand
    new_positions = dms.displaceLigand(positions)
    simulation.context.setPositions(new_positions)
    state=context.getState(getEnergy=True)
    #energy of unbound energy
    U0=state.getPotentialEnergy()
    #print("U0=%f" % U0.value_in_unit(kilocalorie_per_mole))
    #simulation.context.setPositions(positions)
    return U1-U0

def binding_energy_mp(dms_bindingE,simulation_bindingE,positions, temperature, lmbd, q):
    be = binding_energy(dms_bindingE,simulation_bindingE,positions)
    Ut = 0;#total energy, not used for now
    #print("binding energy="+str(be))
    output = "%i %f %12.6e %12.6e %12.6e\n" % (simulation.currentStep, temperature.value_in_unit(kelvin),lmbd, be.value_in_unit(kilocalorie_per_mole), Ut)
    q.put(output)
    
print("Started at: " + str(time.asctime()))
start=datetime.now()

binding_file = '@jobname@_@n@.out'

lmbd = @lambda@
print("lambda = %12.6e" %(lmbd))

rcptfile_input  = '@jobname@_rcpt_@nm1@.dms'
ligfile_input   = '@jobname@_lig_@nm1@.dms'
rcptfile_output = '@jobname@_rcpt_tmp.dms'
ligfile_output  = '@jobname@_lig_tmp.dms'
rcptfile_result = '@jobname@_rcpt_@n@.dms'
ligfile_result  = '@jobname@_lig_@n@.dms'


shutil.copyfile(rcptfile_input, rcptfile_output)
shutil.copyfile(ligfile_input, ligfile_output)

testDes = DesmondDMSFile(ligfile_output, rcptfile_output, BEDAM=True) 
system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.2*nanometer, OPLS = True, implicitSolvent='AGBNP', AGBNPVersion=1)

temperature = @temperature@ * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
natoms_ligand = 43
lig_atom_restr = 12   #index of ligand atom for Vsite (267CA -- central C of ligand)
rcpt_atom_restr = 1535   #index of rcpt atom center of Vsite#267CA -- central C of ligand
kf = 1225.0 #force constant for Vsite in (kj/mol)/nm^2
r0 = 0.5 #radius of Vsite sphere

integrator = LangevinIntegratorBEDAM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond,natoms_ligand,lmbd,lig_atom_restr, rcpt_atom_restr,kf, r0)

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

#simulation context to compute binding energies asymchronously, without image and without a cutoff
dms_bindingE = DesmondDMSFile(ligfile_output, rcptfile_output)
system_bindingE = dms_bindingE.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent='AGBNP', AGBNPVersion=1)
platform_bindingE = Platform.getPlatformByName('Reference')
platform_properties = {}
integrator_bindingE = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
simulation_bindingE = Simulation(dms_bindingE.topology, system_bindingE, integrator_bindingE, platform_bindingE, platform_properties)
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

processes = []
output_queues = []
for i in range(loopStep):
    start_sim = datetime.now()
    simulation.step(stepId)
    end_sim = datetime.now()
    elapsed=end_sim - start_sim
    print("elapsed simulation time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
    if simulation.currentStep%stepId==0:
        positions = simulation.context.getState(getPositions=True).getPositions()
        q = mp.Queue()
        p = mp.Process(target=binding_energy_mp, args=(dms_bindingE,simulation_bindingE,positions[0:natoms_bindingE], temperature, lmbd,q))
        processes.append(p)
        output_queues.append(q)
        p.start()
        
f = open(binding_file, 'w')        
for i in range(0,len(processes)):
    p = processes[i]
    q = output_queues[i]
    if p.is_alive():
        p.join()
    f.write(q.get())
        
f.close()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
print "Updating positions and velocities"
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

shutil.copyfile(rcptfile_output, rcptfile_result)
shutil.copyfile(ligfile_output, ligfile_result)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")

