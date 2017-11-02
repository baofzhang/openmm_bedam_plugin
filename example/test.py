from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, time, shutil
from dmsreader import *
from datetime import datetime
from BEDAMplugin import *

lambdaList = [0.001]
print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = 'binding_hct_lambda_one.dat'
f = open(binding_file, 'w')
for lmbd in lambdaList:
    print("lambda = %f" %(lmbd))
    shutil.copyfile('bcd_recpt.dms','bcd_recpt-out.dms')
    shutil.copyfile('benzene_lig.dms','benzene_lig-out.dms')
    testDes = DesmondDMSFile('bcd_recpt-out.dms','benzene_lig-out.dms',BEDAM='yes') 
    system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=2*nanometer, OPLS = True, implicitSolvent='AGBNP')

    temperature = 300.0 * kelvin
    frictionCoeff = 0.5 / picosecond
    MDstepsize = 0.002 * picosecond
    natoms_ligand = 12
    rcpt_atom_restr = 109
    lig_atom_restr = 154
    kf = 1225.0
    r0 = 0.5
    

    integrator = LangevinIntegratorBEDAM(temperature / kelvin, frictionCoeff/(1/picosecond), MDstepsize / picosecond, natoms_ligand, lmbd, rcpt_atom_restr, lig_atom_restr, kf, r0)
    platform = Platform.getPlatformByName('OpenCL')
    simulation = Simulation(testDes.topology, system, integrator,platform)
    simulation.context.setPositions(testDes.positions)
    simulation.context.setVelocities(testDes.velocities)
    testDes.get_orig_force_parameters(system)
    context=simulation.context
    
    state = simulation.context.getState(getEnergy = True,getForces=True)
    print "Using platform %s" % simulation.context.getPlatform().getName()

    print("Potential energy of doubled system: %f" % state.getPotentialEnergy()._value)
    
    testDes.binding_energy_calculation(system,context,None)
    state=context.getState(getEnergy=True)
    U1=state.getPotentialEnergy()
    #print(state.getPotentialEnergy())
    testDes.set_orig_force_parameters(system,context)
    testDes.binding_energy_calculation(system,context,True)
    state=context.getState(getEnergy=True)
    U0=state.getPotentialEnergy()
    #print(state.getPotentialEnergy())
    print("binding energy="+str(U1-U0))
    print("%f %f %f %f\n" % (0.0, 300.0,lmbd, U1._value - U0._value))	
    testDes.set_orig_force_parameters(system,context)
    state=context.getState(getEnergy=True)
    print("Potential energy of doubled system (should be the same as above): %f" % state.getPotentialEnergy()._value)
    
    simulation.reporters.append(StateDataReporter(stdout, 500, step=True, temperature=True))
    
    print("Minimizing system ...")
    simulation.minimizeEnergy()
    stepId = 500
    loopStep = 10
    
    print "MD with binding energy calculation starting"
    for i in range(loopStep):
        simulation.step(stepId)
        if simulation.currentStep%stepId==0:
            testDes.binding_energy_calculation(system,context,None)
            state=context.getState(getEnergy=True)
            U1=state.getPotentialEnergy()
            #print(state.getPotentialEnergy())
            testDes.set_orig_force_parameters(system,context)
            testDes.binding_energy_calculation(system,context,True)
            state=context.getState(getEnergy=True)
            U0=state.getPotentialEnergy()
            #print(state.getPotentialEnergy())
	    print("binding energy="+str(U1-U0))
	    testDes.set_orig_force_parameters(system,context)


positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
#print "Updating positions and velocities"
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()#
            
f.close()
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
