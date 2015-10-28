from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, time, shutil
from dmsreader import *
from datetime import datetime
from BEDAMplugin import *

lambdaList = [1.0]
print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = 'binding_hct_lambda_one.dat'
f = open(binding_file, 'w')
for lambdaId in lambdaList:
    print("lambda = %f" %(lambdaId))
    testDes = DesmondDMSFile('bcd_recpt.dms','benzene_lig.dms',BEDAM=True) 
    system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.0*nanometer, OPLS = True,implicitSolvent=HCT)
    integrator = LangevinIntegratorBEDAM(300*kelvin, 1/picosecond, 0.002*picoseconds,12,lambdaId,109,154,1225.0,0.5)
    platform = Platform.getPlatformByName('OpenCL')
    simulation = Simulation(testDes.topology, system, integrator,platform)
    simulation.context.setPositions(testDes.positions)
    simulation.context.setVelocities(testDes.velocities)
    testDes.get_orig_force_parameters(system)
    context=simulation.context
    
    state = simulation.context.getState(getEnergy = True,getForces=True)
    print "Using platform %s" % simulation.context.getPlatform().getName()
    print(state.getPotentialEnergy())
    
    simulation.minimizeEnergy()
    stepId = 100
    loopStep = 20
    
    print "binding energy calculation starting"
    for i in range(loopStep):
        simulation.step(stepId)
        if simulation.currentStep%stepId==0:
            print "calculating binding energy"
            testDes.binding_energy_calculation(system,context,None)
            state=context.getState(getEnergy=True)
            U1=state.getPotentialEnergy()
            print(state.getPotentialEnergy())
            testDes.set_orig_force_parameters(system,context)
            testDes.binding_energy_calculation(system,context,True)
            state=context.getState(getEnergy=True)
            U0=state.getPotentialEnergy()
            print(state.getPotentialEnergy())
	    print("binding energy="+str(U1-U0))
	    if i>0:	
	        f.write("%f %f %f %f\n" % (0.0, 300.0,lambdaId, U1._value - U0._value))	
	    testDes.set_orig_force_parameters(system,context)
                            
f.close()
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
