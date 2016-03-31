from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, time, shutil
from dmsreaderwscwrestraint import *
from datetime import datetime
from BEDAMplugin import *
import numpy as np

lambdaList = [@lambda@]
print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = '@jobname@_@n@.out'
f = open(binding_file, 'w')

print("lambda = %f" %(@lambda@))
if @n@ > 1:
    testDes = DesmondDMSFile('@jobname@_rcpt_@nm1@.dms','@jobname@_lig_@nm1@.dms',LIGAND=True,BEDAM=True,CAREST=True,wscore=0.01) 
else:
    testDes = DesmondDMSFile('@jobname@_rcpt.dms','@jobname@_lig.dms',LIGAND=True,BEDAM=True,CAREST=True,wscore=0.01)
system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=2.0*nanometer, OPLS = True,implicitSolvent=@implicitsolvent@,SCORE=True)
integrator = LangevinIntegratorBEDAM(@temperature@*kelvin, 1/picosecond, 0.001*picoseconds,43,@lambda@,884,2443,1254,0.5) #249-CA-C17
platform = Platform.getPlatformByName('@platform@')
platform.setPropertyDefaultValue("OpenCLDeviceIndex", "@pn@")
simulation = Simulation(testDes.topology, system, integrator,platform)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)
testDes.get_orig_force_parameters(system)
context=simulation.context
state = simulation.context.getState(getEnergy = True,getForces=True)
print "Using platform %s" % simulation.context.getPlatform().getName()
print(state.getPotentialEnergy())
simulation.minimizeEnergy()
stepId = @stepgap@
totalSteps = @totalsteps@
loopStep = totalSteps/stepId

print "binding energy calculation starting"    
for i in range(loopStep):
    simulation.step(stepId)
    if simulation.currentStep%stepId==0:
        print "calculating binding energy"
        testDes.binding_energy_calculation(system,context,None)
        state=context.getState(getEnergy=True)
        pe1=state.getPotentialEnergy()
        print(state.getPotentialEnergy())
        testDes.set_orig_force_parameters(system,context)
        testDes.binding_energy_calculation(system,context,True)
        state=context.getState(getEnergy=True)
        pe0=state.getPotentialEnergy()
	ke0=state.getKineticEnergy()._value/4.18/2.0
        print(state.getPotentialEnergy())
	be = (pe1._value-pe0._value)/4.18
	if be>1000000.0 or np.isnan(be):
	    be = 1000000.0
	tpe = pe0._value/4.18 + @lambda@*be
	te = ke0 + tpe
	#if np.isnan(te):
	#    raise Exception("total energy for the system is NaN in step %i" %simulation.currentStep)
	print("binding energy="+str(be))
        f.write("%i %f %f %f %f\n" % (simulation.currentStep,@temperature@,@lambda@,be,te))	
        testDes.set_orig_force_parameters(system,context)
	
f.close()
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
print "Updating positions and velocities"
shutil.copyfile('@jobname@_rcpt.dms','@jobname@_rcpt_@n@.dms') 
shutil.copyfile('@jobname@_lig.dms','@jobname@_lig_@n@.dms') 
output = DesmondDMSFile('@jobname@_rcpt_@n@.dms','@jobname@_lig_@n@.dms',LIGAND=True,BEDAM=True,CAREST=True,wscore=0.01) 
output.setPositions(positions)
output.setVelocities(velocities)
output.close()
testDes.close()

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
