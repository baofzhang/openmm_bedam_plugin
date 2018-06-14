/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/opencl/OpenCLKernels.h"
#include "openmm/opencl/OpenCLArray.h"
#include "OpenCLBEDAMKernels.h"
#include "OpenCLBEDAMKernelSources.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLExpressionUtilities.h"
#include "openmm/opencl/OpenCLIntegrationUtilities.h"
#include "openmm/opencl/OpenCLNonbondedUtilities.h"
#include <algorithm>
#include <cmath>
#include <set>

using namespace BEDAMPlugin;
using namespace OpenMM;
using namespace std;

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
  if (cl.getUseMixedPrecision())
    kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
  else
    kernel.setArg<void*>(index, NULL);
}


const float BOLTZMANN = 1.380658e-23f; // (J/K)
const float AVOGADRO = 6.0221367e23f;
const float RGAS = BOLTZMANN*AVOGADRO; // (J/(mol K))
const float BOLTZ = RGAS/1000;         // (kJ/(mol K))






OpenCLIntegrateLangevinStepBEDAMKernel::~OpenCLIntegrateLangevinStepBEDAMKernel() {
    if (params != NULL)
        delete params;
    
}

void OpenCLIntegrateLangevinStepBEDAMKernel::initialize(const System& system, const LangevinIntegratorBEDAM& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLBEDAMKernelSources::langevin, defines, "");

    kernel1 = cl::Kernel(program, "bedamForce");
    kernel1a = cl::Kernel(program, "CMForceCalcKernel");
    kernel2 = cl::Kernel(program, "integrateLangevinPart1");
    kernel3 = cl::Kernel(program, "integrateLangevinPart2");
    kernel4 = cl::Kernel(program, "copyDataToSecondPart");
    
    params = new OpenCLArray(cl, 3, cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(cl_double) : sizeof(cl_float), "langevinParams");
     
    prevStepSize = -1.0;

    //CM ligand atoms
    //vector<cl_int> atom1;
    int at, nat;
    
    at = integrator.getAtom1Number(0);
    nat = 0;
    while(at >= 0){
      nat += 1;
      atom1.push_back(at);
      at = integrator.getAtom1Number(nat);
    }

    //CM receptor atoms
    //vector<cl_int> atom2;
    int ligId = integrator.getLigandId();
    at = integrator.getAtom2Number(0);
    nat = 0;
    while(at >= 0){
      nat += 1;
      atom2.push_back(at+ligId);
      at = integrator.getAtom2Number(nat);
    }
  
    indexes1 = new OpenCLArray(cl,atom1.size(), sizeof(cl_int), "indexes1");
    indexes1->upload(atom1);

    indexes2 = new OpenCLArray(cl,atom2.size(), sizeof(cl_int), "indexes2");
    indexes2->upload(atom2);
    
}

void OpenCLIntegrateLangevinStepBEDAMKernel::execute(ContextImpl& context, const LangevinIntegratorBEDAM& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    
    int ligId = integrator.getLigandId();
    double lambdaId = integrator.getLamdaId();
    double kf = integrator.getKf();
    double r0 = integrator.getR0();
    
    
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
	
	
	kernel1.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
	kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
	
	kernel1a.setArg<cl::Buffer>(0,indexes1->getDeviceBuffer());
	kernel1a.setArg<cl::Buffer>(1,indexes2->getDeviceBuffer());
	kernel1a.setArg<cl::Buffer>(2,cl.getPosq().getDeviceBuffer());
	kernel1a.setArg<cl::Buffer>(3,cl.getForce().getDeviceBuffer());
	
	kernel2.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, params->getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
	
        kernel3.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel3, 1);
        kernel3.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel3.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel3.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());

	kernel4.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
	kernel4.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());


    }



    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize/tau);
        double fscale = (1-vscale)*tau;
        double noisescale = sqrt(2*kT/tau)*sqrt(0.5*(1-vscale*vscale)*tau);
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            vector<cl_double> p(params->getSize());
            p[0] = vscale;
            p[1] = fscale;
            p[2] = noisescale;
            params->upload(p);
            mm_double2 ss = mm_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            vector<cl_float> p(params->getSize());
            p[0] = (cl_float) vscale;
            p[1] = (cl_float) fscale;
            p[2] = (cl_float) noisescale;
            params->upload(p);
            mm_float2 ss = mm_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }


    
    //call the two atoms restraint force kernel1
    kernel1.setArg<cl_float>(2,lambdaId);
    cl.executeKernel(kernel1, numAtoms);

    kernel1a.setArg<cl_int>(4,atom1.size());
    kernel1a.setArg<cl_int>(5,atom2.size());
    kernel1a.setArg<cl_float>(6,kf);
    kernel1a.setArg<cl_float>(7,r0);
    cl.executeKernel(kernel1a, numAtoms);

    
    // Call the first integration kernel.

    kernel2.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel2, numAtoms); 

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel3, numAtoms); 

    
    // call the third coordinate manipulator kernel
    kernel4.setArg<cl_int>(2,ligId);
    cl.executeKernel(kernel4, numAtoms);
    
    
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateLangevinStepBEDAMKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegratorBEDAM& integrator) {
  return cl.getIntegrationUtilities().computeKineticEnergy(0.0*integrator.getStepSize());//ignore time shift
}
