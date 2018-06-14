/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceBEDAMKernels.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "ReferenceStochasticDynamicsBEDAM.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace BEDAMPlugin;
using namespace OpenMM;
using namespace std;

static int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

static RealOpenMM** allocateRealArray(int length, int width) {
    RealOpenMM** array = new RealOpenMM*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new RealOpenMM[width];
    return array;
}

static void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static void disposeRealArray(RealOpenMM** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    //vector<RealVec>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    //silently ignore time shift (forces are not reliable here because they have not been alchemically hybridized)
    vector<RealVec> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
      shiftedVel[i] = velData[i]; //+forceData[i]*(timeShift/masses[i]);
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}


ReferenceIntegrateLangevinStepBEDAMKernel::~ReferenceIntegrateLangevinStepBEDAMKernel() {
  if (dynamics)
    delete dynamics;
}

void ReferenceIntegrateLangevinStepBEDAMKernel::initialize(const System& system, const LangevinIntegratorBEDAM& integrator) {
int numParticles = system.getNumParticles();
masses.resize(numParticles);
for (int i = 0; i < numParticles; ++i)
  masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateLangevinStepBEDAMKernel::execute(ContextImpl& context, const LangevinIntegratorBEDAM& integrator){

 double temperature = integrator.getTemperature();
 double friction = integrator.getFriction();
 double stepSize = integrator.getStepSize();
 int ligId = integrator.getLigandId();
 double lambdaId = integrator.getLamdaId();
 vector<RealVec>& posData = extractPositions(context);
 vector<RealVec>& velData = extractVelocities(context);
 vector<RealVec>& forceData = extractForces(context);
 int numParticles = context.getSystem().getNumParticles();
 int halfN = numParticles/2;

 //ligand atoms are stored first, followed by receptor atoms
 vector<int> atom1;
 int at, nat;
 at = integrator.getAtom1Number(0);
 nat = 0;
 while(at >= 0){
   nat += 1;
   atom1.push_back(at);
   at = integrator.getAtom1Number(nat);
 }
 
 vector<int> atom2;
 at = integrator.getAtom2Number(0);
 nat = 0;
 while(at >= 0){
   nat += 1;
   atom2.push_back(at+ligId);
   at = integrator.getAtom2Number(nat);
 }
 

 //impose flat-bottom CM restraints (Vsite)
 double kf = integrator.getKf();
 double r0 = integrator.getR0();    

 RealVec rf1 = RealVec(0.,0.,0.);
 RealVec rf2 = RealVec(0.,0.,0.);

 if(atom1.size() > 0 && atom2.size() > 0){//do this only if CM's are defined
   RealVec cm1 = RealVec(0,0,0);
   for (int i = 0; i < atom1.size(); i++) {
     int index = atom1[i];
     cm1 += posData[index]/atom1.size();
    }
   RealVec cm2 = RealVec(0,0,0);
   for (int i = 0; i < atom2.size(); i++) {
     int index = atom2[i];
     cm2 += posData[index]/atom2.size();
    }

   RealVec dist = cm2 - cm1;
   double dr = sqrt(dist.dot(dist));
   double vn1 = 1./atom1.size();
   double vn2 = 1./atom2.size();
   
   if(dr>r0) 
    {
      rf1 = dist *  vn1*kf*(dr-r0)/dr;
      rf2 = dist * -vn2*kf*(dr-r0)/dr;
    }
   
  }
   
 for (int i = 0 ; i < halfN; ++i ){
   forceData[i][0] = lambdaId*forceData[i][0]+(1.0-lambdaId)*forceData[i+halfN][0] ;
   forceData[i][1] = lambdaId*forceData[i][1]+(1.0-lambdaId)*forceData[i+halfN][1] ;
   forceData[i][2] = lambdaId*forceData[i][2]+(1.0-lambdaId)*forceData[i+halfN][2] ; 
 }
   
 for(int i=0;i<atom1.size();i++){
   int index = atom1[i];
   forceData[index][0] += rf1[0];
   forceData[index][1] += rf1[1];
   forceData[index][2] += rf1[2];
 }
 for(int i=0;i<atom2.size();i++){
   int index = atom2[i];
   forceData[index][0] += rf2[0];
   forceData[index][1] += rf2[1];
   forceData[index][2] += rf2[2];
 }
                                                                                                       
 if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
   // Recreate the computation objects with the new parameters.                                                            

   if (dynamics)
     delete dynamics;
   RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
   dynamics = new ReferenceStochasticDynamicsBEDAM(
						   context.getSystem().getNumParticles(),
						   static_cast<RealOpenMM>(stepSize),
						   static_cast<RealOpenMM>(tau),
						   static_cast<RealOpenMM>(temperature) );
   dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
   prevTemp = temperature;
   prevFriction = friction;
   prevStepSize = stepSize;
 }

 dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance(),ligId);                                                                                                                       

  data.time += stepSize;
 data.stepCount++;
}


double ReferenceIntegrateLangevinStepBEDAMKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegratorBEDAM& integrator) {
  return computeShiftedKineticEnergy(context, masses, 0.0*integrator.getStepSize());//ignore time shift
}

