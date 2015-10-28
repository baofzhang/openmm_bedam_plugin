#ifndef OPENMM_REFERENCEBEDAMKERNELS_H_
#define OPENMM_REFERENCEBEDAMKERNELS_H_

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

#include "BEDAMKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/kernels.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "ReferenceStochasticDynamicsBEDAM.h"

using namespace OpenMM;

//class ReferenceStochasticDynamicsBEDAM;

namespace BEDAMPlugin {
//namespace OpenMM {

   /**          
   * This kernel is invoked by LangevinIntegratorBEDAM to take one time step.
   */
class ReferenceIntegrateLangevinStepBEDAMKernel : public IntegrateLangevinStepBEDAMKernel {
public:
  ReferenceIntegrateLangevinStepBEDAMKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateLangevinStepBEDAMKernel(name, platform),
  data(data), dynamics(0) {
  }
  ~ReferenceIntegrateLangevinStepBEDAMKernel();

   
   /**
   * Initialize the kernel, setting up the particle masses.
   * @param system     the System this kernel will be applied to
   * @param integrator the LangevinIntegrator this kernel will be used for
   */
   void initialize(const System& system, const LangevinIntegratorBEDAM& integrator);
   /**
    * Execute the kernel.
    * @param context    the context in which to execute this kernel       
    * @param integrator the LangevinIntegrator this kernel is being used for
    */
   void execute(ContextImpl& context, const LangevinIntegratorBEDAM& integrator);
    /**
     * Compute the kinetic energy.                               
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for    
     */
   double computeKineticEnergy(ContextImpl& context, const LangevinIntegratorBEDAM& integrator);
private:
  ReferencePlatform::PlatformData& data;
  ReferenceStochasticDynamicsBEDAM* dynamics;
  std::vector<RealOpenMM> masses;
  double prevTemp, prevFriction, prevStepSize;
  std::vector<RealVec> restraintForces; 
  };


} // namespace BEDAMPlugin

#endif /*OPENMM_REFERENCEBEDAMKERNELS_H_*/
