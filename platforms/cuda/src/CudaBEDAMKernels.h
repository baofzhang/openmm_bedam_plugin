#ifndef CUDA_BEDAM_KERNELS_H_
#define CUDA_BEDAM_KERNELS_H_

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

#include "BEDAMKernels.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaArray.h"

using namespace OpenMM;
namespace BEDAMPlugin {


/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CudaIntegrateLangevinStepBEDAMKernel : public IntegrateLangevinStepBEDAMKernel {
public:
    CudaIntegrateLangevinStepBEDAMKernel(std::string name, const OpenMM::Platform& platform, OpenMM::CudaContext& cu) : IntegrateLangevinStepBEDAMKernel(name, platform), cu(cu), params(NULL) {
    }
    ~CudaIntegrateLangevinStepBEDAMKernel();
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const LangevinIntegratorBEDAM& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    void execute(OpenMM::ContextImpl& context, const LangevinIntegratorBEDAM& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(OpenMM::ContextImpl& context, const LangevinIntegratorBEDAM& integrator);
private:
    OpenMM::CudaContext& cu;
    double prevTemp, prevFriction, prevStepSize;
    OpenMM::CudaArray* params;
    CUfunction kernel1, kernel2,kernel3,kernel4;
};

} // namespace BEDAMPlugin

#endif /*CUDA_BEDAM_KERNELS_H_*/

