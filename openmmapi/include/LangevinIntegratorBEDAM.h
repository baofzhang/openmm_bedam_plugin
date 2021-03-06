#ifndef OPENMM_LANGEVININTEGRATORBEDAM_H_
#define OPENMM_LANGEVININTEGRATORBEDAM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include <vector>
#include <iostream>
#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExport.h"

using namespace OpenMM;
using namespace std;

namespace BEDAMPlugin {

//namespace OpenMM {
/**
 * This is an Integrator which simulates a System using Langevin dynamics.
 */

class OPENMM_EXPORT LangevinIntegratorBEDAM : public OpenMM::Integrator { 


public:
    /**
     * Create a LangevinIntegrator.
     * 
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     * @param ligId
     * @param lamdaId
     * @param atom1   the number of first atom for restraint force
     * @param atom2   the number of second atom for restraint force
     * @param kf      restraint force constant
     * @param r0      restraint force distance
     
     */
    
    LangevinIntegratorBEDAM(double temperature, double frictionCoeff, double stepSize,int ligId, double lamdaId, double kf, double r0);



    //add a CM atom for the ligand, atom1 is the index of an atom of the ligand
    void addAtom1Number(int atom1) {
      Atom1.push_back(atom1);
	}
    //return the id'th CM atom of the ligand or -1 if out of bounds
    int getAtom1Number(int id) const {
      if(id >= 0 && id < Atom1.size()) {
	return Atom1[id];
      }else{
	return -1;
      }
    }


    //add a CM atom for the receptor, atom2 is the index of an atom in the receptor
    void addAtom2Number(int atom2) {
      Atom2.push_back(atom2);
	}
    //return the id'th CM atom of the ligand or -1 if out of bounds
    int getAtom2Number(int id) const {
      if(id >= 0 && id < Atom2.size()) {
	return Atom2[id];
      }else{
	return -1;
      }
    }



    double getKf() const {
    return Kf;
  	}
    void setKf(double kf) {
      Kf = kf ; 
	}
    double getR0() const {
	return R0;
	}
    void setR0(double r0) {
	R0 = r0;
	}
    //get the ligId                                                                                                                   
    int getLigandId() const {
      return ligId;
    }
    //set the ligId                                                                                                                   
    void setLigandId(int ligandId){
      ligId = ligandId;
    }
    //get the lambdaId                                                                                                                 
    double getLamdaId() const {
      return lamdaid;
    }
    //set the lambdaId                                                                                                                 
    void setLamdaId(double lambdaid){
      lamdaid = lambdaid;
    }
    
    /**
     * Get the temperature of the heat bath (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getTemperature() const {
        return temperature;
    }
    /**
     * Set the temperature of the heat bath (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setTemperature(double temp) {
        temperature = temp;
    }
    /**
     * Get the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @return the friction coefficient, measured in 1/ps
     */
    double getFriction() const {
        return friction;
    }
    /**
     * Set the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @param coeff    the friction coefficient, measured in 1/ps
     */
    void setFriction(double coeff) {
        friction = coeff;
    }
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of random forces will be different.  On
     * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
     * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
     * results on successive runs, even if those runs were initialized identically.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();

private:
    double temperature, friction;
    int randomNumberSeed;
    Kernel kernel;
    int ligId;                                                                                                         
    double lamdaid; 
    vector<int> Atom1;
    vector<int> Atom2;
    double Kf;
    double R0;
    
};

} // namespace BEDAMPlugin

#endif /*OPENMM_LANGEVININTEGRATORBEDAM_H_*/ 
