
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <cstring>
#include <sstream>

//#include "openmm/reference/SimTKOpenMMCommon.h"
//#include "openmm/reference/SimTKOpenMMLog.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"
#include "ReferenceStochasticDynamicsBEDAM.h"
#include "openmm/reference/ReferenceVirtualSites.h"
#include <iostream>
#include <cstdio>

using std::vector;
using OpenMM::RealVec;
using namespace BEDAMPlugin;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param tau            viscosity(?)
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamicsBEDAM::ReferenceStochasticDynamicsBEDAM( int numberOfAtoms,
                                                          RealOpenMM deltaT, RealOpenMM tau,
                                                          RealOpenMM temperature ) : 
           ReferenceDynamics( numberOfAtoms, deltaT, temperature ), _tau( tau ) {

   // ---------------------------------------------------------------------------------------

   static const char* methodName      = "\nReferenceStochasticDynamics::ReferenceStochasticDynamics";
   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   // ensure tau is not zero -- if it is print warning message

#ifdef NOTNOW
   if( _tau == zero ){

      std::stringstream message;
      message << methodName;
      message << " input tau value=" << tau << " is invalid -- setting to 1.";
      SimTKOpenMMLog::printError( message );

      _tau = one;
     
   }
#endif
   xPrime.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceStochasticDynamicsBEDAM::~ReferenceStochasticDynamicsBEDAM( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceStochasticDynamics::~ReferenceStochasticDynamics";

   // ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get tau

   @return tau

   --------------------------------------------------------------------------------------- */

RealOpenMM ReferenceStochasticDynamicsBEDAM::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nReferenceStochasticDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   First SD update; based on code in update.c do_update_sd() Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param inverseMasses       inverse atom masses
   @param xPrime              xPrime

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamicsBEDAM::updatePart1( int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                              vector<RealVec>& velocities,
                                              vector<RealVec>& forces, vector<RealOpenMM>& inverseMasses,
                                              vector<RealVec>& xPrime ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\nReferenceStochasticDynamics::updatePart1";

   // ---------------------------------------------------------------------------------------

   // perform first update

   RealOpenMM tau = getTau();
   const RealOpenMM vscale = EXP(-getDeltaT()/tau);
   const RealOpenMM fscale = (1-vscale)*tau;
   const RealOpenMM kT = BOLTZ*getTemperature();
   const RealOpenMM noisescale = SQRT(2*kT/tau)*SQRT(0.5*(1-vscale*vscale)*tau);

   //   for (int ii = 0; ii < numberOfAtoms; ii++) {
   //  std::cout << "V1:  " << ii << " " << " " << velocities[ii][0] << " " << velocities[ii][1] << " " << velocities[ii][2] << std::endl;
   // }

   
   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (inverseMasses[ii] != 0.0) {
           RealOpenMM sqrtInvMass = SQRT(inverseMasses[ii]);
           for (int jj = 0; jj < 3; jj++) {
               velocities[ii][jj]  = vscale*velocities[ii][jj] + fscale*inverseMasses[ii]*forces[ii][jj] + noisescale*sqrtInvMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
           } 
       }
   }

   //for (int ii = 0; ii < numberOfAtoms; ii++) {
   //  std::cout << "V2:  " << ii << " " << " " << velocities[ii][0] << " " << velocities[ii][1] << " " << velocities[ii][2] << std::endl;
   //}


}

/**---------------------------------------------------------------------------------------

   Second update; based on code in update.c do_update_sd() w/ bFirstHalf = false in Gromacs 3.1.4

   @param numberOfAtoms       number of atoms
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamicsBEDAM::updatePart2( int numberOfAtoms, vector<RealVec>& atomCoordinates,
                                              vector<RealVec>& velocities,
                                              vector<RealVec>& forces, vector<RealOpenMM>& inverseMasses,
                                              vector<RealVec>& xPrime ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\nReferenceStochasticDynamics::updatePart2";

   // ---------------------------------------------------------------------------------------

   // perform second update

   for (int ii = 0; ii < numberOfAtoms; ii++) {
       if (inverseMasses[ii] != 0.0)
           for (int jj = 0; jj < 3; jj++)
               xPrime[ii][jj] = atomCoordinates[ii][jj]+getDeltaT()*velocities[ii][jj];
   }
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceStochasticDynamicsBEDAM::update(const OpenMM::System& system, vector<RealVec>& atomCoordinates,
					 vector<RealVec>& velocities, vector<RealVec>& forces, vector<RealOpenMM>& masses, RealOpenMM tolerance, int ligId) {

   // ---------------------------------------------------------------------------------------

  static const char* methodName      = "\nReferenceStochasticDynamicsBEDAM::update"; 

   static const RealOpenMM zero       =  0.0;
   static const RealOpenMM one        =  1.0;

   // ---------------------------------------------------------------------------------------

   // first-time-through initialization

   int numberOfAtoms = system.getNumParticles();
   if( getTimeStep() == 0 ){
      // invert masses

      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         if (masses[ii] == zero)
             inverseMasses[ii] = zero;
         else
             inverseMasses[ii] = one/masses[ii];
      }
   }

   // 1st update

   updatePart1( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime );

   // 2nd update

   updatePart2( numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime );

   ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
   if (referenceConstraintAlgorithm)
      referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

   // copy xPrime -> atomCoordinates

   RealOpenMM invStepSize = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i)
       if (masses[i] != zero)
           for (int j = 0; j < 3; ++j) {
               velocities[i][j] = invStepSize*(xPrime[i][j]-atomCoordinates[i][j]);
               atomCoordinates[i][j] = xPrime[i][j];
           }

   //ligand atoms are stored first, followed by receptor atoms
   int halfN = numberOfAtoms/2;
   for (int j = halfN; j<halfN+ligId; ++j){//ligand image

     atomCoordinates[j][0] = atomCoordinates[j-halfN][0]+200.0;
     atomCoordinates[j][1] = atomCoordinates[j-halfN][1];
     atomCoordinates[j][2] = atomCoordinates[j-halfN][2];

     velocities[j][0] = velocities[j-halfN][0];
     velocities[j][1] = velocities[j-halfN][1];
     velocities[j][2] = velocities[j-halfN][2];
   }
   for (int j = halfN+ligId; j<numberOfAtoms; ++j ) {//receptor image

     atomCoordinates[j][0] = atomCoordinates[j-halfN][0]+100.0;
     atomCoordinates[j][1] = atomCoordinates[j-halfN][1];
     atomCoordinates[j][2] = atomCoordinates[j-halfN][2];
     
     velocities[j][0] = velocities[j-halfN][0];
     velocities[j][1] = velocities[j-halfN][1];
     velocities[j][2] = velocities[j-halfN][2];
   }
   
   ReferenceVirtualSites::computePositions(system, atomCoordinates);
   incrementTimeStep();
}
