%module BEDAMplugin

%import(module="simtk.openmm") "OpenMMSwigHeaders.i"


/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "LangevinIntegratorBEDAM.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"


%}


/*
 * The code below strips all units before the wrapper
 * functions are called. This code also converts numpy
 * arrays to lists.
*/

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}


/* strip the units off of all input arguments */
%pythonprepend %{
try:
    args=mm.stripUnits(args)
except UnboundLocalError:
    pass
%}


/*
 * Add units to function inputs.

%pythonappend BEDAMPlugin::LangevinIntegratorBEDAM::LangevinIntegratorBEDAM(double temperature, double frictionCoeff, double stepSize, int ligId, double lamdaId, int atom1, int atom2, double kf, double r0) %{
    val[0] = unit.Quantity(val[0], unit.kelvin)
    val[1] = unit.Quantity(val[1], 1/unit.picosecond)
    val[2] = unit.Quantity(val[2], 1/unit.picosecond)
    val[4] = unit.Quantity(val[4], None)
    val[7] = unit.Quantity(val[7], unit.kilojoule_per_mole/(unit.nanometer*unit.nanometer))
    val[8] = unit.Quantity(val[8], unit.nanometer)
%}
*/



namespace BEDAMPlugin{

class LangevinIntegratorBEDAM : public OpenMM::Integrator {
public:
/*
   %apply double INPUT {double temperature};
   %apply double INPUT {double frictionCoeff};
   %apply double INPUT {double stepSize};
   %apply int INPUT {int ligId};
   %apply double INPUT {double lamdaId};
   %apply int INPUT {int atom1};
   %apply int INPUT {int atom2};
   %apply double INPUT {double kf};
   %apply double INPUT {double r0};
*/
   LangevinIntegratorBEDAM(double temperature, double frictionCoeff, double stepSize, int ligId, double lamdaId, int atom1, int atom2, double kf, double r0) ;
   
   int getAtom1Number() const ;
   void setAtom1Number(int atom1) ;
   int getAtom2Number() const ;
   void setAtom2Number(int atom2) ;
   double getKf() const ;
   void setKf(double kf) ;
   double getR0() const ;
   void setR0(double r0) ;
   int getLigandId() const ;
   void setLigandId(int ligandId) ;
   double getLamdaId() const ;
   void setLamdaId(double lambdaid) ;
   double getTemperature() const ;
   void setTemperature(double temp) ;
   double getFriction() const ;
   void setFriction(double coeff) ;
   int getRandomNumberSeed() const ;
   void setRandomNumberSeed(int seed) ;
   virtual void step(int steps) ;
};




}

