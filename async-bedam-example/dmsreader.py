"""
desmonddmsfile.py: Load Desmond dms files

Portions copyright (c) 2013 Stanford University and the Authors
Authors: Robert McGibbon
Contributors: Emilio Gallicchio, Baofeng Zhang, Tony Zhao

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import os
import math


from simtk import openmm as mm
from simtk.openmm.app import forcefield as ff
from simtk.openmm.app import Element, Topology, PDBFile
from simtk.openmm.app.element import hydrogen
from simtk.unit import (nanometer, angstrom, dalton, radian,
                        kilocalorie_per_mole, kilojoule_per_mole,
                        degree, elementary_charge, femtosecond)

from simtk.openmm.app.amberprmtopfile import HCT
from simtk.openmm.app.internal.customgbforces import GBSAHCTForce
import simtk.unit as units
from simtk.openmm import CustomGBForce, Discrete1DFunction
import sqlite3
import numpy as np
#import scipy.integrate
import copy
import time
from abc import ABCMeta, abstractmethod
from simtk import openmm


class DesmondDMSFile(object):
    """DesmondDMSFile parses a Desmond DMS (desmond molecular system) and
    constructs a topology and (optionally) an OpenMM System from it
    """

    def __init__(self, ligfile, rcptfile, BEDAM=None):
        """Load a DMS file

        Parameters:
         - file (string) the name of the file to load
        """

        # sqlite3 is included in the standard lib, but at python
        # compile time, you can disable support (I think), so it's
        # not *guarenteed* to be available. Doing the import here
        # means we only raise an ImportError if people try to use
        # this class, so the module can be safely imported
        import sqlite3
	
        #add three lists for original force parameters
        self.nb_orig_para=[]
        self.cnb_orig_para=[]
        self.cgb_orig_para=[]
        self.implicitSolvent=None
        self.restraint_orig_para=[]
        self.temperature = 300.0 * units.kelvin

        self.rcptfile=rcptfile
        self.ligfile=ligfile
        
        if not os.path.exists(str(self.ligfile)):
            raise IOError("No such file or dirctory: '%'" % str(self.ligfile))
        if not  os.path.exists(str(self.rcptfile)):
            raise IOError("No such file or directory: '%s'" % str(self.rcptfile))

        self._open_lig = False
        self._conn_lig = sqlite3.connect(self.ligfile)
        self._open_lig = True

        self._open_rcpt = False
        self._conn_rcpt = sqlite3.connect(self.rcptfile)
        self._open_rcpt = True
        
	#to save the atom range in ligand and receptor
	self.ligand_atoms=[]
	self.receptor_atoms=[]

        #save the atom range in the imaginary ligand and receptor.
        self.image_ligand_atoms=[]
        self.image_receptor_atoms=[]
	
        self.bedam=BEDAM 
        if self.bedam is not None:
            self.bedamId=0  #flag for BEDAM simulation and initial atom id for imaginary system
            self.file1Id=0
            self.bedamX=1000.0  #distance between the original system and the imaginary system along X axis in angstroms#
            self.bedam1X=2000.0 #ligand
            
        self.tables_lig = None
        self.tables_rcpt = None
        self._readSchemas()

        if len(self.tables_lig) == 0:
            raise IOError('ligand DMS file was not loaded sucessfully. No tables found')
        if 'nbtype' not in self.tables_lig['particle']:
            raise ValueError('No nonbonded parameters associated with ligand '
                             'DMS file. You can add a forcefield with the '
                             'viparr command line tool distributed with desmond')

        if len(self.tables_rcpt) == 0:
            raise IOError('receptor DMS file was not loaded sucessfully. No tables found')
        if 'nbtype' not in self.tables_rcpt['particle']:
            raise ValueError('No nonbonded parameters associated with receptor '
                             'DMS file. You can add a forcefield with the '
                             'viparr command line tool distributed with desmond')
        
        # Build the topology
        self.topology, self.positions, self.velocities = self._createTopology()
        self._topologyAtoms = list(self.topology.atoms())
        self._atomBonds = [{} for x in range(len(self._topologyAtoms))]
        self._angleConstraints = [{} for x in range(len(self._topologyAtoms))]

    def getPositions(self):
        """Get the positions of each atom in the system
        """
        return self.positions

    def getVelocities(self):
        """Get the positions of each atom in the system
        """
        return self.velocities

    def getTopology(self):
        """Get the topology of the system
        """
        return self.topology

    def _addtoTopology(self, top, atoms, positions, velocities, sqlconnection, id_offset=0, pos_offset=mm.Vec3(0,0,0)):
        """Build the topology of the system
        """

        
        #TODO: how to manage more than one global_cell?
        boxVectors = []
        for x, y, z in sqlconnection.execute('SELECT x, y, z FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z))
        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions*angstrom)

        lastChain = None
        lastResId = None
        if id_offset == 0:
            c = top.addChain() #needed?
            
        q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz FROM particle"""
        for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz) in sqlconnection.execute(q):
            newChain = False
            if chain != lastChain:
                lastChain = chain
                c = top.addChain()
                newChain = True
            if resId != lastResId or newChain:
                lastResId = resId
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}

            if atomNumber == 0 and atomName.startswith('Vrt'):
                elem = None
            else:
                elem = Element.getByAtomicNumber(atomNumber)
	    	
            if atomName in atomReplacements:
                atomName = atomReplacements[atomName]

            atoms[atomId+id_offset] = top.addAtom(atomName, elem, r)
            positions.append(mm.Vec3(x, y, z)+pos_offset)

            velocities.append(mm.Vec3(vx, vy, vz))

        for p0, p1 in sqlconnection.execute('SELECT p0, p1 FROM bond'):
            top.addBond(atoms[p0+id_offset], atoms[p1+id_offset])

        return top, positions, velocities

    def _createTopology(self):
        top = Topology()
        atoms = {}
        positions = []
        velocities = []
        self._addtoTopology(top, atoms, positions, velocities, self._conn_lig)
	self.ligand_atoms=range(0,len(atoms))
        self._addtoTopology(top, atoms, positions, velocities, self._conn_rcpt, len(atoms))
	self.receptor_atoms=range(len(self.ligand_atoms),len(atoms))
        if self.bedam is not None:
            self.image_atoms = range(len(atoms),2*len(atoms))
            self.total_atoms = range(0,2*len(atoms))
            self.bedamId = len(atoms) 
            lig_displacement = mm.Vec3(self.bedam1X,0,0)
            rcpt_displacement = mm.Vec3(self.bedamX, 0, 0)
            self._addtoTopology(top, atoms, positions, velocities, self._conn_lig, len(atoms), pos_offset=lig_displacement)
            self._addtoTopology(top, atoms, positions, velocities, self._conn_rcpt, len(atoms), pos_offset=rcpt_displacement)
        positions = positions*angstrom
        velocities = velocities*angstrom/femtosecond
        return top, positions, velocities
        
    def setPositions(self, positions):
        """Update atomic positions in attached DMS file
        """
        q = """UPDATE particle SET x = ?1, y = ?2, z = ?3 WHERE id == ?4"""

        jat = 0
	for j in self.ligand_atoms:
	    (x, y , z) = positions[j].value_in_unit(angstrom)
            self._conn_lig.execute(q, (x,y,z,jat))
            jat += 1
	self._conn_lig.commit()

        iat = 0
        for i in self.receptor_atoms:
            (x, y , z) = positions[i].value_in_unit(angstrom)
            self._conn_rcpt.execute(q, (x,y,z,iat))
            iat += 1
        self._conn_rcpt.commit()

        return iat

    def setVelocities(self, velocities):
        """Update atomic positions in attached DMS file
        """
        q = """UPDATE particle SET vx = ?1, vy = ?2, vz = ?3 WHERE id == ?4"""
        jat = 0
	for j in self.ligand_atoms:
            (vx, vy , vz) = velocities[j].value_in_unit(angstrom/femtosecond)
            self._conn_lig.execute(q, (vx,vy,vz,jat))
            jat += 1
        self._conn_lig.commit()

        iat = 0
        for i in self.receptor_atoms:
            (vx, vy , vz) = velocities[i].value_in_unit(angstrom/femtosecond)
            self._conn_rcpt.execute(q, (vx,vy,vz,iat))
            iat += 1
        self._conn_rcpt.commit()
        
        return iat


    def _get_gb_paramsW(self, sqlconnection, gb_p):
        """
        get charge, radius, screened_radius from hct table in the .dms file to calculate
        --radiusN=radius-0.009  # Offset radius, the radius needs to be converted to nanometer unit
        --screenN=screened_radius*radiusN #Scaled offset radius 
        """
	length_conv = units.angstrom.conversion_factor_to(units.nanometer)

        q="""SELECT charge,radius,screened_radius from hct"""
        try:
            iat = 0
            for charge,radius,screened_radius in sqlconnection.execute(q):
            
                chargeN = charge
                radiusN = radius
                screenN = screened_radius
                radiusN *=length_conv
                radiusN -=0.009
                screenN *=radiusN
                gb_p.append((chargeN,radiusN,screenN))
                iat += 1
            return iat
        except:
            return None

    def _get_gb_params(self):
        gb_p = []
        self._get_gb_paramsW(self._conn_lig, gb_p)
        self._get_gb_paramsW(self._conn_rcpt, gb_p)
        if self.bedam is not None:
            self._get_gb_paramsW(self._conn_lig, gb_p)
            self._get_gb_paramsW(self._conn_rcpt, gb_p)        
        return gb_p

    def _get_agbnp2_paramsW(self, sqlconnection, gb_p):
        """
        get charge, radius, etc. from AGBNP2 table in the .dms file and computes AGBNP2 parameters for each atom
        """
	length_conv = units.angstrom.conversion_factor_to(units.nanometer)
        en_conv = units.kilocalorie_per_mole.conversion_factor_to(units.kilojoule_per_mole)
        gamma_conv = en_conv/(length_conv*length_conv)
        alpha_conv = en_conv*length_conv*length_conv*length_conv;

        # gather AGBNP2 parameters from agbnp2 table
        #   atomic number and charge from particle table matching atom id
        #   sigma and epsilon from nonbonded_param table matching hbtype from particle
        q="""SELECT anum,charge,radius,igamma,ialpha,idelta,sgamma,salpha,sdelta,hbtype,hbw,sigma,epsilon from particle INNER JOIN agbnp2 ON particle.id==agbnp2.id INNER JOIN nonbonded_param ON particle.nbtype==nonbonded_param.id ORDER BY particle.id"""
        
        try:
            iat = 0
            for anum,charge,radius,igamma,ialpha,idelta,sgamma,salpha,sdelta,hbtype,hbw,sigma,epsilon in sqlconnection.execute(q):

                if anum == 1:
                    ishydrogenN = 1
                else:
                    ishydrogenN = 0

                radiusN = length_conv*radius
                chargeN = charge
                gammaN = gamma_conv*(igamma+0.*sgamma)
                alphaN = alpha_conv*(ialpha+salpha)
                # delta parameter is ignored
                hbtypeN = hbtype
                hbwN = en_conv * hbw
                
                gb_p.append([radiusN,chargeN,gammaN,alphaN,hbtype,hbwN,ishydrogenN])
                iat += 1
            return iat
        except:
            print("Warning: unable to retrieve AGBNP parameters")
            return None

    def _get_agbnp2_params(self):
        gb_p = []
        self._get_agbnp2_paramsW(self._conn_lig, gb_p)
        self._get_agbnp2_paramsW(self._conn_rcpt, gb_p)
        if self.bedam is not None:
            self._get_agbnp2_paramsW(self._conn_lig, gb_p)
            self._get_agbnp2_paramsW(self._conn_rcpt, gb_p)        
        return gb_p
    
    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*nanometer,
                     ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None, OPLS=False, implicitSolvent=None, AGBNPVersion=1):
        """Construct an OpenMM System representing the topology described by this dms file

        Parameters:
         - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
           NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
         - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
         - ewaldErrorTolerance (float=0.0005) The error tolerance to use if nonbondedMethod is Ewald or PME.
         - removeCMMotion (boolean=True) If true, a CMMotionRemover will be added to the System
         - hydrogenMass (mass=None) The mass to use for hydrogen atoms bound to heavy atoms.  Any mass added to a hydrogen is
           subtracted from the heavy atom to keep their total mass the same.
         - OPLS (boolean=False) If True, force field parameters are interpreted as OPLS parameters; OPLS variants of 
           torsional and non-bonded forces are constructed.
	 - implicitSolvent (object=None) if not None, the implicit solvent model to use, the only allowed value is HCT
        
        """
        self._checkForUnsupportedTerms()
        sys = mm.System()
        self.OPLS = OPLS
        self.implicitSolvent=implicitSolvent

        # Buld the box dimensions
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')
        
        #TODO: how to manipulate the boxSize for duplicate system, right now, only consider implicit solvent model //02/26/15

        # Create all of the particles
        for mass in self._conn_lig.execute('SELECT mass from particle'):
            sys.addParticle(mass[0]*dalton)
        for mass in self._conn_rcpt.execute('SELECT mass from particle'):
            sys.addParticle(mass[0]*dalton)
        if self.bedam is not None:
            for mass in self._conn_lig.execute('SELECT mass from particle'):
                sys.addParticle(mass[0]*dalton)
            for mass in self._conn_rcpt.execute('SELECT mass from particle'):
                sys.addParticle(mass[0]*dalton)

        # Add all of the forces
        self._addBondsToSystem(sys)
        self._addAnglesToSystem(sys)
        self._addConstraintsToSystem(sys)
        self._addPeriodicTorsionsToSystem(sys, OPLS)
        self._addImproperHarmonicTorsionsToSystem(sys)
        #self._addCMAPToSystem(sys)
        #self._addVirtualSitesToSystem(sys)
        self._addPositionalHarmonicRestraints(sys)
        nb, cnb = self._addNonbondedForceToSystem(sys, OPLS)
        
        # Finish configuring the NonbondedForce.
        methodMap = {ff.NoCutoff:mm.NonbondedForce.NoCutoff,
                     ff.CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     ff.Ewald:mm.NonbondedForce.Ewald,
                     ff.PME:mm.NonbondedForce.PME}
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setUseDispersionCorrection(False)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)
        if cnb is not None:
            cnb.setNonbondedMethod(methodMap[nonbondedMethod])
            cnb.setCutoffDistance(nonbondedCutoff)
            cnb.setUseSwitchingFunction(False)
            cnb.setUseLongRangeCorrection(False);
            
	#add implicit solvent model
	if implicitSolvent is not None:

            #with implicit solvent turn off native reaction field
            #However note that this does not affect the shifted Coulomb potential of the Nonbonded force
            #(it affects the only the energy, not the forces and equation of motion)
            nb.setReactionFieldDielectric(1.0)

            if implicitSolvent is HCT:
                gb_parms = self._get_gb_params()
                if gb_parms:
                    print('Adding HCT GB force ...')
                    gb = GBSAHCTForce(SA='ACE', cutoff=nonbondedCutoff._value)
                    for i in range(len(gb_parms)):	      
                        gb.addParticle(list(gb_parms[i]))#edit 3.13.15

         	# Set cutoff method //2.10.15 for setting up cutoff options in implicit solvent model
                if nonbondedMethod is ff.NoCutoff:
                    gb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
                elif nonbondedMethod is ff.CutoffNonPeriodic:
                    gb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                    gb.setCutoffDistance(nonbondedCutoff)
                elif nonbondedMethod is ff.CutoffPeriodic:
                    gb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
                    gb.setCutoffDistance(nonbondedCutoff)
                else:
                    raise ValueError('Illegal nonbonded method for use with GBSA')
                #gb.setForceGroup(self.GB_FORCE_GROUP)
                gb.finalize()
                sys.addForce(gb)
                self.hct_force = gb;
                self.hct_gp_parms = gb_parms
                
            if implicitSolvent is 'GVolSA':
                #implemented as AGBNP version 0
                implicitSolvent = 'AGBNP'
                AGBNPVersion = 0
                print('Using GVolSA')

            if implicitSolvent is 'AGBNP':
                #load AGBNP plugin if available
                try:
                    from AGBNPplugin import AGBNPForce
                    AGBNPEnabled = True
                except ImportError:
                    AGBNPEnabled = False
                #sets up AGBNP
                if AGBNPEnabled:
                    gb_parms = self._get_agbnp2_params()
                    if gb_parms:
                        gb = AGBNPForce()
                        gb.setNonbondedMethod(methodMap[nonbondedMethod])
                        gb.setCutoffDistance(nonbondedCutoff)
                        gb.setVersion(AGBNPVersion)
                        print('Using AGBNP force version %d' % AGBNPVersion)
                        # add particles
                        for i in range(len(gb_parms)):
                            [radiusN,chargeN,gammaN,alphaN,hbtype,hbwN,ishydrogenN] = gb_parms[i]
                            h_flag = ishydrogenN > 0
                            gb.addParticle(radiusN, gammaN, alphaN, chargeN, h_flag)
                            #print("Adding", radiusN, gammaN, h_flag)
                        sys.addForce(gb)
                        self.agbnp_force = gb
                        self.agbnp_gb_parms = gb_parms
                else:
                    print('Warning: AGBNP is not supported in this version')
            
        # Adjust masses.
        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == hydrogen and atom1.element not in (hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)
		
        # Add a CMMotionRemover.
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        return sys

    def _addBonds(self, sys, bondforce, sqlconnection, id_offset=0):
        q = """SELECT p0, p1, r0, fc, constrained
        FROM stretch_harm_term INNER JOIN stretch_harm_param
        ON stretch_harm_term.param=stretch_harm_param.id"""
        for p0, p1, r0, fc, constrained in sqlconnection.execute(q):
            if constrained:
                sys.addConstraint(p0+id_offset, p1+id_offset, r0*angstrom)
            else:
                # Desmond writes the harmonic bond force without 1/2
                # so we need to to double the force constant
                bondforce.addBond(p0+id_offset, p1+id_offset, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

            # Record information that will be needed for constraining angles.
            self._atomBonds[p0+id_offset][p1+id_offset] = r0*angstrom
            self._atomBonds[p1+id_offset][p0+id_offset] = r0*angstrom
    
    def _addBondsToSystem(self, sys):
        """Create the harmonic bonds
        """
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)
        self._addBonds(sys, bonds, self._conn_lig)
        self._addBonds(sys, bonds, self._conn_rcpt, self.receptor_atoms[0])
        if self.bedam is not None:
            self._addBonds(sys, bonds, self._conn_lig, self.bedamId)
            self._addBonds(sys, bonds, self._conn_rcpt, self.bedamId + self.receptor_atoms[0])
        return bonds

    def _addAngles(self, sys, angleforce, sqlconnection, id_offset=0):
        degToRad = math.pi/180
        q = """SELECT p0, p1, p2, theta0, fc, constrained
        FROM angle_harm_term INNER JOIN angle_harm_param
        ON angle_harm_term.param=angle_harm_param.id"""
        for p0, p1, p2, theta0, fc, constrained in sqlconnection.execute(q):
            if constrained:
                l1 = self._atomBonds[p1+id_offset][p0+id_offset]
                l2 = self._atomBonds[p1+id_offset][p2+id_offset]
                length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                sys.addConstraint(p0, p2, length)
                self._angleConstraints[p1+id_offset][p0+id_offset] = p2
                self._angleConstraints[p1+id_offset][p2+id_offset] = p0
            else:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angleforce.addAngle(p0+id_offset, p1+id_offset, p2+id_offset, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)

    
    def _addAnglesToSystem(self, sys):
        """Create the harmonic angles
        """
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        self._addAngles(sys, angles, self._conn_lig)
        self._addAngles(sys, angles, self._conn_rcpt, self.receptor_atoms[0])
        if self.bedam is not None:
            self._addAngles(sys, angles, self._conn_lig, self.bedamId)
            self._addAngles(sys, angles, self._conn_rcpt, self.bedamId + self.receptor_atoms[0])
        return angles

    
    def _addConstraints(self, sys, sqlconnection, tables, id_offset=0):
        for term_table in [n for n in list(tables.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
            param_table = term_table.replace('term', 'param')
            q = """SELECT p0, p1, r1
            FROM %(term)s INNER JOIN %(param)s
            ON %(term)s.param=%(param)s.id""" % \
                {'term': term_table, 'param': param_table}
            for p0, p1, r1 in sqlconnection.execute(q):
                if not p1+id_offset in self._atomBonds[p0+id_offset]:
                    sys.addConstraint(p0+id_offset, p1+id_offset, r1*angstrom)
                    self._atomBonds[p0+id_offset][p1+id_offset] = r1*angstrom
                    self._atomBonds[p1+id_offset][p0+id_offset] = r1*angstrom

        if 'constraint_hoh_term' in tables:
            degToRad = math.pi/180
            q = """SELECT p0, p1, p2, r1, r2, theta
            FROM constraint_hoh_term INNER JOIN constraint_hoh_param
            ON constraint_hoh_term.param=constraint_hoh_param.id"""
            for p0, p1, p2, r1, r2, theta in sqlconnection.execute(q):
                # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2
                # wihth O-H1 and O-H2 distances r1 and r2
                if not (self._angleConstraints[p0+id_offset].get(p1+id_offset, None) == p2+id_offset):
                    length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                    sys.addConstraint(p1+id_offset, p2+id_offset, length)
        
    
    def _addConstraintsToSystem(self, sys):
        """Add constraints to system. Normally these should already be
        added by the bonds table, but we want to make sure that there's
        no extra information in the constraints table that we're not
        including in the system"""
        self._addConstraints(sys, self._conn_lig, self.tables_lig)
        self._addConstraints(sys, self._conn_rcpt, self.tables_rcpt, self.receptor_atoms[0])
        if self.bedam is not None:
            self._addConstraints(sys, self._conn_lig, self.tables_lig, self.bedamId)
            self._addConstraints(sys, self._conn_rcpt, self.tables_rcpt, self.bedamId + self.receptor_atoms[0])

    def _addPeriodicTorsions(self, sys, pforce, OPLS, sqlconnection, id_offset=0):
        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id"""
        for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in sqlconnection.execute(q):
            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue
                if OPLS:
                    pforce.addTorsion(p0+id_offset, p1+id_offset, p2+id_offset, p3+id_offset,
                                        [order, phi0*degree, fc*kilocalorie_per_mole])
                else:
                    pforce.addTorsion(p0+id_offset, p1+id_offset, p2+id_offset, p3+id_offset,
                                        order, phi0*degree, fc*kilocalorie_per_mole)
            
    def _addPeriodicTorsionsToSystem(self, sys, OPLS):
        """Create the torsion terms
        """
        if OPLS:
            periodic = mm.CustomTorsionForce('f * cos(n * theta - phi0)')
            periodic.addPerTorsionParameter('n')
            periodic.addPerTorsionParameter('phi0')
            periodic.addPerTorsionParameter('f')
        else:
            periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)
        self._addPeriodicTorsions(sys, periodic, OPLS, self._conn_lig)
        self._addPeriodicTorsions(sys, periodic, OPLS, self._conn_rcpt, self.receptor_atoms[0])
        if self.bedam is not None:
            self._addPeriodicTorsions(sys, periodic, OPLS, self._conn_lig,  self.bedamId)
            self._addPeriodicTorsions(sys, periodic, OPLS, self._conn_rcpt, self.bedamId + self.receptor_atoms[0])


    def _addImproperHarmonicTorsions(self, harmonictorsion, sqlconnection, id_offset=0):
        q = """SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id"""
        for p0, p1, p2, p3, phi0, fc in sqlconnection.execute(q):
            harmonictorsion.addTorsion(p0+id_offset, p1+id_offset, p2+id_offset, p3+id_offset,
                                       [phi0*degree, fc*kilocalorie_per_mole])
            
    def _addImproperHarmonicTorsionsToSystem(self, sys):
        """Create the improper harmonic torsion terms
        """
        if ( not self._hasTable(self.tables_lig, 'improper_harm_term') ) and ( not self._hasTable(self.tables_rcpt, 'improper_harm_term') ):
            return

        harmonicTorsion = mm.CustomTorsionForce('k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)

        if self._hasTable(self.tables_lig, 'improper_harm_term'):
            self._addImproperHarmonicTorsions(harmonicTorsion, self._conn_lig)
        if self._hasTable(self.tables_rcpt, 'improper_harm_term'):
            self._addImproperHarmonicTorsions(harmonicTorsion, self._conn_rctp, self.receptor_atoms[0])
        if self.bedam is not None:
            if self._hasTable(self.tables_lig, 'improper_harm_term'):
                self._addImproperHarmonicTorsions(harmonicTorsion, self._conn_lig, self.bedamId)
            if self._hasTable(self.tables_rcpt, 'improper_harm_term'):
                self._addImproperHarmonicTorsions(harmonicTorsion, self._conn_rctp, self.bedamId + self.receptor_atoms[0])

    def _addNonbondedParticles(self, nbforce, cnbforce, sqlconnection, offset = 0):
        nb = nbforce
        cnb = cnbforce
        nb_query = "SELECT charge, sigma, epsilon FROM particle INNER JOIN nonbonded_param ON particle.nbtype=nonbonded_param.id"
        n_atoms_added = 0
        for charge, sigma, epsilon in sqlconnection.execute(nb_query):
            cnb.addParticle([charge, sigma*angstrom, epsilon*kilocalorie_per_mole])
            n_atoms_added += 1
        for charge, sigma, epsilon in sqlconnection.execute(nb_query):
            epsilon = 0
            charge = 0
            nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)
        for p0, p1 in sqlconnection.execute('SELECT p0, p1 FROM exclusion'):
            nb.addException(p0+offset, p1+offset, 0.0, 1.0, 0.0)
            cnb.addExclusion(p0+offset, p1+offset)
        q = "SELECT p0, p1, aij, bij, qij FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param ON pair_12_6_es_term.param=pair_12_6_es_param.id;"
        for p0, p1, a_ij, b_ij, q_ij in sqlconnection.execute(q):
            a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2
            if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                new_epsilon = 0
                new_sigma = 1
            else:
                new_epsilon =  b_ij**2/(4*a_ij)
                new_sigma = (a_ij / b_ij)**(1.0/6.0)                
            nb.addException(p0+offset, p1+offset, q_ij, new_sigma, new_epsilon, True)

        n_total = sqlconnection.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
        n_in_exclusions = sqlconnection.execute("""SELECT COUNT(*)
        FROM exclusion INNER JOIN pair_12_6_es_term
        ON (exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
        OR (exclusion.p1==pair_12_6_es_term.p0 AND exclusion.p0==pair_12_6_es_term.p1)""").fetchone()
        if not n_total == n_in_exclusions:
            raise NotImplementedError('All pair_12_6_es_term must have a corresponding exclusion')

        return n_atoms_added
        
    def _addNonbondedForceToSystem(self, sys, OPLS):
        """Create the nonbonded force
        """
        cnb = None #custom general force

        #need a dummy nonbonded potential to handle exceptions below
        nb = mm.NonbondedForce()
        sys.addForce(nb)
        self.nb_force = nb;

        #custom force with possibly OPLS combination rules 
        ONE_4PI_EPS0 = 138.935456
        umax = 4186 #maximum pair potential energy in kJ/mol (~1000 kcal/mol)
        #step function ensures that negative energies are not capped 
        pair_potential = "step(ueff)*ueff+(1-step(ueff))*u ; ueff = umax*tanh(u/umax) ; u = ONE_4PI_EPS0*charge1*charge2/r + 4.0*epsilon12*sc6*(sc6 - 1); sc6=(sigma12/r)^6"
        if OPLS:
            combination_rules = "sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2)"
        else:
            combination_rules = "sigma12=0.5*(sigma1+sigma2); epsilon12=sqrt(epsilon1*epsilon2)"
        parameters = "ONE_4PI_EPS0 = %f ; umax = %f" % (ONE_4PI_EPS0, umax)
        expression = pair_potential + ";" + combination_rules + ";" + parameters
        cnb = mm.CustomNonbondedForce(expression)
        cnb.addPerParticleParameter("charge")
        cnb.addPerParticleParameter("sigma")
        cnb.addPerParticleParameter("epsilon")
        sys.addForce(cnb)
        self.nb_opls_force = cnb

        self._addNonbondedParticles(nb, cnb, self._conn_lig)
        self._addNonbondedParticles(nb, cnb, self._conn_rcpt, offset = self.receptor_atoms[0])
        #duplicate everything
        if self.bedam is not None:
            self._addNonbondedParticles(nb, cnb, self._conn_lig,  offset = self.bedamId)
            self._addNonbondedParticles(nb, cnb, self._conn_rcpt, offset = self.bedamId + self.receptor_atoms[0])

        return nb, cnb

    def _addPositionalHarmonicRestraintsW(self, cforce, sqlconnection, id_offset=0, pos_offset = mm.Vec3(0,0,0)):
        nrestraint = 0
        q = """SELECT p0, x0, y0, z0, fcx, fcy, fcz FROM posre_harm_term INNER JOIN posre_harm_param ON posre_harm_term.param=posre_harm_param.id"""
        for p0, x0, y0, z0, fcx, fcy, fcz in sqlconnection.execute(q):
            nrestraint += 1
            cforce.addParticle(p0+id_offset,[ (x0+pos_offset[0])*angstrom, (y0+pos_offset[1])*angstrom, (z0++pos_offset[2])*angstrom,
                                   0.5 * fcx * kilocalorie_per_mole/angstrom**2,
                                   0.5 * fcy * kilocalorie_per_mole/angstrom**2,
                                   0.5 * fcz * kilocalorie_per_mole/angstrom**2])
        return nrestraint
    
    def _addPositionalHarmonicRestraints(self, sys):
        if (not self._hasTable(self.tables_lig, 'posre_harm_term')) and (not self._hasTable(self.tables_rcpt, 'posre_harm_term')):
            return
        if self._hasTable(self.tables_lig, 'posre_harm_term') and not self._hasTable(self.tables_lig, 'posre_harm_param'):
            raise IOError('DMS file lacks posre_harm_param table even though posre_harm_term tableis present.')
            return
        if self._hasTable(self.tables_rcpt, 'posre_harm_term') and not self._hasTable(self.tables_rcpt, 'posre_harm_param'):
            raise IOError('DMS file lacks posre_harm_param table even though posre_harm_term tableis present.')
            return
        force = mm.CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        force.addPerParticleParameter("hkx")
        force.addPerParticleParameter("hky")
        force.addPerParticleParameter("hkz")
        nrestraint = 0;
        if self._hasTable(self.tables_lig, 'posre_harm_term'):
            nrestraint += self._addPositionalHarmonicRestraintsW(force, self._conn_lig)
        if self._hasTable(self.tables_rcpt, 'posre_harm_term'):
            nrestraint += self._addPositionalHarmonicRestraintsW(force, self._conn_rcpt, self.receptor_atoms[0])
        if self.bedam is not None:
            lig_displacement = mm.Vec3(self.bedam1X,0,0)
            rcpt_displacement = mm.Vec3(self.bedamX, 0, 0)
            if self._hasTable(self.tables_lig, 'posre_harm_term'):
                nrestraint += self._addPositionalHarmonicRestraintsW(force, self.bedamId, lig_displacement)
            if self._hasTable(self.tables_rcpt, 'posre_harm_term'):
                nrestraint += self._addPositionalHarmonicRestraintsW(force, self._conn_rcpt, self.bedamId + self.receptor_atoms[0], rcpt_displacement)
                
        if nrestraint > 0:
            print("Adding positional restraining force ...")
            sys.addForce(force)
            
    def _hasTable(self, tables, table_name):
        """Does our DMS file contain this table?
        """
        return table_name in tables

    def _readSchemasW(self, sqlconnection):
        """Read the schemas of each of the tables in the dms file, populating
        the `tables` instance attribute
        """
        tables = {}
        for table in sqlconnection.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in sqlconnection.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        return tables
    
    def _readSchemas(self):
        """Read the schemas of each of the tables in the dms files, populating
        the `tables` instance attribute
        """
        self.tables_lig = self._readSchemasW(self._conn_lig)
        self.tables_rcpt = self._readSchemasW(self._conn_rcpt)
        
    def _checkForUnsupportedTermsW(self, tables, sqlconnection):
        """Check the file for forcefield terms that are not currenty supported,
        raising a NotImplementedError
        """
        flat_bottom_potential_terms = ['stretch_fbhw_term', 'angle_fbhw_term',
                                       'improper_fbhw_term', 'posre_fbhw_term']
        if any((t in tables) for t in flat_bottom_potential_terms):
            raise NotImplementedError('Flat bottom potential terms '
                                      'are not implemeneted')

        nbinfo = dict(list(zip(tables['nonbonded_info'],
                          sqlconnection.execute('SELECT * FROM nonbonded_info').fetchone())))

        if nbinfo['vdw_funct'] != 'vdw_12_6':
            raise NotImplementedError('Only Leonard-Jones van der Waals '
                                      'interactions are currently supported')
        if nbinfo['vdw_rule'] != 'arithmetic/geometric':
            """raise NotImplementedError('Only Lorentz-Berthelot nonbonded '
                                      'combining rules are currently supported')"""
            nbinfo['vdw_rule'] = 'arithmetic/geometric'

        """if 'nonbonded_combined_param' in self._tables:
            raise NotImplementedError('nonbonded_combined_param interactions '
                                      'are not currently supported')"""

        if 'alchemical_particle' in tables:
            raise NotImplementedError('Alchemical particles are not supported')
        if 'alchemical_stretch_harm' in tables:
            raise NotImplementedError('Alchemical bonds are not supported')

        if 'polar_term' in tables:
            if sqlconnection.execute("SELECT COUNT(*) FROM polar_term").fetchone()[0] != 0:
                raise NotImplementedError('Drude particles are not currently supported')

    def _checkForUnsupportedTerms(self):
        self._checkForUnsupportedTermsW(self.tables_lig, self._conn_lig)
        self._checkForUnsupportedTermsW(self.tables_rcpt, self._conn_rcpt)
            

    def displaceLigand(self, positions):
        new_positions = list(positions) #a copy
        delta = -200.0 * nanometer;
        for i in self.ligand_atoms:
            (x,y,z) = positions[i]
            new_positions[i] = mm.Vec3(x._value+delta._value,y._value,z._value) * nanometer
        return new_positions
            
    def close(self):
        """Close the SQL connection
        """
        if self._open_lig:
            self._conn_lig.close()
        if self._open_rcpt:
            self._conn_rcpt.close()
            
    def __del__(self):
        self.close()



