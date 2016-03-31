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
ONE_4PI_EPS0 = 138.935485 # OpenMM constant for Coulomb interactions (openmm/platforms/reference/include/SimTKOpenMMRealType.h) in OpenMM units


class DesmondDMSFile(object):
    """DesmondDMSFile parses a Desmond DMS (desmond molecular system) and
    constructs a topology and (optionally) an OpenMM System from it
    """
    GB = "BASIC"
    def __init__(self, rcptfile,ligfile,LIGAND=False,BEDAM=None,Image=None,alphaId=None,CAREST=None,wscore=0.0):
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
	

        #edit on 11.16.15
	#add a soft-core parameter to the nonbondedforce
	#self.softcore=None
	#edit end on 11.16.15

        #edit 3.11.15
        #add four lists for original force parameters for nb,cnb,colnb and cgb
        self.nb_orig_para=[]
        self.cnb_orig_para=[]
	self.cbnb_orig_para=[]
        self.colnb_orig_para=[]
        self.cgb_orig_para=[]
        
        #self.implicitSolvent=None
        #edit 3.17.15
        self.restraint_orig_para=[]
        #edit end 3.17.15
        #edit 3.17.15
        self.carest = CAREST
        self.wscore = wscore
        #edit 12.28.15
        self.calphas={}
        

        self.temperature = 300.0 * units.kelvin

        #edit end 3.17.15
        
        #edit end 3.11.15

        #edit on 4.6.15
        infiniteN = 1e-7
        #self.lambda_restraint = 1.0/lambda_id
        #edit end 4.6.15

        #added on 2.27.15

        self.rcptfile=rcptfile
        self.ligfile=ligfile
        self.ligand = LIGAND
        self.Image=Image
        self.score = None
	
	#edit 3.10.15
	#try to save the atom range in ligand and receptor
	#receptor atoms
	self.receptor_atoms=[]
	#ligand atoms
	self.ligand_atoms=[]
        #save the atom range in the imaginary ligand and receptor. #edit 3.25.15
        self.image_receptor_atoms=[]
        self.image_ligand_atoms=[]
        #edit 3.25.15
	#control the force calculation so that calculating the binding free energy
	#self.alpha_i=1.0
	

	#edit end3.10.15


	#edit 3.9.15
	#alphaId (alphaId=None) If not None, the force control coefficient fcId is set to zero
	#self.fcId = 1.0
	#self.alphaId = alphaId
	#if self.alphaId is not None:
	#    self.fcId = 0.0
        #edit end 3.9.15

	#added on 02/26/15
	#BEDAM (BEDAM=None) If not None, the BEDAM version is used.

        self.bedam=BEDAM 
        if self.bedam is not None:
            self.bedamId=0  #flag for BEDAM simulation and initial atom id for imaginary system #added on 02/26/15
            self.file1Id=0
            self.bedamX=100.0  #distance between the original system and the imaginary system along X axis #added on 02/26/15
            self.bedam1X=200.0
            
        self._open = False
        self._tables = None
        if not  os.path.exists(str(self.rcptfile)):
            raise IOError("No such file or directory: '%s'" % str(self.rcptfile))
        self._conn = sqlite3.connect(self.rcptfile)
        #added on 2.27.15
        #if self.ligfile:
        if self.ligand:             
            self._tables1 = None
            if not os.path.exists(str(self.ligfile)):
                raise IOError("No such file or dirctory: '%'" % str(self.ligfile))
            self._conn1 = sqlite3.connect(self.ligfile) 
        #added end on 2.27.15
        self._open = True
        self._readSchemas()

        if len(self._tables) == 0:
            raise IOError('receptor DMS file was not loaded sucessfully. No tables found')
        
        if 'nbtype' not in self._tables['particle']:
            raise ValueError('No nonbonded parameters associated with this '
                             'receptor DMS file. You can add a forcefield with the '
                             'viparr command line tool distributed with desmond')
        if self.ligand:
            if len(self._tables1) == 0:
                raise IOError('ligand DMS file was not loaded sucessfully. No tables found')
            #if 'nbtype' not in self._tables1['particle']:
            #    raise ValueError('No nonbonded parameters associated with this '
            #                 'ligand DMS file. You can add a forcefield with the '
            #                 'viparr command line tool distributed with desmond')
        #added on 3.1.15
        #if self.ligfile:
        #    if len(self._tables1) == 0:
        #        raise IOError('DMS file was not loaded sucessfully. No tables found')
        #    if 'nbtype' not in self._tables1['particle']:
        #        raise ValueError('No nonbonded parameters associated with this '
        #                         'DMS file. You can add a forcefield with the '
        #                         'viparr command line tool distributed with desmond')

        #added end on 3.1.15
        # build the provenance string
        provenance = []
        q = """SELECT id, user, timestamp, version, workdir, cmdline, executable FROM provenance"""
        try:
            for id, user, timestamp, version, workdir, cmdline, executable in self._conn.execute(q):
                for row in self._conn.execute('SELECT * FROM provenance'):
                    rowdict = dict(list(zip(self._tables['provenance'], row)))
                    provenance.append('%(id)d) %(timestamp)s: %(user)s\n  version: %(version)s\n  '
                                      'cmdline: %(cmdline)s\n  executable: %(executable)s\n' % rowdict)
                    self.provenance = ''.join(provenance)
        except:
            print "DesmondDMSFile: Warning, unable to retrieve provenance."

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

    def getProvenance(self):
        """Get the provenance string of this system
        """
        return self.provenance

    def _createTopology(self):
        """Build the topology of the system
        """
        top = Topology()
        positions = []
        #calphapositions= [] #edit 12/28/15
        #calphalist=[] #edit 12/28/15
        #calphas = {} #edit 12/28/15
        velocities = []

        boxVectors = []
        for x, y, z in self._conn.execute('SELECT x, y, z FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z))
        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions*angstrom)

        #TODO: how to manage more than one global_cell 3.1.15

        atoms = {}
        lastChain = None
        lastResId = None
        c = top.addChain()
        q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz,x0,y0,z0,grp_buffer FROM particle""" #edit on 12/28/15
        for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz,x0,y0,z0,grp_buffer) in self._conn.execute(q):
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

            atoms[atomId] = top.addAtom(atomName, elem, r)
            positions.append(mm.Vec3(x, y, z))
            velocities.append(mm.Vec3(vx, vy, vz))
            #add on 12.28.15
            if(grp_buffer==2):
                self.calphas[atomId]=mm.Vec3(x0,y0,z0)
            
        #print self.calphas

            
	#print len(atoms)
        for p0, p1 in self._conn.execute('SELECT p0, p1 FROM bond'):
            top.addBond(atoms[p0], atoms[p1])

	#collect the number of atoms in receptor
	self.receptor_atoms=range(0,len(atoms))  #edit on 3.10.15
	print len(self.receptor_atoms)

        #edited on 02/27/15
        
        if self.ligand:
            self.fileId = len(atoms) 
            
            #self.bedamX = 10.0   #TODO: need to figure out how to measure the longest length of the system
            lastChain = None
            lastResId = None
                        
            q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz FROM particle"""
            for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz) in self._conn1.execute(q):
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

                atoms[atomId+self.fileId] = top.addAtom(atomName, elem, r)
                #print atomId+self.bedamId
                positions.append(mm.Vec3(x, y, z))

                velocities.append(mm.Vec3(vx, vy, vz))
	
            for p0, p1 in self._conn1.execute('SELECT p0, p1 FROM bond'):
                top.addBond(atoms[p0+self.fileId], atoms[p1+self.fileId])

            
        #edit end on 02/27/15

	#collect the number of atoms in ligands
            self.ligand_atoms=range(len(self.receptor_atoms),len(atoms)) #edit 3.10.15
            print len(self.ligand_atoms)
        self.original_atoms = range(0,len(atoms))
        self.total_atoms = range(0,len(atoms))
        #print self.original_atoms
        if self.bedam is not None:
            self.image_atoms = range(len(atoms),2*len(atoms)) #edit 3.11.15
            #print self.image_atoms
            self.total_atoms = range(0,2*len(atoms))
            print len(self.total_atoms)

        #edited on 02/26/15
        
        if self.bedam is not None:
            self.bedamId = len(atoms) 
            
            lastChain = None
            lastResId = None
            
            q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz FROM particle"""
            for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz) in self._conn.execute(q):
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

                atoms[atomId+self.bedamId] = top.addAtom(atomName, elem, r)
                #print atomId+self.bedamId
                positions.append(mm.Vec3(x+self.bedamX, y, z))

                velocities.append(mm.Vec3(vx, vy, vz))
	
            for p0, p1 in self._conn.execute('SELECT p0, p1 FROM bond'):
                top.addBond(atoms[p0+self.bedamId], atoms[p1+self.bedamId])
                
            #save the id of imaginary receptor #edit 3.25.15
            self.image_receptor_atoms=range(self.bedamId, len(atoms))
            #edit end 3.25.15
            
            

            
        #edit end on 02/26/15
            #edited on 02/27/15
        
            if self.ligand:
                self.file1Id = len(atoms) 
            
                #self.bedamX = 10.0   #TODO: need to figure out how to measure the longest length of the system
                lastChain = None
                lastResId = None
                        
                q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz FROM particle"""
                for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz) in self._conn1.execute(q):
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

                    atoms[atomId+self.file1Id] = top.addAtom(atomName, elem, r)
                    #print atomId+self.bedamId
                    positions.append(mm.Vec3(x+self.bedam1X, y, z))

                    velocities.append(mm.Vec3(vx, vy, vz))
                
                for p0, p1 in self._conn1.execute('SELECT p0, p1 FROM bond'):
                    top.addBond(atoms[p0+self.file1Id], atoms[p1+self.file1Id])

                #save the id of imaginary ligand atoms #edit 3.25.15
                self.image_ligand_atoms=range(self.file1Id, len(atoms))
                #edit end 3.25.15
            
        #edit end on 02/27/15
        #print atoms

        positions = positions*angstrom
        velocities = velocities*angstrom/femtosecond
        
        print positions[0]

        return top, positions, velocities



    #TODO: need to revise setpositions, setvelocities. Ideally they should only update original system. //02/26/15



    def setPositions(self, positions):
        """Update atomic positions in attached DMS file
        """
        q = """UPDATE particle SET x = ?1, y = ?2, z = ?3 WHERE id == ?4"""
        iat = 0
	#edit 3.10.15
        for i in self.receptor_atoms:
            (x, y , z) = positions[i].value_in_unit(angstrom)
            self._conn.execute(q, (x,y,z,iat))
            iat += 1
        self._conn.commit()
        #edit 4.1.15
        
	if self.ligand:
            jat = 0
	    for j in self.ligand_atoms:
		(x, y , z) = positions[j].value_in_unit(angstrom)
            	self._conn1.execute(q, (x,y,z,jat))
            	jat += 1
	    self._conn1.commit()
	
	#edit end 3.10.15
        return iat

    def setVelocities(self, velocities):
        """Update atomic positions in attached DMS file
        """
        q = """UPDATE particle SET vx = ?1, vy = ?2, vz = ?3 WHERE id == ?4"""
        iat = 0
	#edit 3.10.15
        for i in self.receptor_atoms:
            (vx, vy , vz) = velocities[i].value_in_unit(angstrom/femtosecond)
            self._conn.execute(q, (vx,vy,vz,iat))
            iat += 1
        self._conn.commit()
	if self.ligand:
            jat = 0
	    for j in self.ligand_atoms:
            	(vx, vy , vz) = velocities[j].value_in_unit(angstrom/femtosecond)
            	self._conn1.execute(q, (vx,vy,vz,jat))
            	jat += 1
            self._conn1.commit()
	
	#edit end 3.10.15
        return iat


    def _get_gb_params(self):
        """
        get charge, radius, screened_radius from hct table in the .dms file to calculate
        --radiusN=radius-0.009  # Offset radius, the radius needs to be converted to nanometer unit
        --screenN=screened_radius*radiusN #Scaled offset radius 
        """
	length_conv = units.angstrom.conversion_factor_to(units.nanometer)

	gb_p = []

        q="""SELECT charge,radius,screened_radius from hct"""
               

        try:
            for charge,radius,screened_radius in self._conn.execute(q):
            
                chargeN = charge
                radiusN = radius
                screenN = screened_radius
                radiusN *=length_conv
                radiusN -=0.009
                screenN *=radiusN
                gb_p.append((chargeN,radiusN,screenN))
        except:
            return None

        #edited on 02/27/15                                                                             
        if self.ligand:

            q="""SELECT charge,radius,screened_radius from hct"""


            try:
                for charge,radius,screened_radius in self._conn1.execute(q):

                    chargeN = charge
                    radiusN = radius
                    screenN = screened_radius
                    radiusN *=length_conv
                    radiusN -=0.009
                    screenN *=radiusN
                    gb_p.append((chargeN,radiusN,screenN))
            except:
                return None

        #edited end on 02.27.15


        #edited on 02/26/15
        if self.bedam is not None:
            q="""SELECT charge,radius,screened_radius from hct"""
        
            
            try:
                for charge,radius,screened_radius in self._conn.execute(q):
            
                    chargeN = charge
                    radiusN = radius
                    screenN = screened_radius
                    radiusN *=length_conv
                    radiusN -=0.009
                    screenN *=radiusN
                    gb_p.append((chargeN,radiusN,screenN))
            except:
                return None
        
             #edited on 02/27/15                                                                                                                                                  
            if self.ligand:

                q="""SELECT charge,radius,screened_radius from hct"""


                try:
                    for charge,radius,screened_radius in self._conn1.execute(q):

                        chargeN = charge
                        radiusN = radius
                        screenN = screened_radius
                        radiusN *=length_conv
                        radiusN -=0.009
                        screenN *=radiusN
                        gb_p.append((chargeN,radiusN,screenN))
                except:
                    return None

            #edited end on 02.27.15
        #edit end on 02/26/15
        #for i in range(20):
        #    print i,
        #    print gb_p[i]
                
        return gb_p

    
    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=2.0*nanometer,
                     ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None, OPLS=False, implicitSolvent=None,SCORE=None):
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
         - softcore (object=None) if not None, the soft-core model to use.
        
        """
        self._checkForUnsupportedTerms()
        sys = mm.System()
        self.implicitSolvent=implicitSolvent
        #self.softcore = softcore
        self.score = SCORE
	
        # Buld the box dimensions
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        #TODO: how to manipulate the boxSize for duplicate system, right now, only consider implicit solvent model //02/26/15

        # Create all of the particles
        for mass in self._conn.execute('SELECT mass from particle'):
            sys.addParticle(mass[0]*dalton)
        
        if self.ligand:
            for mass in self._conn1.execute('SELECT mass from particle'):
                sys.addParticle(mass[0]*dalton)
        
        #if BEDAM is not None, duplicate the particles #added on 02/25/15
        #record the initial number for imaginary systems
        if self.bedam is not None:
            #self.bedamId = sys.getNumParticles()
            #if nonbondedCutoff is not None:
            #    self.bedamX = 5*nonbondedCutoff
            #else:
            #    raise ValueError('nonbondedCutoff needs to be set up for BEDAM')
            for mass in self._conn.execute('SELECT mass from particle'):
                sys.addParticle(mass[0]*dalton)

            if self.ligand:
                for mass in self._conn1.execute('SELECT mass from particle'):
                    sys.addParticle(mass[0]*dalton)


            #print self.bedamId      #tested on 02/25/15  
            #print self.bedamX


        # Add all of the forces
        self._addBondsToSystem(sys)
        self._addAnglesToSystem(sys)
        self._addConstraintsToSystem(sys)
        self._addPeriodicTorsionsToSystem(sys, OPLS)
        self._addImproperHarmonicTorsionsToSystem(sys)
        self._addCMAPToSystem(sys)
        self._addVirtualSitesToSystem(sys)
        
        nb,cnb,colnb,cbnb = self._addNonbondedForceToSystem(sys, OPLS)

        # Finish configuring the NonbondedForce.
        methodMap = {ff.NoCutoff:mm.NonbondedForce.NoCutoff,
                     ff.CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     ff.Ewald:mm.NonbondedForce.Ewald,
                     ff.PME:mm.NonbondedForce.PME}
                    
        if cnb is not None:
            cnb.setNonbondedMethod(methodMap[nonbondedMethod])
            cnb.setCutoffDistance(nonbondedCutoff)
        if colnb is not None:
            colnb.setNonbondedMethod(methodMap[nonbondedMethod])
            colnb.setCutoffDistance(nonbondedCutoff)
        if cbnb is not None:
            cbnb.setNonbondedMethod(methodMap[nonbondedMethod])
            cbnb.setCutoffDistance(nonbondedCutoff)
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)
        
	#add implicit solvent model. Right now, only HCT model is considered.
	if implicitSolvent is not None and implicitSolvent is not "BASIC":
            if implicitSolvent is HCT:
                gb_parms = self._get_gb_params()
                if gb_parms:
                    print('Adding HCT GB force ...')
                    gb = self.GBSAHCTForce()
                    #if SCORE:
                    #gb = self.GBSAHCTForce(cutoff=nonbondedCutoff/nanometer)
		    #gb = self.GBSAHCTForce(SA='ACE', cutoff=nonbondedCutoff/nanometer)
                    #else:
                    #    gb = GBSAHCTForce(SA='ACE', cutoff=nonbondedCutoff/angstrom)
                    for i in range(len(gb_parms)):	      
                        gb.addParticle(list(gb_parms[i])) #edit 3.13.15
	    # Set cutoff method //2.10.15 for setting up cutoff options in implicit solvent model
            else:
                raise ValueError('Right now, only HCT is supported')
            
            gb.setNonbondedMethod(methodMap[nonbondedMethod])
            gb.setCutoffDistance(nonbondedCutoff)
            #if nonbondedMethod is ff.NoCutoff:
            #    gb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            #    #gb.setNonbondedMethod(cnb.NoCutoff)
            #elif nonbondedMethod is ff.CutoffNonPeriodic:
            #    gb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            #    #gb.setNonbondedMethod(cnb.CutoffNonPeriodic)
            #    gb.setCutoffDistance(cutoff)
            #elif nonbondedMethod is ff.CutoffPeriodic:
            #    gb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            #    #gb.setNonbondedMethod(cnb.CutoffPeriodic)
            #    gb.setCutoffDistance(cutoff)
            #else:
            #    raise ValueError('Illegal nonbonded method for use with GBSA')
            #gb.setForceGroup(self.GB_FORCE_GROUP)
		
            sys.addForce(gb)

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

        #add a c-alpha restraint forces
        if self.carest:
            self._addCalpharestraint(sys)
        

        #edit 4.3.15
        #add restraint force between 109 and 154 atoms
        """
        if self.restraint is not None:
            restraint_force = mm.CustomBondForce("step(r-r0) * 0.5 * K * (r-r0)^2 * lambda_restraints ")
            #restraint_force.addGlobalParameter('lambda_restraints', self.lambda_restraint)
            restraint_force.addPerBondParameter("lambda_restraints")
            restraint_force.addPerBondParameter("K")
            restraint_force.addPerBondParameter("r0")
            restraint_force.addBond(109,154,[self.lambda_restraint,1225,0.5])
            restraint_force.setForceGroup(1)
            print "force group for restraint is"
            print restraint_force.getForceGroup()
            sys.addForce(restraint_force)
        """
        #edit end 4.3.15




        #edit 3.17.15
        """
        if self.restraint is not None:
            restraints = FlatBottomReceptorLigandRestraint(self.temperature, sys, self.positions, self.receptor_atoms, self.ligand_atoms,self.lambda_restraint)
            restraint_force = restraints.getRestraintForce() # Get force object incorporating restraints
            restraint_force.setForceGroup(1)
            print "force group for restraint is"
            print restraint_force.getForceGroup()
            sys.addForce(restraint_force)
            #metadata['standard_state_correction'] = restraints.getStandardStateCorrection() #standard state correction in kT
        """
        #edit end 3.17.15

        #edit 3.11.15
        print sys.getNumForces()
	for i in range(sys.getNumForces()):
	    print sys.getForce(i)
        #for force_index in range(sys.getNumForces()):
        #print sys.getForce(3)
        #print sys.getForce(3).getNumParticles()
        print sys.getForce(3).getCutoffDistance()
        #print sys.getForce(4)
        #print sys.getForce(4).getNumParticles()
        print sys.getForce(4).getCutoffDistance()
        #print sys.getForce(5)
        #print sys.getForce(5).getNumParticles()
        print sys.getForce(5).getCutoffDistance()
        #print sys.getForce(3).getNumParticles()
	#print sys.getForce(6).getCutoffDistance()
        
        #print sys.getForce(0)
        #print sys.getForce(1)
        #print sys.getForce(2)
	#print sys.getForce(3)
	#print sys.getForce(4)
	#print sys.getForce(5)
        #print sys.getForce(6)
        #if self.restraint is not None:
        #    print sys.getForce(6)
        
        #edit end 3.11.15
        #added on 11.23.15
        #add soft-core force between the real ligand atoms and the real receptor atoms
        #if self.softcore is not None:
            #self._modifyNonbondedForce(sys)
	#print sys.getNumForces()
        #edit end on 11.23.15
        
        return sys

    def _addBondsToSystem(self, sys):
        """Create the harmonic bonds
        """
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)

        q = """SELECT p0, p1, r0, fc, constrained FROM stretch_harm_term INNER JOIN stretch_harm_param ON stretch_harm_term.param=stretch_harm_param.id"""
        for p0, p1, r0, fc, constrained in self._conn.execute(q):
            if constrained:
                sys.addConstraint(p0, p1, r0*angstrom)
            else:
                # Desmond writes the harmonic bond force without 1/2
                # so we need to to double the force constant
                bonds.addBond(p0, p1, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

            # Record information that will be needed for constraining angles.
            self._atomBonds[p0][p1] = r0*angstrom
            self._atomBonds[p1][p0] = r0*angstrom

       
        if self.ligand:
            q = """SELECT p0, p1, r0, fc, constrained FROM stretch_harm_term INNER JOIN stretch_harm_param ON stretch_harm_term.param=stretch_harm_param.id"""
            for p0, p1, r0, fc, constrained in self._conn1.execute(q):
                if constrained:
                    sys.addConstraint(p0+self.fileId, p1+self.fileId, r0*angstrom)
                else:
                    # Desmond writes the harmonic bond force without 1/2                                                              
                    # so we need to to double the force constant                                                                      
                    bonds.addBond(p0+self.fileId, p1+self.fileId, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

                # Record information that will be needed for constraining angles.                                                     
                self._atomBonds[p0+self.fileId][p1+self.fileId] = r0*angstrom
                self._atomBonds[p1+self.fileId][p0+self.fileId] = r0*angstrom
       
        if self.bedam is not None:
            q = """SELECT p0, p1, r0, fc, constrained FROM stretch_harm_term INNER JOIN stretch_harm_param ON stretch_harm_term.param=stretch_harm_param.id"""
            for p0, p1, r0, fc, constrained in self._conn.execute(q):
                if constrained:
                    sys.addConstraint(p0+self.bedamId, p1+self.bedamId, r0*angstrom)
                else:
                    # Desmond writes the harmonic bond force without 1/2
                    # so we need to to double the force constant
                    bonds.addBond(p0+self.bedamId, p1+self.bedamId, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

                # Record information that will be needed for constraining angles.
                self._atomBonds[p0+self.bedamId][p1+self.bedamId] = r0*angstrom
                self._atomBonds[p1+self.bedamId][p0+self.bedamId] = r0*angstrom
       
            if self.ligand:
                q = """SELECT p0, p1, r0, fc, constrained FROM stretch_harm_term INNER JOIN stretch_harm_param ON stretch_harm_term.param=stretch_harm_param.id"""
                for p0, p1, r0, fc, constrained in self._conn1.execute(q):
                    if constrained:
                        sys.addConstraint(p0+self.file1Id, p1+self.file1Id, r0*angstrom)
                    else:
                        # Desmond writes the harmonic bond force without 1/2                                                                                                     
                        # so we need to to double the force constant                                                                                                             
                        bonds.addBond(p0+self.file1Id, p1+self.file1Id, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

                    # Record information that will be needed for constraining angles.                                                                                            
                    self._atomBonds[p0+self.file1Id][p1+self.file1Id] = r0*angstrom
                    self._atomBonds[p1+self.file1Id][p0+self.file1Id] = r0*angstrom
        
        return bonds

    def _addAnglesToSystem(self, sys):
        """Create the harmonic angles
        """
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        degToRad = math.pi/180

        q = """SELECT p0, p1, p2, theta0, fc, constrained FROM angle_harm_term INNER JOIN angle_harm_param ON angle_harm_term.param=angle_harm_param.id"""
        for p0, p1, p2, theta0, fc, constrained in self._conn.execute(q):
            if constrained:
                l1 = self._atomBonds[p1][p0]
                l2 = self._atomBonds[p1][p2]
                length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                sys.addConstraint(p0, p2, length)
                self._angleConstraints[p1][p0] = p2
                self._angleConstraints[p1][p2] = p0
            else:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angles.addAngle(p0, p1, p2, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)
        
        if self.ligand:
            q = """SELECT p0, p1, p2, theta0, fc, constrained FROM angle_harm_term INNER JOIN angle_harm_param ON angle_harm_term.param=angle_harm_param.id"""
            for p0, p1, p2, theta0, fc, constrained in self._conn1.execute(q):
                if constrained:
                    l1 = self._atomBonds[p1+self.fileId][p0+self.fileId]
                    l2 = self._atomBonds[p1+self.fileId][p2+self.fileId]
                    length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                    sys.addConstraint(p0+self.fileId, p2+self.fileId, length)
                    self._angleConstraints[p1+self.fileId][p0+self.fileId] = p2
                    self._angleConstraints[p1+self.fileId][p2+self.fileId] = p0
                else:
                    # Desmond writes the harmonic angle force without 1/2                                                             
                    # so we need to to double the force constant                                                                      
                    angles.addAngle(p0+self.fileId, p1+self.fileId, p2+self.fileId, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)
        
        if self.bedam is not None:
            q = """SELECT p0, p1, p2, theta0, fc, constrained FROM angle_harm_term INNER JOIN angle_harm_param ON angle_harm_term.param=angle_harm_param.id"""
            for p0, p1, p2, theta0, fc, constrained in self._conn.execute(q):
                if constrained:
                    l1 = self._atomBonds[p1+self.bedamId][p0+self.bedamId]
                    l2 = self._atomBonds[p1+self.bedamId][p2+self.bedamId]
                    length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                    sys.addConstraint(p0+self.bedamId, p2+self.bedamId, length)
                    self._angleConstraints[p1+self.bedamId][p0+self.bedamId] = p2
                    self._angleConstraints[p1+self.bedamId][p2+self.bedamId] = p0
                else:
                    # Desmond writes the harmonic angle force without 1/2
                    # so we need to to double the force constant
                    angles.addAngle(p0+self.bedamId, p1+self.bedamId, p2+self.bedamId, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)
        
            if self.ligand:
                q = """SELECT p0, p1, p2, theta0, fc, constrained FROM angle_harm_term INNER JOIN angle_harm_param ON angle_harm_term.param=angle_harm_param.id"""
                for p0, p1, p2, theta0, fc, constrained in self._conn1.execute(q):
                    if constrained:
                        l1 = self._atomBonds[p1+self.file1Id][p0+self.file1Id]
                        l2 = self._atomBonds[p1+self.file1Id][p2+self.file1Id]
                        length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                        sys.addConstraint(p0+self.file1Id, p2+self.file1Id, length)
                        self._angleConstraints[p1+self.file1Id][p0+self.file1Id] = p2
                        self._angleConstraints[p1+self.file1Id][p2+self.file1Id] = p0
                    else:
                        # Desmond writes the harmonic angle force without 1/2                                                                                                    
                        # so we need to to double the force constant                                                                                                             
                        angles.addAngle(p0+self.file1Id, p1+self.file1Id, p2+self.file1Id, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)
                    
        return angles

    def _addConstraintsToSystem(self, sys):
        """Add constraints to system. Normally these should already be
        added by the bonds table, but we want to make sure that there's
        no extra information in the constraints table that we're not
        including in the system"""
        for term_table in [n for n in list(self._tables.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
            param_table = term_table.replace('term', 'param')
            q = """SELECT p0, p1, r1 FROM %(term)s INNER JOIN %(param)s ON %(term)s.param=%(param)s.id""" % \
                {'term': term_table, 'param': param_table}
            for p0, p1, r1 in self._conn.execute(q):
                if not p1 in self._atomBonds[p0]:
                    sys.addConstraint(p0, p1, r1*angstrom)
                    self._atomBonds[p0][p1] = r1*angstrom
                    self._atomBonds[p1][p0] = r1*angstrom

        if 'constraint_hoh_term' in self._tables:
            degToRad = math.pi/180
            q = """SELECT p0, p1, p2, r1, r2, theta FROM constraint_hoh_term INNER JOIN constraint_hoh_param ON constraint_hoh_term.param=constraint_hoh_param.id"""
            for p0, p1, p2, r1, r2, theta in self._conn.execute(q):
                # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2
                # wihth O-H1 and O-H2 distances r1 and r2
                if not (self._angleConstraints[p0].get(p1, None) == p2):
                    length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                    sys.addConstraint(p1, p2, length)
        

        if self.ligand: 
            for term_table1 in [n for n in list(self._tables1.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
                param_table1 = term_table1.replace('term', 'param')
                q = """SELECT p0, p1, r1 FROM %(term)s INNER JOIN %(param)s ON %(term)s.param=%(param)s.id""" % \
                    {'term': term_table1, 'param': param_table1}
                for p0, p1, r1 in self._conn1.execute(q):
                    if not p1+self.fileId in self._atomBonds[p0+self.fileId]:
                        sys.addConstraint(p0+self.fileId, p1+self.fileId, r1*angstrom)
                        self._atomBonds[p0+self.fileId][p1+self.fileId] = r1*angstrom
                        self._atomBonds[p1+self.fileId][p0+self.fileId] = r1*angstrom

            if 'constraint_hoh_term' in self._tables1:
                degToRad = math.pi/180
                q = """SELECT p0, p1, p2, r1, r2, theta FROM constraint_hoh_term INNER JOIN constraint_hoh_param ON constraint_hoh_term.param=constraint_hoh_param.id"""
                for p0, p1, p2, r1, r2, theta in self._conn1.execute(q):
                    # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2                                                      
                    # wihth O-H1 and O-H2 distances r1 and r2                                                                         
                    if not (self._angleConstraints[p0+self.fileId].get(p1+self.fileId, None) == p2+self.fileId):
                        length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                        sys.addConstraint(p1+self.fileId, p2+self.fileId, length)
        
        if self.bedam is not None:
            for term_table in [n for n in list(self._tables.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
                param_table = term_table.replace('term', 'param')
                q = """SELECT p0, p1, r1 FROM %(term)s INNER JOIN %(param)s ON %(term)s.param=%(param)s.id""" % \
                    {'term': term_table, 'param': param_table}
                for p0, p1, r1 in self._conn.execute(q):
                    if not p1+self.bedamId in self._atomBonds[p0+self.bedamId]:
                        sys.addConstraint(p0+self.bedamId, p1+self.bedamId, r1*angstrom)
                        self._atomBonds[p0+self.bedamId][p1+self.bedamId] = r1*angstrom
                        self._atomBonds[p1+self.bedamId][p0+self.bedamId] = r1*angstrom

            if 'constraint_hoh_term' in self._tables:
                degToRad = math.pi/180
                q = """SELECT p0, p1, p2, r1, r2, theta FROM constraint_hoh_term INNER JOIN constraint_hoh_param ON constraint_hoh_term.param=constraint_hoh_param.id"""
                for p0, p1, p2, r1, r2, theta in self._conn.execute(q):
                    # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2
                    # wihth O-H1 and O-H2 distances r1 and r2
                    if not (self._angleConstraints[p0+self.bedamId].get(p1+self.bedamId, None) == p2+self.bedamId):
                        length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                        sys.addConstraint(p1+self.bedamId, p2+self.bedamId, length)
            
            if self.ligand:
                for term_table1 in [n for n in list(self._tables1.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
                    param_table1 = term_table1.replace('term', 'param')
                    q = """SELECT p0, p1, r1 FROM %(term)s INNER JOIN %(param)s ON %(term)s.param=%(param)s.id""" % \
                        {'term': term_table1, 'param': param_table1}
                    for p0, p1, r1 in self._conn1.execute(q):
                        if not p1+self.file1Id in self._atomBonds[p0+self.file1Id]:
                            sys.addConstraint(p0+self.file1Id, p1+self.file1Id, r1*angstrom)
                            self._atomBonds[p0+self.file1Id][p1+self.file1Id] = r1*angstrom
                            self._atomBonds[p1+self.file1Id][p0+self.file1Id] = r1*angstrom

                if 'constraint_hoh_term' in self._tables1:
                    degToRad = math.pi/180
                    q = """SELECT p0, p1, p2, r1, r2, theta FROM constraint_hoh_term INNER JOIN constraint_hoh_param ON constraint_hoh_term.param=constraint_hoh_param.id"""
                    for p0, p1, p2, r1, r2, theta in self._conn1.execute(q):
                        # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2  
                        # wihth O-H1 and O-H2 distances r1 and r2                                                                                         
                        if not (self._angleConstraints[p0+self.file1Id].get(p1+self.file1Id, None) == p2+self.file1Id):
                            length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                            sys.addConstraint(p1+self.file1Id, p2+self.file1Id, length)
     
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

        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 FROM dihedral_trig_term INNER JOIN dihedral_trig_param ON dihedral_trig_term.param=dihedral_trig_param.id"""
        for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn.execute(q):
            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue
                if OPLS:
                    periodic.addTorsion(p0, p1, p2, p3, [order, phi0*degree, fc*kilocalorie_per_mole])
                else:
                    periodic.addTorsion(p0, p1, p2, p3, order, phi0*degree, fc*kilocalorie_per_mole)

        if self.ligand:
            q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 FROM dihedral_trig_term INNER JOIN dihedral_trig_param ON dihedral_trig_term.param=dihedral_trig_param.id"""
            for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn1.execute(q):
                for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                    if fc == 0:
                        continue
                    if OPLS:
                        periodic.addTorsion(p0+self.fileId, p1+self.fileId, p2+self.fileId, p3+self.fileId, [order, phi0*degree, fc*kilocalorie_per_mole])
                    else:
                        periodic.addTorsion(p0+self.fileId, p1+self.fileId, p2+self.fileId, p3+self.fileId, order, phi0*degree, fc*kilocalorie_per_mole)

        
        if self.bedam is not None:
            q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 FROM dihedral_trig_term INNER JOIN dihedral_trig_param ON dihedral_trig_term.param=dihedral_trig_param.id"""
            for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn.execute(q):
                for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                    if fc == 0:
                        continue
                    if OPLS:
                        periodic.addTorsion(p0+self.bedamId, p1+self.bedamId, p2+self.bedamId, p3+self.bedamId, [order, phi0*degree, fc*kilocalorie_per_mole])
                    else:
                        periodic.addTorsion(p0+self.bedamId, p1+self.bedamId, p2+self.bedamId, p3+self.bedamId, order, phi0*degree, fc*kilocalorie_per_mole)
            
            if self.ligand:
                q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 FROM dihedral_trig_term INNER JOIN dihedral_trig_param ON dihedral_trig_term.param=dihedral_trig_param.id"""
                for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn1.execute(q):
                    for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                        if fc == 0:
                            continue
                        if OPLS:
                            periodic.addTorsion(p0+self.file1Id, p1+self.file1Id, p2+self.file1Id, p3+self.file1Id, [order, phi0*degree, fc*kilocalorie_per_mole])
                        else:
                            periodic.addTorsion(p0+self.file1Id, p1+self.file1Id, p2+self.file1Id, p3+self.file1Id, order, phi0*degree, fc*kilocalorie_per_mole)

    def _addImproperHarmonicTorsionsToSystem(self, sys):
        """Create the improper harmonic torsion terms
        """
        if not self._hasTable('improper_harm_term'):
            return

        harmonicTorsion = mm.CustomTorsionForce('k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)

        q = """SELECT p0, p1, p2, p3, phi0, fc FROM improper_harm_term INNER JOIN improper_harm_param ON improper_harm_term.param=improper_harm_param.id"""
        for p0, p1, p2, p3, phi0, fc in self._conn.execute(q):
            harmonicTorsion.addTorsion(p0, p1, p2, p3, [phi0*degree, fc*kilocalorie_per_mole])
        
        if self.ligand:
            if not self._hasTable1('improper_harm_term'):
                return
            q = """SELECT p0, p1, p2, p3, phi0, fc FROM improper_harm_term INNER JOIN improper_harm_param ON improper_harm_term.param=improper_harm_param.id"""
            for p0, p1, p2, p3, phi0, fc in self._conn1.execute(q):
                harmonicTorsion.addTorsion(p0+self.fileId, p1+self.fileId, p2+self.fileId, p3+self.fileId, [phi0*degree, fc*kilocalorie_per_mole])
     
        if self.bedam is not None:
            q = """SELECT p0, p1, p2, p3, phi0, fc FROM improper_harm_term INNER JOIN improper_harm_param ON improper_harm_term.param=improper_harm_param.id"""
            for p0, p1, p2, p3, phi0, fc in self._conn.execute(q):
                harmonicTorsion.addTorsion(p0+self.bedamId, p1+self.bedamId, p2+self.bedamId, p3+self.bedamId, [phi0*degree, fc*kilocalorie_per_mole])
       
            if self.ligand:
                if not self._hasTable1('improper_harm_term'):
                    return
                q = """SELECT p0, p1, p2, p3, phi0, fc FROM improper_harm_term INNER JOIN improper_harm_param ON improper_harm_term.param=improper_harm_param.id"""
                for p0, p1, p2, p3, phi0, fc in self._conn1.execute(q):
                    harmonicTorsion.addTorsion(p0+self.file1Id, p1+self.file1Id, p2+self.file1Id, p3+self.file1Id, [phi0*degree, fc*kilocalorie_per_mole])
         #edit end on 02/27/15 

    def _addCMAPToSystem(self, sys):
        """Create the CMAP terms
        """
        if not self._hasTable('torsiontorsion_cmap_term'):
            return

        # Create CMAP torsion terms
        cmap = mm.CMAPTorsionForce()
        sys.addForce(cmap)
        cmap_indices = {}

        for name in [k for k in list(self._tables.keys()) if k.startswith('cmap')]:
            size2 = self._conn.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
            fsize = math.sqrt(size2)
            if fsize != int(fsize):
                raise ValueError('Non-square CMAPs are not supported')
            size = int(fsize)

            map = [0 for i in range(size2)]
            for phi, psi, energy in self._conn.execute("SELECT phi, psi, energy FROM %s" % name):
                i = int((phi % 360) / (360.0 / size))
                j = int((psi % 360) / (360.0 / size))
                map[i+size*j] = energy
            index = cmap.addMap(size, map*kilocalorie_per_mole)
            cmap_indices[name] = index

        q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
        for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn.execute(q):
            cmap.addTorsion(cmap_indices[cmapid], p0, p1, p2, p3, p4, p5, p6, p7)

        if self.ligand:
            if not self._hasTable1('torsiontorsion_cmap_term'):
                return
            for name in [k for k in list(self._tables1.keys()) if k.startswith('cmap')]:
                size2 = self._conn1.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
                fsize = math.sqrt(size2)
                if fsize != int(fsize):
                    raise ValueError('Non-square CMAPs are not supported')
                size = int(fsize)

                map = [0 for i in range(size2)]
                for phi, psi, energy in self._conn1.execute("SELECT phi, psi, energy FROM %s" % name):
                    i = int((phi % 360) / (360.0 / size))
                    j = int((psi % 360) / (360.0 / size))
                    map[i+size*j] = energy
                index = cmap.addMap(size, map*kilocalorie_per_mole)
                cmap_indices[name] = index

            q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
            for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn1.execute(q):
                cmap.addTorsion(cmap_indices[cmapid], p0+self.fileId, p1+self.fileId, p2+self.fileId, p3+self.fileId, p4+self.fileId, p5+self.fileId, p6+self.fileId, p7+self.fileId)

        if self.bedam is not None:
            for name in [k for k in list(self._tables.keys()) if k.startswith('cmap')]:
                size2 = self._conn.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
                fsize = math.sqrt(size2)
                if fsize != int(fsize):
                    raise ValueError('Non-square CMAPs are not supported')
                size = int(fsize)

                map = [0 for i in range(size2)]
                for phi, psi, energy in self._conn.execute("SELECT phi, psi, energy FROM %s" % name):
                    i = int((phi % 360) / (360.0 / size))
                    j = int((psi % 360) / (360.0 / size))
                    map[i+size*j] = energy
                index = cmap.addMap(size, map*kilocalorie_per_mole)
                cmap_indices[name] = index

            q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
            for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn.execute(q):
                cmap.addTorsion(cmap_indices[cmapid], p0+self.bedamId, p1+self.bedamId, p2+self.bedamId, p3+self.bedamId, p4+self.bedamId, p5+self.bedamId, p6+self.bedamId, p7+self.bedamId)

            if self.ligand:
                if not self._hasTable1('torsiontorsion_cmap_term'):
                    return
                for name in [k for k in list(self._tables1.keys()) if k.startswith('cmap')]:
                    size2 = self._conn1.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
                    fsize = math.sqrt(size2)
                    if fsize != int(fsize):
                        raise ValueError('Non-square CMAPs are not supported')
                    size = int(fsize)

                    map = [0 for i in range(size2)]
                    for phi, psi, energy in self._conn1.execute("SELECT phi, psi, energy FROM %s" % name):
                        i = int((phi % 360) / (360.0 / size))
                        j = int((psi % 360) / (360.0 / size))
                        map[i+size*j] = energy
                    index = cmap.addMap(size, map*kilocalorie_per_mole)
                    cmap_indices[name] = index

                q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
                for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn1.execute(q):
                    cmap.addTorsion(cmap_indices[cmapid], p0+self.file1Id, p1+self.file1Id, p2+self.file1Id, p3+self.file1Id, p4+self.file1Id, p5+self.file1Id, p6+self.file1Id, p7+self.file1Id)

    def _addNonbondedForceToSystem(self, sys, OPLS):
        """Create the nonbonded force
        """
        cnb = None
	w = self.wscore
        
        colnb = None
        cbnb = None
 
        nb = mm.NonbondedForce()
        sys.addForce(nb)

	if self.implicitSolvent is not "BASIC":
            
	    if self.score is None:
                cnb = mm.CustomNonbondedForce("4.0*epsilon12*((sigma12/r)^12 - (sigma12/r)^6); sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2)")
                cnb.addPerParticleParameter("sigma")
                cnb.addPerParticleParameter("epsilon")
                sys.addForce(cnb)
	    else:
                cnb = mm.CustomNonbondedForce("4.0*epsilon12*((sigma12/reff)^12 - (sigma12/reff)^6);sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2); reff = sqrt(r*r+w*w); w = %f" %w)
                cnb.addPerParticleParameter("sigma")
                cnb.addPerParticleParameter("epsilon")
                sys.addForce(cnb)
                colnb = mm.CustomNonbondedForce("ONE_4PI_EPS0*chargeprod/reff; chargeprod = charge1*charge2; reff = sqrt(r*r+w*w); ONE_4PI_EPS0 = %f; w=%f" %(ONE_4PI_EPS0,w))
                colnb.addPerParticleParameter("charge")
                sys.addForce(colnb)

        else:
            if self.score is None:
                cnb = mm.CustomNonbondedForce("4.0*epsilon12*((sigma12/r)^12 - (sigma12/r)^6); sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2)")
                cnb.addPerParticleParameter("sigma")
                cnb.addPerParticleParameter("epsilon")
                sys.addForce(cnb)
                cbnb = mm.CustomNonbondedForce("ONE_4PI_EPS0*chargeprod/(40.0*r*r); chargeprod = charge1*charge2; ONE_4PI_EPS0 = %f" %(ONE_4PI_EPS0))
                cbnb.addPerParticleParameter("charge")
                sys.addForce(cbnb)

            else:
                cnb = mm.CustomNonbondedForce("4.0*epsilon12*((sigma12/reff)^12 - (sigma12/reff)^6);sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2); reff = sqrt(r*r+w*w); w = %f" %w)
                cnb.addPerParticleParameter("sigma")
                cnb.addPerParticleParameter("epsilon")
                sys.addForce(cnb)
                colnb = mm.CustomNonbondedForce("ONE_4PI_EPS0*chargeprod/(40.0*reff*reff); chargeprod = charge1*charge2; reff = sqrt(r*r+w*w); ONE_4PI_EPS0 = %f; w=%f" %(ONE_4PI_EPS0,w))
                colnb.addPerParticleParameter("charge")
                sys.addForce(colnb)

        q = """SELECT charge, sigma, epsilon FROM particle INNER JOIN nonbonded_param ON particle.nbtype=nonbonded_param.id"""
        for charge, sigma, epsilon in self._conn.execute(q):
            cnb.addParticle([sigma*angstrom, epsilon*kilocalorie_per_mole])   
            if colnb is not None:
                colnb.addParticle([charge])
            epsilon = 0
            if self.score or self.implicitSolvent is "BASIC":
                charge = 0
            nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)
            if cbnb is not None:
                cbnb.addParticle([charge]) 
          
        for p0, p1 in self._conn.execute('SELECT p0, p1 FROM exclusion'):
            cnb.addExclusion(p0, p1)
            nb.addException(p0, p1, 0.0, 1.0, 0.0)
            if colnb is not None:
                colnb.addExclusion(p0, p1)
            if cbnb is not None:
                cbnb.addExclusion(p0, p1)
  
        q = """SELECT p0, p1, aij, bij, qij FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param ON pair_12_6_es_term.param=pair_12_6_es_param.id;"""
        for p0, p1, a_ij, b_ij, q_ij in self._conn.execute(q):
            a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2

            if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                new_epsilon = 0
                new_sigma = 1
            else:
                new_epsilon =  b_ij**2/(4*a_ij)
                new_sigma = (a_ij / b_ij)**(1.0/6.0)
            
            nb.addException(p0, p1, q_ij, new_sigma, new_epsilon, True)
            
        n_total = self._conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
        n_in_exclusions = self._conn.execute("""SELECT COUNT(*) FROM exclusion INNER JOIN pair_12_6_es_term
        ON (exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
        OR (exclusion.p1==pair_12_6_es_term.p0 AND exclusion.p0==pair_12_6_es_term.p1)""").fetchone()
        if not n_total == n_in_exclusions:
            raise NotImplementedError('All pair_12_6_es_term must have a corresponding exclusion')

        if self.ligand:
            q = """SELECT charge, sigma, epsilon FROM particle INNER JOIN nonbonded_param ON particle.nbtype=nonbonded_param.id"""
            for charge, sigma, epsilon in self._conn1.execute(q):
                cnb.addParticle([sigma*angstrom, epsilon*kilocalorie_per_mole])
                if colnb is not None:
                    colnb.addParticle([charge])
                if cbnb is not None:
                    cbnb.addParticle([charge])
                epsilon = 0
                if self.score or self.implicitSolvent is "BASIC":
                    charge = 0
                nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)

            for p0, p1 in self._conn1.execute('SELECT p0, p1 FROM exclusion'):
                cnb.addExclusion(p0+self.fileId, p1+self.fileId)
                nb.addException(p0+self.fileId, p1+self.fileId, 0.0, 1.0, 0.0)
                if colnb is not None:
                    colnb.addExclusion(p0+self.fileId, p1+self.fileId)
                if cbnb is not None:
                    cbnb.addExclusion(p0+self.fileId, p1+self.fileId)
                

            q = """SELECT p0, p1, aij, bij, qij FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param ON pair_12_6_es_term.param=pair_12_6_es_param.id;"""
            for p0, p1, a_ij, b_ij, q_ij in self._conn1.execute(q):
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2

                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
                #if not SCORE and self.implicitSolvent is not "BASIC":
                nb.addException(p0+self.fileId, p1+self.fileId, q_ij, new_sigma, new_epsilon, True)

            n_total = self._conn1.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = self._conn1.execute("""SELECT COUNT(*) FROM exclusion INNER JOIN pair_12_6_es_term 
            ON (exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1) 
            OR (exclusion.p1==pair_12_6_es_term.p0 AND exclusion.p0==pair_12_6_es_term.p1)""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_term must have a corresponding exclusion')

        
        if self.bedam is not None:
            q = """SELECT charge, sigma, epsilon FROM particle INNER JOIN nonbonded_param ON particle.nbtype=nonbonded_param.id"""        
            for charge, sigma, epsilon in self._conn.execute(q):
                cnb.addParticle([sigma*angstrom, epsilon*kilocalorie_per_mole])
                if colnb is not None:
                    colnb.addParticle([charge])
                if cbnb is not None:
                    cbnb.addParticle([charge])
                epsilon = 0
                if self.score or self.implicitSolvent is "BASIC" :
                    charge = 0
                nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)
      
      
            for p0, p1 in self._conn.execute('SELECT p0, p1 FROM exclusion'):                
                cnb.addExclusion(p0+self.bedamId, p1+self.bedamId)
                nb.addException(p0+self.bedamId, p1+self.bedamId, 0.0, 1.0, 0.0)
                if colnb is not None:
                    colnb.addExclusion(p0+self.bedamId, p1+self.bedamId)
                if cbnb is not None:
                    cbnb.addExclusion(p0+self.bedamId, p1+self.bedamId)
                
                
            q = """SELECT p0, p1, aij, bij, qij
            FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
            ON pair_12_6_es_term.param=pair_12_6_es_param.id;"""
            for p0, p1, a_ij, b_ij, q_ij in self._conn.execute(q):
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2

                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
  
                nb.addException(p0+self.bedamId, p1+self.bedamId, q_ij, new_sigma, new_epsilon, True)
            
            n_total = self._conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = self._conn.execute("""SELECT COUNT(*)
            FROM exclusion INNER JOIN pair_12_6_es_term
            ON (exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
            OR (exclusion.p1==pair_12_6_es_term.p0 AND exclusion.p0==pair_12_6_es_term.p1)""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_term must have a corresponding exclusion')


            if self.ligand:
                q = """SELECT charge, sigma, epsilon FROM particle INNER JOIN nonbonded_param ON particle.nbtype=nonbonded_param.id"""
                for charge, sigma, epsilon in self._conn1.execute(q):
                    cnb.addParticle([sigma*angstrom, epsilon*kilocalorie_per_mole])
                    if colnb is not None:
                        colnb.addParticle([charge])
                    if cbnb is not None:
                        cbnb.addParticle([charge])
                    epsilon = 0
                    if self.score or self.implicitSolvent is "BASIC":
                        charge = 0
                    nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)
                    
                
                for p0, p1 in self._conn1.execute('SELECT p0, p1 FROM exclusion'):                    
                    cnb.addExclusion(p0+self.file1Id, p1+self.file1Id)
                    nb.addException(p0+self.file1Id, p1+self.file1Id, 0.0, 1.0, 0.0)
                    if colnb is not None:
                        colnb.addExclusion(p0+self.file1Id, p1+self.file1Id)
                    if cbnb is not None:
                        cbnb.addExclusion(p0+self.file1Id, p1+self.file1Id)


                q = """SELECT p0, p1, aij, bij, qij FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param ON pair_12_6_es_term.param=pair_12_6_es_param.id;"""
                for p0, p1, a_ij, b_ij, q_ij in self._conn1.execute(q):
                    a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                    b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                    q_ij = q_ij*elementary_charge**2

                    if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                        new_epsilon = 0
                        new_sigma = 1
                    else:
                        new_epsilon =  b_ij**2/(4*a_ij)
                        new_sigma = (a_ij / b_ij)**(1.0/6.0)
                    #if not SCORE and self.implicitSolvent is not "BASIC":
                    nb.addException(p0+self.file1Id, p1+self.file1Id, q_ij, new_sigma, new_epsilon, True)

                n_total = self._conn1.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
                n_in_exclusions = self._conn1.execute("""SELECT COUNT(*)                                                                   
                FROM exclusion INNER JOIN pair_12_6_es_term                                                                               
                ON (exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)                                            
                OR (exclusion.p1==pair_12_6_es_term.p0 AND exclusion.p0==pair_12_6_es_term.p1)""").fetchone()
                if not n_total == n_in_exclusions:
                    raise NotImplementedError('All pair_12_6_es_term must have a corresponding exclusion')

        return nb,cnb,colnb,cbnb

    def _addVirtualSitesToSystem(self, sys):
        """Create any virtual sites in the systempy
        """
        if not any(t.startswith('virtual_') for t in list(self._tables.keys())):
            return

        if 'virtual_lc2_term' in self._tables:
            q = """SELECT p0, p1, p2, c1 FROM virtual_lc2_term INNER JOIN virtual_lc2_param ON virtual_lc2_term.param=virtual_lc2_param.id"""
            for p0, p1, p2, c1 in self._conn.execute(q):
                vsite = mm.TwoParticleAverageSite(p1, p2, (1-c1), c1)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_lc3_term' in self._tables:
            q = """SELECT p0, p1, p2, p3, c1, c2 FROM virtual_lc3_term INNER JOIN virtual_lc3_param ON virtual_lc3_term.param=virtual_lc3_param.id"""
            for p0, p1, p2, p3, c1, c2 in self._conn.execute(q):
                vsite = mm.ThreeParticleAverageSite(p1, p2, p3, (1-c1-c2), c1, c2)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_out3_term' in self._tables:
            q = """SELECT p0, p1, p2, p3, c1, c2, c3 FROM virtual_out3_term INNER JOIN virtual_out3_param ON virtual_out3_term.param=virtual_out3_param.id"""
            for p0, p1, p2, p3, c1, c2, c3 in self._conn.execute(q):
                vsite = mm.OutOfPlaneSite(p1, p2, p3, c1, c2, c3)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_fdat3_term' in self._tables:
            raise NotImplementedError('OpenMM does not currently support '
                                      'fdat3-style virtual sites')

        #TODO: edit addVirtual to include BEDAM #02/26/15

    def _hasTable(self, table_name):
        """Does our DMS file contain this table?
        """
        return table_name in self._tables
    
    def _hasTable1(self, table_name):
        """Does our DMS file contain this table?                                                                                      
        """
        return table_name in self._tables1

    def _readSchemas(self):
        """Read the schemas of each of the tables in the dms file, populating
        the `_tables` instance attribute
        """
        tables = {}
        for table in self._conn.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in self._conn.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        self._tables = tables

        if self.ligand:
            tables1 = {}
            for table1 in self._conn1.execute("SELECT name FROM sqlite_master WHERE type='table'"):
                names1 = []
                for e in self._conn1.execute('PRAGMA table_info(%s)' % table):
                    names1.append(str(e[1]))
                tables1[str(table[0])] = names1
            self._tables1 = tables1


    def _checkForUnsupportedTerms(self):
        """Check the file for forcefield terms that are not currenty supported,
        raising a NotImplementedError
        """
        if 'posre_harm_term' in self._tables:
            raise NotImplementedError('Position restraints are not implemented.')
        flat_bottom_potential_terms = ['stretch_fbhw_term', 'angle_fbhw_term',
                                       'improper_fbhw_term', 'posre_fbhw_term']
        if any((t in self._tables) for t in flat_bottom_potential_terms):
            raise NotImplementedError('Flat bottom potential terms '
                                      'are not implemeneted')

        nbinfo = dict(list(zip(self._tables['nonbonded_info'],
                          self._conn.execute('SELECT * FROM nonbonded_info').fetchone())))

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

        if 'alchemical_particle' in self._tables:
            raise NotImplementedError('Alchemical particles are not supported')
        if 'alchemical_stretch_harm' in self._tables:
            raise NotImplementedError('Alchemical bonds are not supported')

        if 'polar_term' in self._tables:
            if self._conn.execute("SELECT COUNT(*) FROM polar_term").fetchone()[0] != 0:
                raise NotImplementedError('Drude particles are not currently supported')

        #add on 2.27.15
        #if self.ligfile:
        #    if 'posre_harm_term' in self._tables1:
        #        raise NotImplementedError('Position restraints are not implemented.')
        #    flat_bottom_potential_terms1 = ['stretch_fbhw_term', 'angle_fbhw_term',
        #                                   'improper_fbhw_term', 'posre_fbhw_term']
        #    if any((t in self._tables1) for t in flat_bottom_potential_terms1):
        #        raise NotImplementedError('Flat bottom potential terms '
        #                                  'are not implemeneted')
        #
        #    nbinfo1 = dict(list(zip(self._tables1['nonbonded_info'],
        #                           self._conn1.execute('SELECT * FROM nonbonded_info').fetchone())))
        #
        #    if nbinfo1['vdw_funct'] != 'vdw_12_6':
        #        raise NotImplementedError('Only Leonard-Jones van der Waals '
        #                                  'interactions are currently supported')
        #    if nbinfo1['vdw_rule'] != 'arithmetic/geometric':
        #        """#raise NotImplementedError('Only Lorentz-Berthelot nonbonded '                                                           #
        #        'combining rules are currently supported')"""
        #        nbinfo1['vdw_rule'] = 'arithmetic/geometric'
        #
        #        """#if 'nonbonded_combined_param' in self._tables:                                                                          #   
        #             raise NotImplementedError('nonbonded_combined_param interactions '                                                     #   
        #  'are not currently supported')"""
        #
        #    if 'alchemical_particle' in self._tables1:
        #        raise NotImplementedError('Alchemical particles are not supported')
        #    if 'alchemical_stretch_harm' in self._tables1:
        #        raise NotImplementedError('Alchemical bonds are not supported')
        #
        #    if 'polar_term' in self._tables1:
        #        if self._conn1.execute("SELECT COUNT(*) FROM polar_term").fetchone()[0] != 0:
        #            raise NotImplementedError('Drude particles are not currently supported')


    def get_orig_force_parameters(self,sys):
	
	force_nb = sys.getForce(3)
        force_cnb = sys.getForce(4)
	if self.implicitSolvent is HCT:
	    if self.score is None:
                force_cgb = sys.getForce(5)
            else:
                force_colnb = sys.getForce(5)
		force_cgb = sys.getForce(6)          
        elif self.implicitSolvent is "BASIC":
            if self.score is None:
                force_cbnb = sys.getForce(5)
            else:
                force_colnb = sys.getForce(5)
	elif self.implicitSolvent is None:
            if self.score:
                force_colnb = sys.getForce(5)

              
        for i in self.total_atoms:
	    val = force_nb.getParticleParameters(i) 
	    self.nb_orig_para.append((val[0]._value,val[1]._value,val[2]._value))
            self.cnb_orig_para.append(force_cnb.getParticleParameters(i))
            if self.implicitSolvent is HCT:
                self.cgb_orig_para.append(force_cgb.getParticleParameters(i))
                if self.score:
                    self.colnb_orig_para.append(force_colnb.getParticleParameters(i))
            elif self.implicitSolvent is "BASIC":
                if self.score:
                    self.colnb_orig_para.append(force_colnb.getParticleParameters(i))
                else:
                    self.cbnb_orig_para.append(force_cbnb.getParticleParameters(i))
            elif self.implicitSolvent is None:
                if self.score:
                    self.colnb_orig_para.append(force_colnb.getParticleParameters(i))
                    

    def set_orig_force_parameters(self,sys,context):
	
	force_nb = sys.getForce(3)
        force_cnb = sys.getForce(4)
        if self.implicitSolvent is HCT:
	    if self.score is None:
                force_cgb = sys.getForce(5)
            else:
                force_colnb = sys.getForce(5)
		force_cgb = sys.getForce(6)          
        elif self.implicitSolvent is "BASIC":
            if self.score:
                force_colnb = sys.getForce(5)
            else:
                force_cbnb = sys.getForce(5)
	elif self.implicitSolvent is None:
            if self.score:
                force_colnb = sys.getForce(5)


        for i in self.total_atoms:
            val = list(self.nb_orig_para[i])
            force_nb.setParticleParameters(i,val[0],val[1],val[2])
            force_cnb.setParticleParameters(i,list(self.cnb_orig_para[i]))
            if self.implicitSolvent is HCT:
                force_cgb.setParticleParameters(i,list(self.cgb_orig_para[i]))
                if self.score:
                    force_colnb.setParticleParameters(i,list(self.colnb_orig_para[i]))
            elif self.implicitSolvent is "BASIC":
                if self.score:
                    force_colnb.setParticleParameters(i,list(self.colnb_orig_para[i]))
                else:
                    force_cbnb.setParticleParameters(i,list(self.cbnb_orig_para[i]))
            elif self.implicitSolvent is None:
                if self.score:
                    force_colnb.setParticleParameters(i,list(self.colnb_orig_para[i]))

        force_nb.updateParametersInContext(context)
        force_cnb.updateParametersInContext(context)
        if self.implicitSolvent is HCT:
            force_cgb.updateParametersInContext(context)
            if self.score:
                force_colnb.updateParametersInContext(context)
        elif self.implicitSolvent is "BASIC":
            if self.score:
                force_colnb.updateParametersInContext(context)
            else:
                force_cbnb.updateParametersInContext(context)
        elif self.implicitSolvent is None:
            if self.score:
                force_colnb.updateParametersInContext(context)


    def binding_energy_calculation(self,sys,context,original=None):
        
        #change the parameters for forces(including GB force, nonbonded force right now)
        #first get forces from sys, then change the parameters of the forces to zero them.

	force_nb = sys.getForce(3)
        force_cnb = sys.getForce(4)
        if self.implicitSolvent is HCT:
	    if self.score is None:
                force_cgb = sys.getForce(5)
            else:
                force_colnb = sys.getForce(5)
		force_cgb = sys.getForce(6)          
        elif self.implicitSolvent is "BASIC":
            if self.score:
                force_colnb = sys.getForce(5)
            else:
                force_cbnb = sys.getForce(5)
	elif self.implicitSolvent is None:
            if self.score:
                force_colnb = sys.getForce(5)
	

        if original is not None:
            for i in self.original_atoms:
                force_nb.setParticleParameters(i,0.0,0.0,0.0)
                force_cnb.setParticleParameters(i,[0.0,0.0])
                if self.implicitSolvent is HCT:
                    force_cgb.setParticleParameters(i,[0.0,-0.009])
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])
                elif self.implicitSolvent is "BASIC":
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])
                    else:
                        force_cbnb.setParticleParameters(i,[0.0])
                elif self.implicitSolvent is None:
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])

        else:
            for i in self.image_atoms:
                force_nb.setParticleParameters(i,0.0,0.0,0.0)
                force_cnb.setParticleParameters(i,[0.0,0.0])
                if self.implicitSolvent is HCT:
                    force_cgb.setParticleParameters(i,[0.0,-0.009])
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])
                elif self.implicitSolvent is "BASIC":
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])
                    else:
                        force_cbnb.setParticleParameters(i,[0.0])
                elif self.implicitSolvent is None:
                    if self.score:
                        force_colnb.setParticleParameters(i,[0.0])
                

        force_nb.updateParametersInContext(context)
        force_cnb.updateParametersInContext(context)
        if self.implicitSolvent is HCT:
            force_cgb.updateParametersInContext(context)
            if self.score:
                force_colnb.updateParametersInContext(context)
        elif self.implicitSolvent is "BASIC":
            if self.score:
                force_colnb.updateParametersInContext(context)
            else:
                force_cbnb.updateParametersInContext(context)
        elif self.implicitSolvent is None:
            if self.score:
                force_colnb.updateParametersInContext(context)
     
    #added on 11.25.15 for gb force blow out
    def GBSAHCTForce(self,solventDielectric=78.5, soluteDielectric=1, SA=None,
                 cutoff=None, kappa=0.0):

        w=self.wscore
        custom = CustomGBForce()

        custom.addPerParticleParameter("q")
        custom.addPerParticleParameter("or") # Offset radius
        custom.addPerParticleParameter("sr") # Scaled offset radius
        if self.score:
            custom.addComputedValue("I", "step(reff+sr2-or1)*0.5*(1/L-1/U+0.25*(reff-sr2^2/reff)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/reff);"
                                  "U=reff+sr2;"
                                  "L=max(or1, D);"
                                  "D=abs(reff-sr2);"
                                  "reff = sqrt(r*r+w*w);"
                                  "w = %f;"%w, CustomGBForce.ParticlePairNoExclusions)
        else:
            custom.addComputedValue("I", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                                  "U=r+sr2;"
                                  "L=max(or1, D);"
                                  "D=abs(r-sr2);", CustomGBForce.ParticlePairNoExclusions)
        custom.addComputedValue("B", "1/(1/or-I)", CustomGBForce.SingleParticle)
        #custom.addComputedValue("B", "0.15", CustomGBForce.SingleParticle)
        self._createEnergyTerms(custom, solventDielectric, soluteDielectric, SA, cutoff, kappa, 0.009)
        return custom


    def _createEnergyTerms(self,force, solventDielectric, soluteDielectric, SA, cutoff, kappa, offset):
        # Add the energy terms to the CustomGBForce.  These are identical for all the GB models.
    
        w = self.wscore

        params = "; solventDielectric=%.16g; soluteDielectric=%.16g; kappa=%.16g; offset=%.16g" % (solventDielectric, soluteDielectric, kappa, offset)
        if cutoff is not None:
            params += "; cutoff=%.16g" % cutoff
        if kappa > 0:
            force.addEnergyTerm("-0.5*138.935485*(1/soluteDielectric-exp(-kappa*B)/solventDielectric)*q^2/B"+params,
                CustomGBForce.SingleParticle)
        elif kappa < 0:
            # Do kappa check here to avoid repeating code everywhere
            raise ValueError('kappa/ionic strength must be >= 0')
        else:
            force.addEnergyTerm("-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B"+params,
                CustomGBForce.SingleParticle)
	    #force.addEnergyTerm("-0.0*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B"+params,
            #    CustomGBForce.SingleParticle)
        if SA=='ACE':
            force.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
        elif SA is not None:
            raise ValueError('Unknown surface area method: '+SA)
        if cutoff is None:
            if kappa > 0:
                force.addEnergyTerm("-138.935485*(1/soluteDielectric-exp(-kappa*f)/solventDielectric)*q1*q2/f;"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
            else:
                if self.score:
                    force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                                "f=sqrt(reff^2+B1*B2*exp(-reff^2/(4*B1*B2)));"
                                "reff = sqrt(r*r+w*w);"
                                "w=%f" %w+params, CustomGBForce.ParticlePairNoExclusions)
                else:
                    force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
        else:
            if kappa > 0:
                force.addEnergyTerm("-138.935485*(1/soluteDielectric-exp(-kappa*f)/solventDielectric)*q1*q2*(1/f-"+str(1/cutoff)+");"
                                "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
            else:
		if self.score:
                    force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2*(1/f-"+str(1/cutoff)+");"
                                "f=sqrt(reff^2+B1*B2*exp(-reff^2/(4*B1*B2)));"
                                "reff = sqrt(r*r+w*w);"
                                "w=%f" %w+params, CustomGBForce.ParticlePairNoExclusions)
		else:
		    force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2*(1/f-"+str(1/cutoff)+");"
                               "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
	            #force.addEnergyTerm("-0.0*138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                    #            "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)

    #added on 112315
    
    #added on 122815
    def _addCalpharestraint(self,sys):
        """add a c-alpha restraint, applying to all residues in protein"""
        force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 5.0*kilocalorie_per_mole/angstrom**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i in self.calphas:
            force.addParticle(i, self.calphas[i]*angstrom)
	    #print i,self.calphas[i]
        sys.addForce(force)


    """
    def _modifyNonbondedForce(self,sys):
        
        #Create a custom nonbonded force for soft-core application
        
        real_ligand_atoms = copy.deepcopy(self.ligand_atoms)
        real_receptor_atoms = copy.deepcopy(self.receptor_atoms)
        alchemical_atom_indices = real_ligand_atoms
        
        wv = 0.5#*angstrom #soft-core parameters
        
        
        #Create a copy of the NonbondedForce to handle non-alchemical interactions
        #nonbonded_force = copy.deepcopy(reference_force)
	#nb = mm.NonbondedForce()
        #reference_force =mm.NonbondedForce()
        #nonbonded_force = copy.deepcopy(reference_force)
	nonbonded_force = mm.NonbondedForce()
	reference_force = mm.NonbondedForce()
        #sys.addForce(nonbonded_force)
        
	#print("here is good")
        #create CustomNonbondedForce objects to handle softcore interactions between alchemically-modified system and rest of system
        
        #create atom groups
        natoms = sys.getNumParticles()
        atomset1 = set(alchemical_atom_indices) #only alchemically-modified atoms
        atomset2 = set(range(sys.getNumParticles()))# all atoms, including alchemical region
        #atomset2 = set(real_receptor_atoms)
        #CustomNonbondedForce energy expression
        sterics_energy_expression = ""
        electrostatics_energy_expression = ""
	
	
        #select functional form based on nonbonded method
        method = reference_force.getNonbondedMethod()
        if method in [openmm.NonbondedForce.NoCutoff]:
            #soft-core Lennard-Jones
            sterics_energy_expression += "U_sterics = 4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;"
            # soft-core Coulomb
            electrostatics_energy_expression += "U_electrostatics = ONE_4PI_EPS0*lambda_electrostatics*chargeprod/reff_electrostatics;"
	
        elif method in [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic]:
            # soft-core Lennard-Jones
            sterics_energy_expression += "U_sterics = lambda_sterics*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;"
            # reaction-field electrostatics
            epsilon_solvent = reference_force.getReactionFieldDielectric()
            r_cutoff = reference_force.getCutoffDistance()
            electrostatics_energy_expression += "U_electrostatics = lambda_electrostatics*ONE_4PI_EPS0*chargeprod*(reff_electrostatics^(-1) + k_rf*reff_electrostatics^2 - c_rf);"
            k_rf = r_cutoff**(-3) * ((epsilon_solvent - 1) / (2*epsilon_solvent + 1))
            c_rf = r_cutoff**(-1) * ((3*epsilon_solvent) / (2*epsilon_solvent + 1))          
            electrostatics_energy_expression += "k_rf = %f;" % (k_rf / k_rf.in_unit_system(unit.md_unit_system).unit)
            electrostatics_energy_expression += "c_rf = %f;" % (c_rf / c_rf.in_unit_system(unit.md_unit_system).unit)
        else:
            raise Exception("Nonbonded method %s not supported yet." % str(method))
        
        # Add additional definitions common to all methods
        sterics_energy_expression += "reff_sterics = sqrt(r*r+w*w);"#sigma*((softcore_alpha*(1.-lambda_sterics) + (r/sigma)^6))^(1/6);" # effective softcore distance for sterics
        #sterics_energy_expression += "softcore_alpha = %f;" % softcore_alpha
        sterics_energy_expression += "w = %f;" %wv
        electrostatics_energy_expression += "reff_electrostatics = sqrt(r*r+w*w);" #sqrt(softcore_beta*(1.-lambda_electrostatics) + r^2);" # effective softcore distance for electrostatics
        #electrostatics_energy_expression += "softcore_beta = %f;" % (softcore_beta / softcore_beta.in_unit_system(unit.md_unit_system).unit)
        electrostatics_energy_expression += "ONE_4PI_EPS0 = %f;" % ONE_4PI_EPS0 # already in OpenMM units
        electrostatics_energy_expression += "w=%f;" %wv
        # Define mixing rules.
        sterics_mixing_rules = ""
        sterics_mixing_rules += "epsilon = sqrt(epsilon1*epsilon2);" # mixing rule for epsilon
        sterics_mixing_rules += "sigma = sqrt(sigma1*sigma2);" # mixing rule for sigma
        electrostatics_mixing_rules = ""
        electrostatics_mixing_rules += "chargeprod = charge1*charge2;" # mixing rule for charges
        
        # Create CustomNonbondedForce to handle interactions between alchemically-modified atoms and rest of system.
        electrostatics_custom_nonbonded_force = openmm.CustomNonbondedForce("U_electrostatics;" + electrostatics_energy_expression + electrostatics_mixing_rules)
        electrostatics_custom_nonbonded_force.addGlobalParameter("lambda_electrostatics", 1.0)
        electrostatics_custom_nonbonded_force.addPerParticleParameter("charge") # partial charge
        sterics_custom_nonbonded_force = openmm.CustomNonbondedForce("U_sterics;" + sterics_energy_expression + sterics_mixing_rules)
        sterics_custom_nonbonded_force.addGlobalParameter("lambda_sterics", 1.0)
        sterics_custom_nonbonded_force.addPerParticleParameter("sigma") # Lennard-Jones sigma
        sterics_custom_nonbonded_force.addPerParticleParameter("epsilon") # Lennard-Jones epsilon
        
        # Set parameters to match reference force.
        #sterics_custom_nonbonded_force.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction())
        #electrostatics_custom_nonbonded_force.setUseSwitchingFunction(False)
        #sterics_custom_nonbonded_force.setCutoffDistance(nonbonded_force.getCutoffDistance())
        #electrostatics_custom_nonbonded_force.setCutoffDistance(nonbonded_force.getCutoffDistance())
        #sterics_custom_nonbonded_force.setSwitchingDistance(nonbonded_force.getSwitchingDistance())
        #sterics_custom_nonbonded_force.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())
        #electrostatics_custom_nonbonded_force.setUseLongRangeCorrection(False)

        # Set periodicity and cutoff parameters corresponding to reference Force.
        if nonbonded_force.getNonbondedMethod() in [openmm.NonbondedForce.Ewald, openmm.NonbondedForce.PME, openmm.NonbondedForce.CutoffPeriodic]:
            sterics_custom_nonbonded_force.setNonbondedMethod( openmm.CustomNonbondedForce.CutoffPeriodic )
            electrostatics_custom_nonbonded_force.setNonbondedMethod( openmm.CustomNonbondedForce.CutoffPeriodic )
        else:
            sterics_custom_nonbonded_force.setNonbondedMethod( nonbonded_force.getNonbondedMethod() )
            electrostatics_custom_nonbonded_force.setNonbondedMethod( nonbonded_force.getNonbondedMethod() )
        
        # Restrict interaction evaluation to be between alchemical atoms and rest of environment.
        # TODO: Exclude intra-alchemical region if we are separately handling that through a separate CustomNonbondedForce for decoupling.
        sterics_custom_nonbonded_force.addInteractionGroup(atomset1, atomset2)
        electrostatics_custom_nonbonded_force.addInteractionGroup(atomset1, atomset2)

        # Add custom forces.
        sys.addForce(sterics_custom_nonbonded_force)
        sys.addForce(electrostatics_custom_nonbonded_force)
        # Create CustomBondForce to handle exceptions for both kinds of interactions.
        #custom_bond_force = openmm.CustomBondForce("U_sterics + U_electrostatics;" + sterics_energy_expression + electrostatics_energy_expression)
        #custom_bond_force.addGlobalParameter("lambda_electrostatics", 1.0);
        #custom_bond_force.addGlobalParameter("lambda_sterics", 1.0);
        #custom_bond_force.addPerBondParameter("chargeprod") # charge product
        #custom_bond_force.addPerBondParameter("sigma") # Lennard-Jones effective sigma
        #custom_bond_force.addPerBondParameter("epsilon") # Lennard-Jones effective epsilon
        #sys.addForce(custom_bond_force)

        # Move NonbondedForce particle terms for alchemically-modified particles to CustomNonbondedForce.
        for particle_index in range(nonbonded_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
            # Add parameters to custom force handling interactions between alchemically-modified atoms and rest of system.
            sterics_custom_nonbonded_force.addParticle([sigma, epsilon])
            electrostatics_custom_nonbonded_force.addParticle([charge])
            # Turn off Lennard-Jones contribution from alchemically-modified particles.
            if particle_index in alchemical_atom_indices:
                nonbonded_force.setParticleParameters(particle_index, 0*charge, sigma, 0*epsilon)
	    #if particle_index in real_receptor_atoms:
		#nonbonded_force.setParticleParameters(particle_index, 0*charge, sigma, 0*epsilon)
        # Move NonbondedForce exception terms for alchemically-modified particles to CustomNonbondedForce/CustomBondForce.
        for exception_index in range(nonbonded_force.getNumExceptions()):
            # Retrieve parameters.
            [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
            # Exclude this atom pair in CustomNonbondedForce.
            sterics_custom_nonbonded_force.addExclusion(iatom, jatom)
            electrostatics_custom_nonbonded_force.addExclusion(iatom, jatom)
            # Move exceptions involving alchemically-modified atoms to CustomBondForce.
            #if self.annihilate_sterics and (iatom in alchemical_atom_indices) and (jatom in alchemical_atom_indices):
                # Add special CustomBondForce term to handle alchemically-modified Lennard-Jones exception.
            #    custom_bond_force.addBond(iatom, jatom, [chargeprod, sigma, epsilon])
                # Zero terms in NonbondedForce.
            #    nonbonded_force.setExceptionParameters(exception_index, iatom, jatom, 0*chargeprod, sigma, 0*epsilon)
        # TODO: Add back NonbondedForce terms for alchemical system needed in case of decoupling electrostatics or sterics via second CustomBondForce.
        # TODO: Also need to change current CustomBondForce to not alchemically disappear system.
	#print sys.getNumParticles()
	#print sys.getNumForces()
        return
        """


        
        #create a dict of parameters, right now, only alpha_i is created
        #binding_parameters = {
        #    #'binding_alpha' : binding_state
        #    'sigma': binding_state,
        #    'epsilon': binding_state
        #    }
        #set parameters in context
        #for parameter in binding_parameters:
        #    #right now only one parameter
        #    try:
        #        context.setParameter(parameter,binding_parameters[parameter])
        #    except Exception as e:
        #        pass
        #
        #return
    #edit end 3.10.15

    #edit on 3.16.15

    #edit on 4.3.15
    """
    def energy_decomposition(self,sys):
        
        all_names = dict()
        force_group_names = dict()
        energy_components = dict()
        kcal = u.kilocalories_per_mole

        for attr in dir(sys):
            if attr.endswith('_FORCE_GROUP'):
                val = getattr(sys,attr)
                all_names[val] = attr.replace('_FORCE_GROUP','').lower()
        for force in context.getSystem().getForces():
            gp = force.getForceGroup()
            force_group_names[gp] = all_names[gp]
        
        for grp, name in force_group_names.iteritems():
            state = context.getState(getEnergy=True, groups=1<<grp)
            energy_components[name] = state.getPotentialEnergy().value_in_unit(kcal)
        
        e = context.getState(getEnergy = True).getPotentialEnergy()
        energy_components['total'] = e.value_in_unit(kcal)

        return energy_components
    """
    #edit end on 4.3.15

    #added on 10.26.15
    def _closestAtomToCentroid(self,sys, positions, indices = None, masses=None):
        """
        Identify the closest atom to the centroid of the given coordinate set
        Parameters
        ----------
        coordinates: units.Quantity of natoms * 3 with units compatible with nanometers coordinates of object to identify atom closes to centroid
        List of atoms indices for which closest atom to centroid is to be computed.
        masses: simtk.unit.Quantity of natoms with units compatible with amu
            Masses of particles used to weight distance calculation, if not None (default : None)
        
        Returns
        ------
        closest_atom: int
            Index of atom closest to centroid of specified atoms
        """
        system = sys
        coordinates = units.Quantity(np.array(positions/positions.unit),positions.unit)
        indice = list(indices)

        if indices is not None:
            coordinates = coordinates[indice,:]

        # Get dimentionless coordinates
        x_unit = coordinates.unit
        x = coordinates / x_unit

        """
        # Get the masses of the system
        all_masses = []
        for atom_index in range(sys.getNumParticles()):
            mass = sys.getParticleMass(atom_index)
            xmass = mass/mass.unit
            all_masses.append(xmass)
        
        if indice1 and indice2 is not None:
            masses = all_masses[indice1:(indice2+indice1)]
        else:
            masses = all_masses
        """
        # Determine number of atoms
        natoms = x.shape[0]
        print natoms
        
        # Compute (natoms,1) array of normalized weights
        w = np.ones([natoms,1])
        #w = masses #/ masses.unit # (natoms,)array
        #w = np.reshape(w, (natoms,1)) #(natoms,1) array
        #w /=w.sum()
        
        # Compute centroid (still in dimensionless units)
        centroid = (np.tile(w, (1,3)) * x).sum(0) # (3,) array
        
        # Compute distances from centroid
        distances = np.sqrt(((x - np.tile(centroid,(natoms,1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i
        
        #Determine closest atom
        closest_atom = int(np.argmin(distances))
        
        if indices is not None:
            closest_atom = indices[closest_atom]
        
        return closest_atom

    def computeCOM(self,sys):
        
        receptor_atom = self._closestAtomToCentroid(sys,self.positions,self.receptor_atoms)
        ligand_atom = self._closestAtomToCentroid(sys,self.positions,self.ligand_atoms)
        print (receptor_atom,ligand_atom)
        return (receptor_atom,ligand_atom)
    #added end on 10.26.15

    def close(self):
        """Close the SQL connection
        """
        if self._open:
            self._conn.close()

    def __del__(self):
        self.close()



#=============================================================================================
# MODULE CONSTANTS
#=============================================================================================
#COPY AND EDIT FROM YANK

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant

class ReceptorLigandRestraint(object):
    """
    Impose a single restraint between ligand and protein to prevent ligand from drifting too far
    from protein in implicit solvent calculations.
    This restraint strength is controlled by a global context parameter called 'lambda_restraints'.
    NOTES
    To create a subclass that uses a different restraint energy function, follow these steps:
    * Redefine class variable 'energy_function' with the energy function of choice.
    * Redefine class variable 'bond_parameter_names' to list the parameters in 'energy_function'
    that must be computed for each system.
    * Redefine the _determineBondParameters() member function to compute these parameters and
    return them in a list in the same order as 'bond_parameter_names'.
    """
    def __init__(self, temperature, system, coordinates, receptor_atoms, ligand_atoms,lambda_restraint):
        """
        Initialize a receptor-ligand restraint class.
        Parameters
        ----------
        temperature : state's temperature
        The thermodynamic state specifying temperature, pressure, etc. to which restraints are to be added
        coordinates : simtk.unit.Quantity of natoms x 3 with units compatible with nanometers
        Reference coordinates to use for imposing restraints
        receptor_atoms : list of int
        A complete list of receptor atoms
        ligand_atoms : list of int
        A complete list of ligand atoms
        """
        #self.state = state
        self.system = system
        self.coordinates = units.Quantity(np.array(coordinates / coordinates.unit),coordinates.unit)
        self.receptor_atoms = list(receptor_atoms)
        self.ligand_atoms = list(ligand_atoms)

        self.temperature = temperature
        self.kT = kB * self.temperature #thermal energy
        self.beta = 1.0 / self.kT # inverse temperature
        self.lambda_restraint = lambda_restraint

    
        #Determine atoms closet to centroids on ligand and receptor.
        self.restrained_receptor_atom = self._closestAtomToCentroid(self.coordinates,self.receptor_atoms)
        self.restrained_ligand_atom = self._closestAtomToCentroid(self.coordinates,self.ligand_atoms)
        
        print ("restrained receptor atom: %d" % self.restrained_receptor_atom)
        print ("restrained ligand atom: %d" % self.restrained_ligand_atom)
        
        #Determine radius of gyration of receptor.
        self.radius_of_gyration = self._computeRadiusOfGyration(self.coordinates[self.receptor_atoms,:])

        #Determine parameters
        self.bond_parameters = self._determineBondParameters()
        
        #Determine standard state correction
        self.standard_state_correction = self._computeStandardStateCorrection()
        
        return

    @abstractmethod
    def _determineBondParameters(self):
        """
        Determine bond parameters for CustomBondForce between protein and ligand.
        Returns
        -------
        parameters : list
        List of parameters for CustomBondForce
        Notes
        -----
        The spring constant is selected to give 1 kT at one standard deviation of receptor atoms about the receptor
        restrained atom.
        """
        pass


    def _computeRadiusOfGyration(self,coordinates):
        
        """
        Compute the radius of gyration of the specified coordinate set.
        Parameters
        ----------
        coordinates : simtk.unit.Quantity with units compatible with angstrom
        The coordinate set (natoms x 3) for which the radius of gyration is to be computed.
        Returns
        -------
        radius_of_gyration : simtk.unit.Quantity with units compatible with angstrom
        The radius of gyration
        """


        unit = coordinates.unit

        # Get dimensionless receptor coordinates
        x = coordinates / unit

        # Get dimensionless restrained atom coordinate
        xref = x.mean(0)
        xref = np.reshape(xref,(1,3)) #(1,3) array

        # Compute distances from restrained atom.
        natoms = x.shape[0]
        distances = np.sqrt(((x - np.tile(xref, (natoms,1)))**2).sum(1)) # distance[i] is the distance from the centroid to particle i
        
        #Compute std dev of distances from restrained atom
        radius_of_gyration = distances.std() * unit

        return radius_of_gyration

    def _createRestraintForce(self, particle1, particle2, mm = None):
        """
        Create a new copy of the receptor-ligand restraint force.
        Parameters
        ----------
        particle1 : int
            Index of first particle for which restraint is to be applied.
        particle2 : int
            Index of second particle for which restraint is to be applied
        mm : simtk.openmm compliant interface, optional, default=None
        If specified, use an alternative OpenMM API implementation.
        Otherwise, use simtk.openmm.

        Returns
        -------
        force : simtk.openmm.CustomBondForce
        A restraint force object
        """

        if mm is None: mm = openmm
        
        force = openmm.CustomBondForce(self.energy_function)
        force.addGlobalParameter('lambda_restraints', self.lambda_restraint) #edit on 4.6.15
        for parameter in self.bond_parameter_names:
            force.addPerBondParameter(parameter)
        restraint_index = force.addBond(particle1,particle2,self.bond_parameters)
        print restraint_index
        return force

    def _computeStandardStateCorrection(self):
        """
        Compute the standard state correction for the arbitrary restraint energy function.
        Returns
        -------
        DeltaG : float
            Computed standard-state correction in dimensionless units (kT)
        Notes
        -----
        Equivalent to the free energy of releasing restraints and confining into a box of standard state size.
        """
        initial_time = time.time()

        r_min = 0 * units.nanometers
        r_max = 100 * units.nanometers 
        
        # Create a System object containing two particles connected by the reference force
        system = openmm.System()
        system.addParticle(1.0 * units.amu)
        system.addParticle(1.0 * units.amu)
        force = self._createRestraintForce(0,1)
        system.addForce(force)

        # Create a Reference context to evaluate energies on the CPU
        integrator = openmm.VerletIntegrator(1.0 * units.femtoseconds)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)

        # Set default positions
        positions = units.Quantity(np.zeros([2,3]),units.nanometers)
        context.setPositions(positions)

        # Create a function to compute integrand as a function of interparticle separation.
        beta = self.beta

        def integrand(r):
            """
            Parameters
            ----------
            r : float
                Inter-particle separation in nanometers

            Returns
            -------
            dI : float
                Contribution to integrand(in nm^2)
            """
            positions[1,0] = r * units.nanometers
            context.setPositions(positions)
            state = context.getState(getEnergy=True)
            potential = state.getPotentialEnergy()
            dI = 4.0 * math.pi * r**2 * math.exp(-beta * potential)
            return dI

        (shell_volume,shell_volume_error) = scipy.integrate.quad(lambda r : integrand(r), r_min / units.nanometers, r_max / units.nanometers) * units.nanometers**3 #integrate shell volume
        print("shell_volume = %f nm^3" % (shell_volume / units.nanometers**3))
        
        # Compute standard-state volume for a single molecule in a box of size (1 L) / (avogadros number)
        liter = 1000.0 * units.centimeters**3 # one liter
        box_volume = liter / (units.AVOGADRO_CONSTANT_NA * units.mole) # standard state volume
        print("box_volume = %f nm^3" % (box_volume / units.nanometers**3))

        # Compute standard state correction for releasing shell restraints into standard-state box ( in units of kT)
        DeltaG = - math.log(box_volume / shell_volume)
        print ("Standard state correction: %.3f kT" % DeltaG)

        final_time = time.time()
        elapsed_time = final_time - initial_time
        print ("restraints: _computeStandardStateCorrection: %.3f s elapsed" % elapsed_time)

        # Return standard state correction (in kT)

    def getRestraintForce (self, mm=None):
        """
        Returns a new force object that imposes the receptor-ligand restraint
        
        Parameters
        ----------
        mm : simtk.openmm compliant interface
        
        Returns
        -------
        force: simtk.openmm.HarmonicBondForce
        """
        
        return self._createRestraintForce(self.restrained_receptor_atom,self.restrained_ligand_atom,mm=mm)

    def getRestrainedSystemCopy(self):
        """
        Returns a copy of the restrained system.

        Returns
        -------
        system: simtk.openmm.System
            a copy of the restrained system
        """

        system = copy.deepcopy(self.system)
        force = self.getRestraintForce()
        system.addForce(force)

        return system
    
    def getStandardStateCorrection(self):
        """
        Return the standard state correction.

        Returns
        -------
        correction: float
            the standard-state correction, in kT
        """
        return self.standard_state_correction

    def getReceptorRadiusOfGyration(self):
        """
        Returns the radius of gyration of the receptor

        Returns
        -------
        radius_of_gyration: simtk.unit.Quantity with units compatible with angstrom Radius of gyration of the receptor
        """
        return self.radius_of_gyration

    def _closestAtomToCentroid(self,sys,coordinates, indices = None, masses = None):
        """
        Identify the closest atom to the centroid of the given coordinate set
        Parameters
        ----------
        coordinates: units.Quantity of natoms * 3 with units compatible with nanometers coordinates of object to identify atom closes to centroid
        List of atoms indices for which closest atom to centroid is to be computed.
        masses: simtk.unit.Quantity of natoms with units compatible with amu
            Masses of particles used to weight distance calculation, if not None (default : None)
        
        Returns
        ------
        closest_atom: int
            Index of atom closest to centroid of specified atoms
        """
        
        if indices is not None:
            coordinates = coordinates[indices,:]

        # Get dimentionless coordinates
        x_unit = coordinates.unit
        x = coordinates / x_unit
        
        # Determine number of atoms
        natoms = x.shape[0]
        
        # Compute (natoms,1) array of normalized weights
        w = np.ones([natoms,1])
        if masses:
            w = masses / masses.unit # (natoms,)array
            w = np.reshape(w, (natoms,1)) #(natoms,1) array
            w /=w.sum()
        
        # Compute centroid (still in dimensionless units)
        centroid = (np.tile(w, (1,3)) * x).sum(0) # (3,) array
        
        # Compute distances from centroid
        distances = np.sqrt(((x - np.tile(centroid,(natoms,1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i
        
        #Determine closest atom
        closest_atom = int(np.argmin(distances))
        
        if indices is not None:
            closest_atom = indices[closest_atom]
        
        return closest_atom

    


class FlatBottomReceptorLigandRestraint(ReceptorLigandRestraint):
        
    energy_function = 'step(r-r0) * (K/2) * (r-r0)^2 / lambda_restraint' #flat-bottom restraint
    bond_parameter_names = ['K', 'r0'] # list of bond parameters that appear in energy function above
        
    def _determineBondParameters(self):
        """
        Determine bond parameters for CustomBondForce between protein and ligand
        
        Returns

        parameters(list) - list of parameters for CustomBondForce
        
        Note

        r0, the distance at which the harmonic restraint is imposed, is set at twice the robust estimate of the standard deviation (from mean absolute deviation) plus 5 A
        K, the spring constant, is set to 5.923 kcal/mol/A**2, for comparing with impact, K is set to 3.0 kcal/mol/A**2
        """

        x_unit = self.coordinates.unit

        # Get dimensionless receptor coordinates
        x = self.coordinates[self.receptor_atoms,:] / x_unit

        # Determine number of atoms
        natoms = x.shape[0]

        if (natoms > 3):
            # Compute median absolute distance to central atom
            # (Working in non-unit-bearing floats for speed)
            xref = np.reshape(x[self.restrained_receptor_atom,:], (1,3)) #(1,3) array
            distances = np.sqrt(((x - np.tile(xref, (natoms,1)))**2).sum(1)) # distances[i] is the distance from the centroid to particle i
            median_absolute_distance = np.median(abs(distances))

            # Convert back to unit-bearing quantity
            median_absolute_distance *= x_unit

            #Convert to estimator of standard deviation for normal distribution
            sigma = 1.4826 * median_absolute_distance

            # Calculate r0, which is a multiple of sigma plus 5 A
            #r0 = 2*sigma + 5.0 * units.angstroms
            r0 = 5.0 * units.angstroms
        else:
            DEFAULT_DISTANCE = 15.0 * units.angstroms
            print "WARNING: receptor only contains %d atoms; using default distance of %s" % (natoms, str(DEFAULT_DISTANCE))
            r0 = DEFAULT_DISTANCE
        print ("restraint distance r0 = %.1f A" % (r0 / units.angstroms))
        
        # Set spring constant
        # K = (2.0 * 0.0083144621 *50 * 298.0 * 100) * units.kilojoules_per_mole/units.nanometers**2
        K = 3.0 * units.kilocalories_per_mole / units.angstroms**2
        print ("K = %.1f kcal/mol/A^2" % (K / (units.kilocalories_per_mole / units.angstroms**2)))

        bond_parameters = [K, r0]
        
        return bond_parameters
    #edit end 3.16.15


    


    
