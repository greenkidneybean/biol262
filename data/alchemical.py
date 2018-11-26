from simtk import openmm as openmm
from simtk.unit import *
from simtk.openmm import app
import numpy as np

"""
Alchemical host-guest


"""


def main():
  print("Reading the PSF file");

  # Read the PSF file
  psf = app.CharmmPsfFile('holo_neutr.psf');

  boxsize = 4.895883 # Boxsize in nm
  psf.setBox(boxsize*nanometer, boxsize*nanometer, boxsize*nanometer);

  print("Reading the pdb file")
  pdb  = app.PDBFile('holo_neutr.pdb')

  # Load the parameter set
  lib = 'toppar/'
  params = app.CharmmParameterSet(
                            'toppar/cgenff3.0.1/top_all36_cgenff.rtf',
                            'toppar/cgenff3.0.1/par_all36_cgenff.prm',
                            'g1_new.str','oah_groups.str',
                            'toppar_water_ions.str'
  )


  platform = openmm.Platform.getPlatformByName('CUDA');

  #isPeriodic = False
  #isShortSim = False
  isShortSim = True
  isPeriodic = False
  #if not isPeriodic :
  #  isShortSim = True;

  # Creating the system
  if (isPeriodic) :
    system = psf.createSystem(params, 
                              nonbondedMethod=app.PME, 
                              nonbondedCutoff=1.2*nanometer, 
                              switchDistance=1.0*nanometer, 
                              ewaldErrorTolerance= 0.0001, 
                              constraints=app.HBonds
    );
  else :
    system = psf.createSystem(params, 
    #                          nonbondedMethod=app.PME, 
                              nonbondedCutoff=1.2*nanometer, 
                              switchDistance=1.0*nanometer, 
                              ewaldErrorTolerance= 0.0001, 
                              constraints=app.HBonds
    );

  #Retrieve the nonbonded force 
  forces = { force.__class__.__name__  : for force in system.getForces()}:
  nbforce = forces['NonbondedForce']

  
if __name__ == '__main__' :
  main();  
