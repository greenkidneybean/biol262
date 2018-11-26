from simtk import openmm as openmm
from simtk.openmm import app
from simtk.unit import *

"""
only python3

module list
1) anaconda3/4.2.0  2) cuda/8.0.44


"""


def main():
  print("Reading the PSF file");

  # Read the PSF file
  psf = app.CharmmPsfFile('g1_25mm.psf');

  boxsize = 5 # Boxsize in nm
  psf.setBox(boxsize*nanometer, boxsize*nanometer, boxsize*nanometer);

  print("Reading the pdb file")
  pdb = app.PDBFile('g1_25mm.pdb');

  # Load the parameter set
  lib = 'toppar/'
  #params = app.CharmmParameterSet('toppar/top_all36_cgenff.rtf','toppar/par_all36_cgenff.prm','toppar/cgenff3.0.1/top_all36_cgenff.rtf','toppar/cgenff3.0.1/par_all36_cgenff.prm','g1_new.str','toppar_water_ions.str')
  params = app.CharmmParameterSet('toppar/cgenff3.0.1/top_all36_cgenff.rtf','toppar/cgenff3.0.1/par_all36_cgenff.prm','g1_new.str','toppar_water_ions.str')


  #platform = openmm.Platform.getPlatformByName('CUDA');
  platform = openmm.Platform.getPlatformByName('Reference');

  # Creating the system
  system = psf.createSystem(params, 
                            nonbondedMethod=app.PME, 
                            nonbondedCutoff=1.2*nanometer, 
                            switchDistance=1.0*nanometer, 
                            ewaldErrorTolerance= 0.0001, 
                            constraints=app.HBonds);

  # Thermostat @ 298 K
  system.addForce(openmm.AndersenThermostat(298*kelvin, 1/picosecond)) 
  # adding the barostat for now
  system.addForce(openmm.MonteCarloBarostat(1*bar, 298*kelvin));

  integrator = openmm.VerletIntegrator(0.001*picoseconds)

  simulation = app.Simulation(psf.topology, system, integrator, platform)
  simulation.context.setPositions(pdb.getPositions())
  simulation.minimizeEnergy(maxIterations = 500)

  #nsavcrd = 10000 # save coordinates every 10 ps
  #nstep  = 2000000 # write dcd files every 2 ps
  #nprint = 2000 # report every 2 ps
  nsavcrd = 10 # save coordinates every 10 ps
  nstep  = 200 # write dcd files every 2 ps
  nprint = 20 # report every 2 ps
  firstDcdStep = nsavcrd ;


  # Reporters
  dcd = app.DCDReporter('g1_new.dcd', nsavcrd)
  dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstDcdStep, nsavcrd)
  simulation.reporters.append(dcd)

  simulation.reporters.append(app.StateDataReporter('g1_new.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
  simulation.step(nstep)
  simulation.reporters.pop()
  simulation.reporters.pop()
  dcd._out.close()

if __name__ == '__main__' :
  main();
