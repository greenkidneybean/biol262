from simtk import openmm as openmm
from simtk.openmm import app
from simtk.unit import *
import mdtraj as md
import sys

"""
only python3

module list
1) anaconda3/4.2.0  2) cuda/8.0.44


"""


def main():
  print("Reading the PSF file");

  # Read the PSF file
  #psf = app.CharmmPsfFile('g1_25mm.psf');
  psf = app.CharmmPsfFile('holo_neutr.psf');

  boxsize = 4.895883 # Boxsize in nm
  psf.setBox(boxsize*nanometer, boxsize*nanometer, boxsize*nanometer);

  print("Reading the pdb file")
  #pdb = app.PDBFile('g1_25mm.pdb');
  #pdb  = app.PDBFile('md_298k_100ns.pdb')
  pdb  = app.PDBFile('holo_neutr.pdb')

  # Load the parameter set
  lib = 'toppar/'
  #params = app.CharmmParameterSet('toppar/top_all36_cgenff.rtf','toppar/par_all36_cgenff.prm','toppar/cgenff3.0.1/top_all36_cgenff.rtf','toppar/cgenff3.0.1/par_all36_cgenff.prm','g1_new.str','toppar_water_ions.str')
  #params = app.CharmmParameterSet('toppar/cgenff3.0.1/top_all36_cgenff.rtf','toppar/cgenff3.0.1/par_all36_cgenff.prm','g1_new.str','toppar_water_ions.str')
  params = app.CharmmParameterSet(
                            'toppar/cgenff3.0.1/top_all36_cgenff.rtf',
                            'toppar/cgenff3.0.1/par_all36_cgenff.prm',
                            'g1_new.str','oah_groups.str',
                            'toppar_water_ions.str'
  )


  platform = openmm.Platform.getPlatformByName('CUDA');

  isShortSim = False
  #isShortSim = True
  isPeriodic = True
  #isPeriodic = False
  #if not isPeriodic :
  #  isShortSim = True;

  # Creating the system
  if (isPeriodic) :
    print("PME is being used")
    system = psf.createSystem(params, 
                              nonbondedMethod=app.PME, 
                              nonbondedCutoff=1.2*nanometer, 
                              switchDistance=1.0*nanometer, 
                              ewaldErrorTolerance= 0.0001, 
                              constraints=app.HBonds
    );
  else :
    print("PME is not being used")
    system = psf.createSystem(params, 
    #                          nonbondedMethod=app.PME, 
                              nonbondedCutoff=1.2*nanometer, 
                              switchDistance=1.0*nanometer, 
                              ewaldErrorTolerance= 0.0001, 
                              constraints=app.HBonds
    );

  #for force in system.getForces():
  #  print(force, force.usesPeriodicBoundaryConditions())
  # Thermostat @ 298 K
  #system.addForce(openmm.AndersenThermostat(298*kelvin, 1/picosecond)) 
  # adding the barostat for now
  #system.addForce(openmm.MonteCarloBarostat(1*bar, 298*kelvin));

  # adding positional restriants
  if isPeriodic :
    force = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2");
  else :
    force = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)");
    #force = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2");

  force.addGlobalParameter("k", 0.1*kilocalories_per_mole/angstroms**2);
  force.addPerParticleParameter("x0");
  force.addPerParticleParameter("y0");
  force.addPerParticleParameter("z0");

  topology = pdb.getTopology();
  positions = pdb.getPositions();
  #for atom in  topology.atoms():
  #  print(atom)

  host = [];
  guest = [];
  for res in topology.residues():
    if res.name == 'OAH' :
      for atom in res.atoms():
        host.append(atom.index);
        force.addParticle(atom.index, positions[atom.index])
        #force.addParticle(atom.index, (0, 0, 0)*nanometers)
    if res.name == 'GOA' :
      for atom in res.atoms() :
        guest.append(atom.index);
  system.addForce(force)

  print("Does customExternalForce use periodic boundary condition : ", force.usesPeriodicBoundaryConditions())
  # adding restraint between the host and guest
  # this will be inside the umbrella sampling loop
  host_guest_centroid_force = openmm.CustomCentroidBondForce(2,"0.5*k*(distance(g1,g2)-d0)^2");
  host_guest_centroid_force.addGlobalParameter("k", 10.0*kilocalories_per_mole/angstroms**2);
  #d0_global_parameter_index = force.addGlobalParameter("d0", 2.0*angstroms);
  #d0_perBond_parameter_index = force.addPerBondParameter("d0", 2.0*angstroms);
  d0_perBond_parameter_index = host_guest_centroid_force.addPerBondParameter("d0");
  group1 = host_guest_centroid_force.addGroup(host)
  group2 = host_guest_centroid_force.addGroup(guest);
  host_guest_bond = host_guest_centroid_force.addBond([group1, group2] )
  host_guest_centroid_force.setBondParameters(host_guest_bond, [group1, group2], [20*angstroms]);
  system.addForce(host_guest_centroid_force);
  #sys.exit(0)

  """
  # Restrain along z axis
  # adding positional restriants
  if isPeriodic :
    z_force = openmm.CustomExternalForce("k*periodicdistance(x, x0)^2");
  else :
    z_force = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2)");
    #force = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2");

  z_force.addGlobalParameter("k", 0.1*kilocalories_per_mole/angstroms**2);
  z_force.addPerParticleParameter("x0");
  z_force.addPerParticleParameter("y0");

  for res in topology.residues():
    if res.name == 'GOA' :
      for atom in res.atoms():
        pos = list(positions[atom.index])
        print(pos[2])
        z_force.addParticle(atom.index, [pos[0], pos[1]])
  system.addForce(z_force)
  """


  # langevin integrator
  integrator = openmm.LangevinIntegrator(
                              298*kelvin,         # Temperature 
                              1.0/picoseconds,    # Friction coefficient
                              0.001*picoseconds   # Time step
  )
  #integrator = openmm.VerletIntegrator(0.001*picoseconds);
  simulation = app.Simulation(psf.topology, system, integrator, platform)
  simulation.context.setPositions(pdb.getPositions())
  
  #currentState = simulation.context.getState(getPositions=True)
  #pdbr = app.PDBReporter("pdbreport.pdb",1);
  #pdbr.report(simulation, currentState)
  
  print("Minimizing...")
  simulation.minimizeEnergy(maxIterations = 2000)
  print("Equilibrating...")
  simulation.context.setVelocitiesToTemperature(298*kelvin);
  #currentState = simulation.context.getState(getPositions=True)
  #pdbr = app.PDBReporter("pdbreport_after_min.pdb",1);
  #pdbr.report(simulation, currentState)
  #pdbr = app.PDBReporter("pdbreport_after_step1.pdb",1);


  
  nsavcrd = 10000 # save coordinates every 10 ps
  nprint = 2000 # report every 2 ps
  if isShortSim :
    nstep  = 20000 # run the simulation for 0.2ns 
  else :
    nstep  = 2000000 # run the simulation for 2ns 
  firstDcdStep = nsavcrd ;

  # Reporters
  dcdReportInterval = 10000 # save coordinates every 10 ps
  dcd = app.DCDReporter('umb_3.dcd', dcdReportInterval)
  dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstDcdStep, nsavcrd)

  stateReporter = app.StateDataReporter(
                      'umb_3.out', 
                      nprint, 
                      step=True, 
                      kineticEnergy=True, 
                      potentialEnergy=True, 
                      totalEnergy=True, 
                      temperature=True, 
                      volume=True, 
                      speed=True
  )

  simulation.reporters.append(dcd)
  simulation.reporters.append(stateReporter)
  #simulation.reporters.append(pdbr);
  #simulation.reporters.append(app.StateDataReporter('umb.out', nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))

  # Run the simulation 
  print("Simulation begins ...")
  #simulation.step(1)
  #sys.exit(0)
  for i in range(15) : 
    print("Simulation for umbrella : ", i)
    host_guest_centroid_force.setBondParameters(host_guest_bond, [group1, group2], [(5+3*i)*angstroms]);
    #force.setGlobalParameterDefaultValue(d0_global_parameter_index, (2+i)*angstroms)
    host_guest_centroid_force.updateParametersInContext(simulation.context);
    print("host guest bond parameters", host_guest_centroid_force.getBondParameters(host_guest_bond));
    #serialized = openmm.XmlSerializer.serialize(simulation.system)
    #of = open("serial_sim_"+str(i)+".xml",'w');
    #of.write(serialized);
    #of.close();
    #simulation.saveState("state_"+str(i)+".xml");
    simulation.step(nstep)

  simulation.reporters.pop()
  simulation.reporters.pop()
  dcd._out.close()

if __name__ == '__main__' :
  main();
