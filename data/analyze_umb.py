from __future__ import print_function
import mdtraj as md
import numpy as np

def analyze(dcdFileName, topologyFileName):
  traj = md.load_dcd(dcdFileName, topologyFileName);
  print(traj);
  print("shape : ", traj.xyz.shape)
  host_atoms = [atom.index for atom in traj.topology.atoms if atom.residue.name=='OAH']
  guest_atoms = [atom.index for atom in traj.topology.atoms if atom.residue.name=='GOA']
  for frame_id in range(traj.n_frames) :
    host_centroid = np.mean(traj.xyz[frame_id, host_atoms, :],axis=0);
    guest_centroid = np.mean(traj.xyz[frame_id, guest_atoms, :],axis=0);
    print(np.linalg.norm(host_centroid - guest_centroid))    
def main():
  analyze('umb_2.dcd', 'holo_neutr.psf')
if __name__ == '__main__' :
  main();
