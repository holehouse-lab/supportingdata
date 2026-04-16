# quick script to combine the VMD-derived frames into a single PDB ensemble
import mdtraj as md
x=md.load_pdb('frame1.pdb')
for i in range (2,7):
    x = x + md.load_pdb(f'frame{i}.pdb')
x.save_pdb('bound_ensemble.pdb')