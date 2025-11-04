
# Prepare the system

## NOTE - AF2 structures were first renamed to be "start.pdb"
gmx pdb2gmx -f start.pdb -o complex.gro -ff oplsaa -water tip3p -ignh # convert to gro file and pick ff
gmx editconf -f complex.gro -o complex_newbox.gro -c -d 1.0 -bt cubic # prepare to solvate
gmx solvate -cp complex_newbox.gro -cs spc216.gro -o complex_solv.gro -p topol.top # solvate
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr
echo 13 | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral # Add ions
	
# Minimize

gmx grompp -f minim.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em # minimize

# Post Processing

echo 1 1 | gmx trjconv -f em.gro -s em.tpr -o centered.pdb -pbc nojump -center # convert to pdb
## NOTE - the output structure from GROMACS minimization, "centered.pdb", was renamed
