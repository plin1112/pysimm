from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

# create an empty system
s = system.System()

# add a molecule to the system
m = s.molecules.add(system.Molecule())

# create a reference to the Gaff2 forcefield object
f = forcefield.Gaff2()

# set interaction styles of our system to those of the force field
s.pair_style = f.pair_style
s.bond_style = f.bond_style
s.angle_style = f.angle_style
s.dihedral_style = f.dihedral_style
s.improper_style = f.improper_style

# get references of particle types from force field
gaff_c3 = s.particle_types.add(f.particle_types.get('c3')[0].copy())
gaff_hc = s.particle_types.add(f.particle_types.get('hc')[0].copy())

# add particles one at a time starting with the carbons
c1 = s.particles.add(system.Particle(type=gaff_c3, x=0, y=0, z=0, charge=0, molecule=m))
c2 = s.add_particle_bonded_to(system.Particle(type=gaff_c3, charge=0, molecule=m), c1, f)

h1 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
h2 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c1, f)
h3 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)
h4 = s.add_particle_bonded_to(system.Particle(type=gaff_hc, charge=0, molecule=m), c2, f)

s.apply_charges(f, charges='gasteiger')

# set up the simulation box with a padding of 10 angstroms
s.set_box(padding=10)

# identify linker atoms
c1.linker = 'head'
c2.linker = 'tail'

# ensure we are optimizing using the LJ potential
s.pair_style = 'lj'

# optimize structure
lmps.quick_min(s, min_style='fire')

# write a few different file formats
s.write_xyz('pe_monomer.xyz')
s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')
s.write_chemdoodle_json('pe_monomer.json')

# run the random_walk polymerization method making a chain of 10 repeat units
# the forcefield object is supplied to get new forcefield types not in the monomer system
polymer = random_walk(s, 10, forcefield=f)

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
