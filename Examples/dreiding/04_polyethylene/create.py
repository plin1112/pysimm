from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

s = system.System()
m = s.molecules.add(system.Molecule())
f = forcefield.Dreiding()

s.pair_style = f.pair_style
s.bond_style = f.bond_style
s.angle_style = f.angle_style
s.dihedral_style = f.dihedral_style
s.improper_style = f.improper_style

dreiding_C_3 = s.particle_types.add(f.particle_types.get('C_3')[0].copy())
dreiding_H_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())

c1 = s.particles.add(system.Particle(type=dreiding_C_3, x=0, y=0, z=0, charge=0, molecule=m))
c2 = s.add_particle_bonded_to(system.Particle(type=dreiding_C_3, charge=0, molecule=m), c1, f)

h1 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h2 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h3 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)
h4 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)

s.apply_charges(f, charges='gasteiger')

s.set_box(padding=10)

c1.linker = 'head'
c2.linker = 'tail'

s.pair_style = 'lj'

# do a quick minimization of the monomer
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
