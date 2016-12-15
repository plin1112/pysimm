Example 4: Creating a polyethylene chain using the random walk application
==========================================================================

### Importing pysimm modules/packages

This example code creates a polyethylene chain using the random walk application. To begin we need to import three modules/packages from pysimm: system, lmps, forcefield, and the random_walk application from the pysimm.apps package.

```
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
```

If you encounter an error **"ImportError: No module named pysimm"** make sure that the pysimm directory you cloned from github is in your PYTHONPATH. See installation instructions for further directions.

### Creating a polyethylene monomer

Following similar steps as Example 1, we'll build our polyethylene monomer starting from scratch. First we build an empty system and add a molecule.

```
s = system.System()
m = s.molecules.add(system.Molecule())
```

Then we'll create a reference to a Dreiding forcefield object, and set all of the interaction styles in our system to those in the force field.

```
f = forcefield.Dreiding()
s.pair_style = f.pair_style
s.bond_style = f.bond_style
s.angle_style = f.angle_style
s.dihedral_style = f.dihedral_style
s.improper_style = f.improper_style
```

Next we'll need references to ParticleType objects for the C_3 and H_ atom types in Dreiding.

```
dreiding_C_3 = s.particle_types.add(f.particle_types.get('C_3')[0].copy())
dreiding_H_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())
```

Then we'll build our monomer one atom at a time, starting with the carbon atoms.

```
c1 = s.particles.add(system.Particle(type=dreiding_C_3, x=0, y=0, z=0, charge=0, molecule=m))
c2 = s.add_particle_bonded_to(system.Particle(type=dreiding_C_3, charge=0, molecule=m), c1, f)

h1 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h2 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c1, f)
h3 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)
h4 = s.add_particle_bonded_to(system.Particle(type=dreiding_H_, charge=0, molecule=m), c2, f)
```

Lastly, we'll derive Gasteiger charges, set the simulation box, identify linker atoms, and optimize our structure using the Lennard-Jones potential.

```
s.apply_charges(f, charges='gasteiger')
s.set_box(padding=10)
c1.linker = 'head'
c2.linker = 'tail'
s.pair_style = 'lj'
lmps.quick_min(s, min_style='fire')
```

### Writing the monomer system to various file formats

The **system.System** class has various method to format our system data into different file types. Here we write xyz, yaml, lammps data, and chemdoodle json file formats.

```
s.write_xyz('pe_monomer.xyz')
s.write_yaml('pe_monomer.yaml')
s.write_lammps('pe_monomer.lmps')
s.write_chemdoodle_json('pe_monomer.json')
```

### Creating a polymer chain

The **random_walk** application requires a reference monomer, the number of repeat units in the desired chain, and the force field from which new force field types will be acquired. During polymerization, new force field terms will be determined automatically, and LAMMPS simulations will be performed to relax new polymer bonds. The random_walk function returns the newly created polymer **system.System** object.

`polymer = random_walk(s, 10, forcefield=f)`

### Writing the polymer system to various file formats

Like above, we can write our system to various file formats.

```
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')
```