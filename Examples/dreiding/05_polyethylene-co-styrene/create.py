from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
from pysimm.models.monomers.dreiding.pe import monomer as pe_monomer
from pysimm.models.monomers.dreiding.ps import monomer as ps_monomer

# get a pe monomer from the pysimm.models database
pe = pe_monomer()

# get a ps monomer from the pysimm.models database
ps = ps_monomer()

# we'll instantiate a Dreiding forcefield object for use later
f = forcefield.Dreiding()

# we will determine partial charges using the gasteiger algorithm
pe.apply_charges(f, charges='gasteiger')

# we will determine partial charges using the gasteiger algorithm
ps.apply_charges(f, charges='gasteiger')

# do a quick minimization of the monomer
lmps.quick_min(pe, min_style='fire')

# do a quick minimization of the monomer
lmps.quick_min(ps, min_style='fire')

# the buckingham potential isn't great at small distances, and therefore we use the LJ potential while growing the polymer
pe.pair_style = 'lj'
ps.pair_style = 'lj'

# run the copolymer random walk method with 10 total repeat units, using an alternating pattern
# here, settings is passed to the LAMMPS simulations during the polymerization, and we will use 2 processors
polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f, settings={'np': 2})

# write a few different file formats
polymer.write_xyz('polymer.xyz')
polymer.write_yaml('polymer.yaml')
polymer.write_lammps('polymer.lmps')
polymer.write_chemdoodle_json('polymer.json')