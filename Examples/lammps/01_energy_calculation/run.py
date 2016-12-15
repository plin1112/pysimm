from pysimm import system, lmps, forcefield

s = system.read_pubchem_smiles('CCCCCC')

s.apply_forcefield(forcefield.Dreiding())

initial_energy = lmps.energy(s)

lmps.quick_min(s)

final_energy_dreiding = lmps.energy(s, 'all')

s = system.read_pubchem_smiles('CCCCCC')

s.apply_forcefield(forcefield.Gaff2())

initial_energy = lmps.energy(s)

lmps.quick_min(s, min_style='sd')
lmps.quick_min(s, min_style='cg')

final_energy_gaff = lmps.energy(s, 'all')

print('Initial energy: {}\nFinal energy (Dreiding): {}\nFinal energy (Dreiding): {}'.format(initial_energy, final_energy_dreiding, final_energy_gaff))