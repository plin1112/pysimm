# ******************************************************************************
# pysimm.apps.equilibrate module
# ******************************************************************************
#
# 21-step equilibration algorithm written using pysimm tools
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from itertools import izip

from pysimm import system, lmps

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


def equil(s, **kwargs):
    """pysimm.apps.equilibrate.equil

    Runs a 21-step compression/decompression equilibration algorithm

    Args:
        s: pysimm.system.System object
        tmax: maximum temperature during equilibration
        pmax: maximum pressure during equilibration
        tfinal: desired final temperature of final system
        pfinal: desired final pressure of final system
        np: number of processors to use during equilibration simulations
        p_steps: list of pressures to use during equilibration (must match length of length_list)
        length_list: list of simulation durations to use during equilibration (must match length of p_steps)

    Returns:
        None
    """
    nb_cutoff = kwargs.get('nb_cutoff') or 12.0
    tmax = kwargs.get('tmax') or 1000
    pmax = kwargs.get('pmax') or 50000
    tfinal = kwargs.get('tfinal') or 300
    pfinal = kwargs.get('pfinal') or 1
    
    thermo = kwargs.get('thermo') if 'thermo' in kwargs else 1000
    thermo_style = kwargs.get('thermo_style') if 'thermo_style' in kwargs else 'custom step time temp density vol press etotal emol epair'

    nanohub = kwargs.get('nanohub')
    np = kwargs.get('np')

    p_list = kwargs.get('p_steps')

    if not p_list:
        p_list = [0.02*pmax, 0.6*pmax, pmax, 0.5*pmax, 0.1*pmax, 0.01*pmax, pfinal]

    length_list = kwargs.get('length_list')

    if not length_list:
        length_list = [100000, 100000, 100000, 100000, 100000, 100000, 100000]

    dump = kwargs.get('dump') or 1000
    dump_append = kwargs.get('dump_append') if kwargs.get('dump_append') is not None else True

    scale_v = kwargs.get('scale_v') if kwargs.get('scale_v') is not None else True

    settings = {
        'dump': dump, 'cutoff': nb_cutoff, 'dump_append': dump_append,
        'scale_v': scale_v, 'thermo': thermo, 'thermo_style': thermo_style
    }

    sim = lmps.Simulation(s, name='equil')
    sim.add_min(
        nanohub=nanohub, np=np, dump=dump, dump_append=dump_append, 
        thermo=thermo, thermo_style=thermo_style
    )
    sim.add_md(temp=tfinal, nanohub=nanohub, np=np, **settings)

    step = 0
    for p, l in izip(p_list, length_list):
        step += 1
        if l:
            sim.add_md(length=l, temp=tmax, **settings)
            sim.add_md(length=l, temp=tfinal, **settings)
            sim.add_md(length=l, ensemble='npt', temp=tfinal, pressure=p, **settings)

    sim.run(np=np, nanohub=nanohub)

    s.write_lammps('equil.lmps')
    s.write_xyz('equil.xyz')
