# ******************************************************************************
# pysimm.lmps module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
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

import shlex
import shutil
from subprocess import call, Popen, PIPE
from Queue import Queue, Empty
from threading import Thread
import os
import sys
import json
from random import randint
from time import strftime

from pysimm.system import read_lammps
from pysimm import error_print
from pysimm import warning_print
from pysimm import verbose_print
from pysimm import debug_print
from pysimm.utils import PysimmError, Item, ItemContainer

try:
    from Rappture.tools import getCommandOutput as RapptureExec
except ImportError:
    pass

try:
    from matplotlib import pyplot as plt
except ImportError:
    pass

LAMMPS_EXEC = os.environ.get('LAMMPS_EXEC')
verbose = False
templates = {}

def check_lmps_exec():
    if LAMMPS_EXEC is None:
        print 'you must set environment variable LAMMPS_EXEC'
        return False
    else:
        try:
            stdout, stderr = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                                   stdin=PIPE, stdout=PIPE,
                                   stderr=PIPE).communicate()
            if verbose:
                print 'using %s LAMMPS machine' % LAMMPS_EXEC
            return True
        except OSError:
            print 'LAMMPS is not configured properly for one reason or another'
            return False
            
check_lmps_exec()


class Thermo(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
    
    def write(self, inp):
        if self.freq is not None:
            inp = 'thermo {}\n'.format(self.freq) + inp
        if self.style:
            inp = 'thermo_style {}\n'.format(self.style) + inp
        return inp
        

class Ensemble(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        
    def write(self, inp):
        t_inp = ''
        if self.ensemble == 'nvt':
            t_inp += 'fix pysimm_md all nvt temp {} {} {}\n'.format(
                self.temp.start, self.temp.stop, self.temp.damp
            )
        elif self.ensemble == 'npt':
            iso = 'iso' if self.pressure.iso else 'aniso'
            t_inp += 'fix pysimm_md all npt temp {} {} {} {} {} {} {}\n'.format(
                self.temp.start, self.temp.stop, self.temp.damp,
                iso, self.pressure.start, self.pressure.stop, 
                self.pressure.damp
            )
        elif self.ensemble == 'nve':
            if self.limit is not None:
                t_inp += 'fix pysimm_md all nve/limit {}\n'.format(self.limit)
            else:
                t_inp += 'fix pysimm_md all nve\n'
        if t_inp:
            return t_inp + inp + 'unfix pysimm_md\n'
        else:
            return inp

class Velocity(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        
    def write(self, inp):
        t_inp = ''
        if self.new_v:
            t_inp += 'velocity all create {} {}\n'.format(
                self.temp.start, self.seed
            )
        elif self.scale_v:
            t_inp += 'velocity all scale {}\n'.format(
                self.temp.start
            )
        return t_inp + inp
            
            
class Shake(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        self.tol = 0.0001 if 'tol' not in kwargs else self.tol
        self.itr = 20 if 'itr' not in kwargs else self.itr
        self.stats = 10 if 'stats' not in kwargs else self.stats
        
    def write(self, inp):
        t_inp = ''
        if self.bond_types or self.angle_types:
            t_inp += 'fix pysimm_shake all shake {tol} {itr} {stats}'.format(
                tol=self.tol,
                itr=self.itr,
                stats=self.stats
            )
            if self.bond_types:
                t_inp += ' b {bt_tags}'.format(
                    bt_tags = ' '.join(map(str, [bt.tag for bt in self.bond_types]))
                )
            if self.angle_types:
                t_inp += ' a {at_tags}'.format(
                    at_tags = ' '.join(map(str, [at.tag for at in self.angle_types]))
                )
            t_inp +='\n'
            
            if self.velocity['new_v'] or self.velocity['scale_v']:
                t_inp += 'run 0\n'
                t_inp += 'velocity all scale {}\n'.format(self.temp.start)
            
            inp = t_inp + inp + 'unfix pysimm_shake\n'
        return inp
        
        
class Dump(Item):
    def __init__(self, **kwargs):
        Item.__init__(self, **kwargs)
        
    def write(self, inp):
        t_inp = ''
        name = '{}.lammpstrj'.format(self.name) or 'pysimm_dump.lammpstrj'
        if self.freq is not None:
            t_inp += 'dump pysimm_dump all custom {} {} id q x y z vx vy vz\n'.format(
                self.freq, name
            )
            if self.append:
                t_inp += 'dump_modify pysimm_dump append yes\n'
            inp = t_inp + inp + 'undump pysimm_dump\n'
        return inp
        
    
class MolecularDynamics(object):
    """pysimm.lmps.MolecularDynamics

    Template object to contain LAMMPS MD settings

    Attributes:
        timestep: timestep value to use during MD
        ensemble: 'nvt' or 'npt' or 'nve'
        limit: numerical value to use with nve when limiting particle displacement
        temp: temperature dictionary {start, stop, damp}
        pressure: pressure dictionary {start, stop, iso, damp}
        pressure_iso: True to require box dimensions change isotropically, False to allow dimensions to change independently (False by default)
        new_v: True to have LAMMPS generate new velocities
        seed: seed value for RNG (random by default)
        scale_v: True to scale velocities to given temperature default=False
        length: length of MD simulation in number of timesteps
        thermo: frequency to print thermodynamic data default=1000
        thermo_style: LAMMPS formatted input for thermo_style
        dump: frequency to dump trajectory
        dump_name: prefix of trajectory dump file
        dump_append: True to append to previous dump file is it exists
    """
    def __init__(self, **kwargs):

        self.timestep = kwargs.get('timestep', 1)
        self.ensemble = kwargs.get('ensemble', 'nvt')
        self.limit = kwargs.get('limit')
        self.length = kwargs.get('length', 2000)
        
        if isinstance(kwargs.get('thermo'), int):
            self.thermo = {'freq': kwargs.get('thermo')}
        else:
            self.thermo = kwargs.get('thermo', {
                'freq': 1000,
                'style': 'custom step time temp density vol press etotal emol epair'
            })
        
        self.dump = {
            'name': 'pysimm.dump.tmp',
            'append': False
        }
        self.dump.update(kwargs.get('dump', {}))
        
        if isinstance(kwargs.get('temp'), int) or isinstance(kwargs.get('temp'), float):
            t_value = kwargs.get('temp')
        else:
            t_value = 300
        
        self.temp = Item(**{
            'start': t_value,
            'stop': t_value,
            'damp': 100
        })
        
        if isinstance(kwargs.get('temp'), dict):
            self.temp.set(**kwargs.get('temp'))
        
        if isinstance(kwargs.get('pressure'), int) or isinstance(kwargs.get('pressure'), float):
            p_value = kwargs.get('pressure')
        else:
            p_value = 1
        
        self.pressure = Item(**{
            'start': p_value,
            'stop': p_value,
            'iso': False,
            'damp': 1000
        })
        if isinstance(kwargs.get('pressure'), dict):
            self.pressure.set(**kwargs.get('pressure', {}))
        
        self.seed = kwargs.get('seed', randint(10000, 99999))
        self.velocity = {'seed': self.seed}
        self.velocity.update(kwargs.get('velocity', {}))
        
        self.shake = kwargs.get('shake', {})
        self.input = ''

    def write(self):
        """pysimm.lmps.MolecularDynamics.write

        Create LAMMPS input for a molecular dynamics simulation.

        Args:
            None

        Returns:
            input string
        """
        self.input = ''

        self.input += 'run {}\n'.format(self.length)
        
        self.input = Shake(
            temp=self.temp, velocity=self.velocity, **self.shake
        ).write(self.input)
        
        self.input = Ensemble(
            ensemble=self.ensemble, temp=self.temp,
            pressure=self.pressure, limit=self.limit
        ).write(self.input)
        
        self.input = Dump(
            **self.dump
        ).write(self.input)
        
        self.input = 'timestep {}\n'.format(self.timestep) + self.input
        
        self.input = Velocity(
            temp=self.temp, **self.velocity
        ).write(self.input)
        
        self.input = Thermo(
            **self.thermo
        ).write(self.input)
        
        return self.input
        
        
class SteeredMolecularDynamics(MolecularDynamics):
    def __init__(self, **kwargs):
        MolecularDynamics.__init__(self, **kwargs)
        self.p1 = kwargs.get('p1')
        self.p2 = kwargs.get('p2')
        self.k = kwargs.get('k', 20.0)
        self.v = kwargs.get('v', 0.001)
        self.d = kwargs.get('d', 3.0)
    
    def write(self):
        """pysimm.lmps.SteeredMolecularDynamics.write

        Create LAMMPS input for a steered molecular dynamics simulation.

        Args:
            sim: pysimm.lmps.Simulation object reference

        Returns:
            input string
        """
        
        self.input += MolecularDynamics.write(self)
        
        self.input = 'fix steer p1 smd cvel {} {} couple p2 auto auto auto {}\n'.format(self.k, self.v, self.d) + self.input
        self.input = 'group p2 id {}\n'.format(self.p2.tag) + self.input
        self.input = 'group p1 id {}\n'.format(self.p1.tag) + self.input
        
        self.input += 'unfix steer\n'

        return self.input
        

class Minimization(object):
    """pysimm.lmps.Minimization

    Template object to contain LAMMPS energy minimization settings.

    Attributes:
        min_style: LAMMPS minimization style default='sd'
        etol: energy tolerance default=1e-3
        ftol: force tolerance default=1e-3
        maxiter: maximum iterations default=10000
        max eval: maximum force evaluations default=100000
        thermo: frequency to print thermodynamic data default=1000
        thermo_style: LAMMPS formatted input for thermo_style
        dump: frequency to dump trajectory
        dump_name: prefix of trajectory dump file
        dump_append: True to append to previous dump file is it exists
    """
    def __init__(self, **kwargs):

        self.min_style = kwargs.get('min_style', 'fire')
        self.dmax = kwargs.get('dmax')
        self.etol = kwargs.get('etol', 1.0e-3)
        self.ftol = kwargs.get('ftol', 1.0e-3)
        self.maxiter = kwargs.get('maxiter', 10000)
        self.maxeval = kwargs.get('maxeval', 100000)
        
        if isinstance(kwargs.get('thermo'), int):
            self.thermo = {'freq': kwargs.get('thermo')}
        else:
            self.thermo = kwargs.get('thermo', {
                'freq': 1000,
                'style': 'custom step time temp density vol press etotal emol epair'
            })
        
        self.dump = {
            'name': 'pysimm.dump.tmp',
            'append': False
        }
        self.dump.update(kwargs.get('dump', {}))
        
        if isinstance(kwargs.get('temp'), int) or isinstance(kwargs.get('temp'), float):
            t_value = kwargs.get('temp')
        else:
            t_value = 300
        
        self.temp = Item(**{
            'start': t_value,
            'stop': t_value,
            'damp': 100
        })
        
        if isinstance(kwargs.get('temp'), dict):
            self.temp.set(**kwargs.get('temp'))
        
        self.seed = kwargs.get('seed') or randint(10000, 99999)
        self.velocity = {'seed': self.seed}
        self.velocity.update(kwargs.get('velocity', {}))

        self.input = ''

    def write(self):
        """pysimm.lmps.Minimization.write

        Create LAMMPS input for an energy minimization simulation.

        Args:
            sim: pysimm.lmps.Simulation object reference

        Returns:
            input string
        """
        self.input = ''

        self.input += 'min_style %s\n' % self.min_style
        if self.dmax:
            self.input += 'min_modify dmax %s\n' % self.dmax
        self.input += ('minimize %s %s %s %s\n' % (self.etol, self.ftol,
                                                   self.maxiter, self.maxeval))
            
        self.input = Velocity(
            temp=self.temp, **self.velocity
        ).write(self.input)
        
        self.input = Thermo(
            **self.thermo
        ).write(self.input)

        return self.input


class CustomInput(object):
    """pysimm.lmps.CustomInput

    Template object to contain custom LAMMPS input.

    Attributes:
        custom_input: custom input string
    """
    def __init__(self, custom_input):
        self.input = '{}\n'.format(custom_input)

    def write(self):
        """pysimm.lmps.CustomInput.write

        Create LAMMPS input for a custom simulation.

        Args:
            sim: pysimm.lmps.Simulation object reference

        Returns:
            input string
        """
        return self.input


class Simulation(object):
    """pysimm.lmps.Simulation

    Organizational object for LAMMPS simulation. Should contain combination of
    pysimm.lmps.MolecularDynamics, pysimm.lmps.Minimization, and/or pysimm.lmps.CustomInput object.

    Attributes:
        atom_style: LAMMPS atom_style default=full
        kspace_style: LAMMPS kspace style default='pppm 1e-4'
        units: LAMMPS set of units to use default=real
        special_bonds: LAMMPS special bonds input
        nonbond_mixing: type of mixing rule for nonbonded interactions default=arithmetic
        cutoff: cutoff for nonbonded interactions default=12
        print_to_screen: True to have LAMMPS output printed to stdout
        name: name id for simulations
        log: prefix for LAMMPS log file
        write: file name to write final LAMMPS data file default=None
    """
    def __init__(self, s, **kwargs):

        self.system = s

        self.atom_style = kwargs.get('atom_style', 'full')
        self.kspace_style = kwargs.get('kspace_style', 'pppm 1e-4')
        self.units = kwargs.get('units', 'real')
        self.special_bonds = kwargs.get('special_bonds')
        self.nonbond_mixing = kwargs.get('nonbond_mixing')
        self.cutoff = kwargs.get('cutoff', 12.0)
        self.inner_cutoff = kwargs.get('inner_cutoff', 8.0)

        self.print_to_screen = kwargs.get('print_to_screen', False)
        self.name = kwargs.get('name')
        self.log = kwargs.get('log')
        self.write = kwargs.get('write')

        self.input = ''
        self.custom = kwargs.get('custom')

        self.sim = kwargs.get('sim', [])
        
        self.results = []

    def add_md(self, template=None, **kwargs):
        """pysimm.lmps.Simulation.add_md

        Add pysimm.lmps.MolecularDyanmics template to simulation

        Args:
            template: pysimm.lmps.MolecularDynamics object reference
            **kwargs: if template is None these are passed to
            pysimm.lmps.MolecularDynamics constructor to create new template
        """
        if template is None:
            self.sim.append(MolecularDynamics(**kwargs))
        elif isinstance(template, MolecularDynamics):
            self.sim.append(template)
        else:
            error_print('you must add an object of type MolecularDynamics to Simulation')
            
        return self

    def add_min(self, template=None, **kwargs):
        """pysimm.lmps.Simulation.add_min

        Add pysimm.lmps.Minimization template to simulation

        Args:
            template: pysimm.lmps.Minimization object reference
            **kwargs: if template is None these are passed to
            pysimm.lmps.Minimization constructor to create new template
        """
        if template is None:
            self.sim.append(Minimization(**kwargs))
        elif isinstance(template, Minimization):
            self.sim.append(template)
        else:
            error_print('you must add an object of type Minimization to Simulation')
        
        return self

    def add_custom(self, custom=''):
        """pysimm.lmps.Simulation.add_custom

        Add custom input string to simulation

        Args:
            custom: custom LAMMPS input string to add to Simulation
        """
        self.sim.append(CustomInput(custom))
        
        return self

    def write_input(self, init=True):
        """pysimm.lmps.Simulation.write_input

        Creates LAMMPS input string including initialization and input from templates/custom input

        Args:
            None

        Returns:
            None
        """
        self.input = ''

        if init:
            self.input += write_init(
                self.system, atom_style=self.atom_style, kspace_style=self.kspace_style,
                special_bonds=self.special_bonds, units=self.units,
                nonbond_mixing=self.nonbond_mixing,
                nb_cut=self.cutoff, inner_cutoff=self.inner_cutoff
            )

        if self.log:
            self.input += 'log %s append\n' % self.log
        elif self.name:
            self.input += 'log %s.log append\n' % '_'.join(self.name.split())
        else:
            self.input += 'log log.lammps append\n'

        for template in self.sim:
            self.input += template.write()
            
        self.input += 'write_dump all custom pysimm.dump.tmp id q x y z vx vy vz\n'
        self.input += 'quit\n'
        return self

    def run(self, np=None, kokkos_gpus=None, gpu_gpus=None, nanohub=None, rewrite=True, init=True, write_input=True):
        """pysimm.lmps.Simulation.run

        Begin LAMMPS simulation.

        Args:
            np: number of threads to use (serial by default) default=None
            nanohub: dictionary containing nanohub resource information default=None
            rewrite: True to rewrite input before running default=True
            init: True to write initialization part of LAMMPS input script (set to False if using complete custom input)
            write_input: write lammps input to pysimm.sim.in
        """
        if self.custom:
            rewrite = False
            self.input += '\nwrite_data pysimm_md.lmps\n'
        if rewrite:
            self.write_input(init)
        if isinstance(write_input, str):
            with file(write_input, 'w') as f:
                f.write(self.input)
        elif write_input:
            with file('pysimm.sim.in', 'w') as f:
                f.write(self.input)
        try:
            call_lammps(self, np, kokkos_gpus, gpu_gpus, nanohub)
        except OSError as ose:
            raise PysimmError('There was a problem calling LAMMPS with mpiexec'), None, sys.exc_info()[2]
        except IOError as ioe:
            if check_lmps_exec():
                raise PysimmError('There was a problem running LAMMPS. The process started but did not finish successfully. Check the log file, or rerun the simulation with print_to_screen=True to debug issue from LAMMPS output'), None, sys.exc_info()[2]
            else:
                raise PysimmError('There was a problem running LAMMPS. LAMMPS is not configured properly. Make sure the LAMMPS_EXEC environment variable is set to the correct LAMMPS executable path. The current path is set to:\n\n{}'.format(LAMMPS_EXEC)), None, sys.exc_info()[2]
                
        try:
            if self.log:
                log_name = self.log
            elif self.name:
                log_name = '{}.log'.format('_'.join(self.name.split()))
            else:
                log_name = 'log.lammps'
            r = read_log(log_name)
            self.results = r
        except Exception as e:
            print('Failed reading log file')
            print(e)
            
    def plot_results(self):
        steps = []
        last_step = 0
        for sim in self.results:
            sts = sim['Step']
            for st in sts:
                steps.append(st+last_step)
            last_step += sts[-1]
        for k in self.results[-1]:
            if k != 'Step':
                data = [sm[k] for sm in self.results]
                data = [d for sublist in data for d in sublist]
                plt.plot(steps, data)
                plt.xlabel('Step')
                plt.ylabel(k)
                plt.show()


def enqueue_output(out, queue):
    """pysimm.lmps.enqueue_output

    Helps queue output for printing to screen during simulation.
    """
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()


def call_lammps(simulation, np, kokkos_gpus, gpu_gpus, nanohub):
    """pysimm.lmps.call_lammps

    Wrapper to call LAMMPS using executable name defined in pysimm.lmps module.

    Args:
        simulation: pysimm.lmps.Simulation object reference
        np: number of threads to use
        nanohub: dictionary containing nanohub resource information default=None

    Returns:
        None
    """
    if nanohub:
        with file('temp.in', 'w') as f:
            f.write(simulation.input)
        if simulation.name:
            print('%s: sending %s simulation to computer cluster at nanoHUB' % (strftime('%H:%M:%S'), simulation.name))
        else:
            print('%s: sending simulation to computer cluster at nanoHUB' % strftime('%H:%M:%S'))
        sys.stdout.flush()
        cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
               'lammps-09Dec14-parallel -e both -l none -i temp.in'
               % (nanohub.get('cores'), nanohub.get('walltime')))
        cmd = shlex.split(cmd)
        exit_status, stdo, stde = RapptureExec(cmd)
    else:
        if simulation.name:
            print('%s: starting %s LAMMPS simulation'
                  % (strftime('%H:%M:%S'), simulation.name))
        else:
            print('%s: starting LAMMPS simulation'
                  % strftime('%H:%M:%S'))
        if kokkos_gpus:
            if np:
                p = Popen(['mpiexec', '-np', str(np),
                           LAMMPS_EXEC, '-e', 'both', '-l', 'none',
                           '-k', 'on g', str(kokkos_gpus), '-sf' 'kk', '-pk', 'kokkos'],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
            else:
                p = Popen(['mpiexec',
                           LAMMPS_EXEC, '-e', 'both', '-l', 'none',
                           '-k', 'on g', str(kokkos_gpus), '-sf' 'kk', '-pk', 'kokkos'],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
        elif gpu_gpus:
            if np:
                p = Popen(['mpiexec', '-np', str(np),
                           LAMMPS_EXEC, '-e', 'both', '-l', 'none',
                           '-sf', 'gpu', '-pk', 'gpu', str(gpu_gpus)],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
            else:
                p = Popen(['mpiexec',
                           LAMMPS_EXEC, '-e', 'both', '-l', 'none',
                           '-sf', 'gpu', '-pk', 'gpu', str(gpu_gpus)],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
        elif np:
            p = Popen(['mpiexec', '-np', str(np),
                       LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=PIPE, stdout=PIPE, stderr=PIPE)
        else:
            p = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.stdin.write(simulation.input)
        q = Queue()
        t = Thread(target=enqueue_output, args=(p.stdout, q))
        t.daemon = True
        t.start()

        while t.isAlive() or not q.empty():
            try:
                line = q.get_nowait()
            except Empty:
                pass
            else:
                if simulation.print_to_screen:
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    
    simulation.system.read_lammps_dump('pysimm.dump.tmp')

    try:
        os.remove('temp.lmps')
    except OSError as e:
        print e
        
    try:
        os.remove('pysimm.dump.tmp')
        if simulation.name:
            print('%s: %s simulation using LAMMPS successful'
                  % (strftime('%H:%M:%S'), simulation.name))
        else:
            print('%s: molecular dynamics using LAMMPS successful'
                  % (strftime('%H:%M:%S')))
    except OSError as e:
        if simulation.name:
            raise PysimmError('%s simulation using LAMMPS UNsuccessful' % simulation.name)
        else:
            raise PysimmError('molecular dynamics using LAMMPS UNsuccessful')


def quick_md(s, np=None, nanohub=None, **kwargs):
    """pysimm.lmps.quick_md

    Convenience function to call an individual MD simulation. kwargs are passed to MD constructor

    Args:
        s: system to perform simulation on
        np: number of threads to use
        nanohub: dictionary containing nanohub resource information default=None

    Returns:
        None
    """
    sim = Simulation(s, **kwargs)
    sim.add_md(**kwargs)
    sim.run(np, nanohub)
    
    return sim


def quick_min(s, np=None, nanohub=None, **kwargs):
    """pysimm.lmps.quick_min

    Convenience function to call an individual energy minimization simulation. kwargs are passed to min constructor

    Args:
        s: system to perform simulation on
        np: number of threads to use
        nanohub: dictionary containing nanohub resource information default=None

    Returns:
        None
    """
    sim = Simulation(s, **kwargs)
    sim.add_min(**kwargs)
    sim.run(np, nanohub)
    
    return sim
    
    
def energy(s, all=False, np=None, **kwargs):
    """pysimm.lmps.energy

    Convenience function to calculate energy of a given System object.

    Args:
        s: system to calculate energy
        all: returns decomposition of energy if True (default: False)
        np: number of threads to use for simulation

    Returns:
        total energy or disctionary of energy components
    """
    sim = Simulation(s, log='pysimm_calc.tmp.log', **kwargs)
    sim.add_md(length=0, thermo={'freq': 1, 'style': 'custom step etotal epair emol evdwl ecoul ebond eangle edihed eimp'}, **kwargs)
    sim.run(np)
    with file('pysimm_calc.tmp.log') as f:
        line = f.next()
        while line.split()[0] != 'Step':
            line = f.next()
        line = f.next()
        step, etotal, epair, emol, evdwl, ecoul, ebond, eangle, edihed, eimp = map(float, line.split())
    try:
        os.remove('pysimm_calc.tmp.log')
    except:
        error_print('error likely occurred during simulation')
    if all:
        return {
                'step': int(step),
                'etotal': etotal,
                'epair': epair,
                'emol': emol,
                'evdwl': evdwl,
                'ecoul': ecoul,
                'ebond': ebond,
                'eangle': eangle,
                'edihed': edihed,
                'eimp': eimp
               }
    else:
        return etotal
        
        
def read_log(fname):
    sims = []
    with open(fname) as f:
        for line in f:
            if 'Step' in line:
                sims.append({})
                headers = line.split()
                for h in headers:
                    sims[-1][h] = []
                line = f.next()
                while 'Loop' not in line:
                    for d, h in zip(map(float, line.split()), headers):
                        sims[-1][h].append(d)
                    line = f.next()
    return sims


def md(s, template=None, **kwargs):
    """pysimm.lmps.md

    Convenience function for performing LAMMPS MD

    *** WILL BE DEPRECATED - USE QUICK_MD INSTEAD ***
    """
    global LAMMPS_EXEC

    if template:
        template.update(kwargs)
        kwargs = template

    name = kwargs.get('name') or False
    log = kwargs.get('log')
    write = kwargs.get('write') or False
    print_to_screen = kwargs.get('print_to_screen') if kwargs.get(
        'print_to_screen') is not None else False
    special_bonds = kwargs.get('special_bonds') or 'amber'
    cutoff = kwargs.get('cutoff') or 12.0
    timestep = kwargs.get('timestep') or 1
    ensemble = kwargs.get('ensemble') or 'nvt'
    temp = kwargs.get('temp')
    pressure = kwargs.get('pressure') or 1.
    new_v = kwargs.get('new_v')
    seed = kwargs.get('seed') or randint(10000, 99999)
    scale_v = kwargs.get('scale_v')
    length = kwargs.get('length') or 2000
    thermo = kwargs.get('thermo') or 1000
    thermo_style = kwargs.get('thermo_style')
    nonbond_mixing = kwargs.get('nonbond_mixing')
    kspace_style = kwargs.get('kspace_style') or 'pppm 1e-4'

    nanohub = kwargs.get('nanohub') or {}

    pbs = kwargs.get('pbs')
    np = kwargs.get('np')
    kokkos = kwargs.get('kokkos')

    dump = kwargs.get('dump') or False
    dump_name = kwargs.get('dump_name')
    dump_append = kwargs.get('dump_append')

    if temp is None:
        t_start = kwargs.get('t_start')
        t_stop = kwargs.get('t_stop')
        if t_start is None:
            t_start = 1000.
        if t_stop is None:
            t_stop = t_start
    else:
        t_start = temp
        t_stop = temp

    command = write_init(s, nb_cut=cutoff, special_bonds=special_bonds,
                         nonbond_mixing=nonbond_mixing, kspace_style=kspace_style)
    if log:
        command += 'log %s append\n' % log
    elif name:
        command += 'log %s.log append\n' % '_'.join(name.split())
    else:
        command += 'log log.lammps append\n'
    if thermo:
        command += 'thermo %s\n' % int(thermo)
    if thermo_style:
        command += 'thermo_style %s\n' % thermo_style
    command += 'timestep %s\n' % timestep
    if ensemble == 'nvt':
        command += 'fix 1 all %s temp %s %s 100\n' % (ensemble, t_start, t_stop)
    elif ensemble == 'npt':
        command += ('fix 1 all %s temp %s %s 100 iso %s %s 100\n'
                    % (ensemble, t_start, t_stop, pressure, pressure))
    if new_v:
        command += 'velocity all create %s %s\n' % (t_start, seed)
    elif scale_v:
        command += 'velocity all scale %s\n' % t_start

    if dump:
        if dump_name:
            command += ('dump pysimm_dump all atom %s %s.lammpstrj\n'
                        % (dump, dump_name))
        elif name:
            command += ('dump pysimm_dump all atom %s %s.lammpstrj\n'
                        % (dump, '_'.join(name.split())))
        else:
            command += ('dump pysimm_dump all atom %s pysimm_dump.lammpstrj\n'
                        % dump)
        if dump_append:
            command += 'dump_modify pysimm_dump append yes\n'
    command += 'run %s\n' % length
    command += 'unfix 1\n'
    if write:
        command += 'write_data %s\n' % write
    else:
        command += 'write_data pysimm_md.lmps\n'

    with open('temp.in', 'w') as f:
        f.write(command)

    if name:
        print('%s: starting %s simulation using LAMMPS'
              % (strftime('%H:%M:%S'), name))
    else:
        print('%s: starting molecular dynamics using LAMMPS'
              % strftime('%H:%M:%S'))

    if nanohub:
        if name:
            print('%s: sending %s simulation to computer cluster' % (strftime('%H:%M:%S'), name))
        sys.stdout.flush()
        cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
               'lammps-09Dec14-parallel -e both -l none -i temp.in'
               % (nanohub.get('cores'), nanohub.get('walltime')))
        cmd = shlex.split(cmd)
        exit_status, stdo, stde = RapptureExec(cmd)
    elif pbs:
        call('mpiexec %s -e both -l log' % LAMMPS_EXEC, shell=True,
             stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
    else:
        if np:
            p = Popen(['mpiexec', '-np', str(np),
                       LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
        elif kokkos:
            p = Popen([LAMMPS_EXEC, '-k', 'on', '-sf', 'kk', '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
        else:
            p = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)

        while True:
            out = p.stdout.read(1)
            if out == '' and p.poll() is not None:
                break
            if out != '' and print_to_screen:
                sys.stdout.write(out)
                sys.stdout.flush()

    if write:
        n = read_lammps(write, quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    else:
        n = read_lammps('pysimm_md.lmps', quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    for p in n.particles:
        p_ = s.particles[p.tag]
        p_.x = p.x
        p_.y = p.y
        p_.z = p.z
        p_.vx = p.vx
        p_.vy = p.vy
        p_.vz = p.vz
    s.dim = n.dim
    os.remove('temp.in')

    try:
        os.remove('temp.lmps')
    except OSError as e:
        print e

    if not write:
        try:
            os.remove('pysimm_md.lmps')
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: molecular dynamics using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        except OSError:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: molecular dynamics using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))
            return False

    else:
        if os.path.isfile(write):
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: molecular dynamics using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        else:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: molecular dynamics using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))
            return False


def minimize(s, template=None, **kwargs):
    """pysimm.lmps.minimize

    Convenience function for performing LAMMPS energy minimization

    *** WILL BE DEPRECATED - USE QUICK_MIN INSTEAD ***
    """
    global LAMMPS_EXEC

    if template:
        template.update(kwargs)
        kwargs = template

    name = kwargs.get('name') or False
    log = kwargs.get('log') or 'log.lammps'
    write = kwargs.get('write') or False
    print_to_screen = kwargs.get('print_to_screen') if kwargs.get(
        'print_to_screen') is not None else False
    special_bonds = kwargs.get('special_bonds') or 'amber'
    cutoff = kwargs.get('cutoff') or 12.0
    min_style = kwargs.get('min_style')
    fire_etol = kwargs.get('sd_etol') or 1.0e-3
    fire_ftol = kwargs.get('sd_ftol') or 1.0e-3
    fire_maxiter = kwargs.get('sd_maxiter') or 10000
    fire_maxeval = kwargs.get('sd_maxeval') or 100000
    sd_etol = kwargs.get('sd_etol') or 1.0e-3
    sd_ftol = kwargs.get('sd_ftol') or 1.0e-3
    sd_maxiter = kwargs.get('sd_maxiter') or 10000
    sd_maxeval = kwargs.get('sd_maxeval') or 100000
    cg_etol = kwargs.get('cg_etol') or 1.0e-6
    cg_ftol = kwargs.get('cg_ftol') or 1.0e-6
    cg_maxiter = kwargs.get('cg_maxiter') or 10000
    cg_maxeval = kwargs.get('cg_maxeval') or 100000
    thermo = kwargs.get('thermo') or 1000
    thermo_style = kwargs.get('thermo_style')
    nonbond_mixing = kwargs.get('nonbond_mixing')
    kspace_style = kwargs.get('kspace_style') or 'pppm 1e-4'

    nanohub = kwargs.get('nanohub') or {}

    pbs = kwargs.get('pbs')
    np = kwargs.get('np')

    command = write_init(s, nb_cut=cutoff, special_bonds=special_bonds,
                         nonbond_mixing=nonbond_mixing, kspace_style=kspace_style)
    if log:
        command += 'log %s append\n' % log
    elif name:
        command += 'log %s.log append\n' % '_'.join(name.split())
    else:
        command += 'log log.lammps append\n'

    if thermo:
        command += 'thermo %s\n' % int(thermo)
    if thermo_style:
        command += 'thermo_style %s\n' % thermo_style

    if not min_style or min_style == 'sd':
        command += 'min_style sd\n'
        command += ('minimize %s %s %s %s\n'
                    % (sd_etol, sd_ftol, sd_maxiter, sd_maxeval))

        command += 'min_style cg\n'
        command += ('minimize %s %s %s %s\n'
                    % (cg_etol, cg_ftol, cg_maxiter, cg_maxeval))

    elif min_style == 'fire':
        command += 'timestep 1\n'
        command += 'min_style fire\n'
        command += ('minimize %s %s %s %s\n'
                    % (fire_etol, fire_ftol, fire_maxiter, fire_maxeval))

    if write:
        command += 'write_data %s\n' % write
    else:
        command += 'write_data pysimm_min.lmps\n'

    with open('temp.in', 'w') as f:
        f.write(command)

    if name:
        print('%s: starting %s simulation using LAMMPS'
              % (strftime('%H:%M:%S'), name))
    else:
        print('%s: starting minimization using LAMMPS'
              % strftime('%H:%M:%S'))

    if nanohub:
        if name:
            print('%s: sending %s simulation to computer cluster' % (strftime('%H:%M:%S'), name))
        sys.stdout.flush()
        cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
               'lammps-09Dec14-parallel -e both -l none -i temp.in'
               % (nanohub.get('cores'), nanohub.get('walltime')))
        cmd = shlex.split(cmd)
        exit_status, stdo, stde = RapptureExec(cmd)
    if pbs:
        call('mpiexec %s -e both -l log' % LAMMPS_EXEC, shell=True,
                 stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
    else:
        if np:
            p = Popen(['mpiexec', '-np', str(np),
                       LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
        else:
            p = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)

        while True:
            out = p.stdout.read(1)
            if out == '' and p.poll() is not None:
                break
            if out != '' and print_to_screen:
                sys.stdout.write(out)
                sys.stdout.flush()

    if write:
        n = read_lammps(write, quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    else:
        n = read_lammps('pysimm_min.lmps', quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    for p in n.particles:
        s.particles[p.tag].x = p.x
        s.particles[p.tag].y = p.y
        s.particles[p.tag].z = p.z
    os.remove('temp.in')

    try:
        os.remove('temp.lmps')
    except OSError as e:
        print e

    if not write:
        try:
            os.remove('pysimm_min.lmps')
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: minimization using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        except OSError:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: minimization using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))
            return False

    else:
        if os.path.isfile(write):
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: minimization using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        else:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: minimization using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))
            return False


def relax(s, template=None, **kwargs):
    """pysimm.lmps.md

    Convenience function for performing LAMMPS MD

    *** WILL BE DEPRECATED - USE QUICK_MD INSTEAD ***
    """
    global LAMMPS_EXEC

    if template:
        template.update(kwargs)
        kwargs = template

    name = kwargs.get('name') or False
    log = kwargs.get('log') or 'log.lammps'
    write = kwargs.get('write') or False
    print_to_screen = kwargs.get('print_to_screen') if kwargs.get(
        'print_to_screen') is not None else False
    special_bonds = kwargs.get('special_bonds') or 'amber'
    cutoff = kwargs.get('cutoff') or 12.0
    xmax = kwargs.get('xmax') or 0.1
    temp = kwargs.get('temp')
    new_v = kwargs.get('new_v')
    seed = kwargs.get('seed') or randint(10000, 99999)
    scale_v = kwargs.get('scale_v')
    length = kwargs.get('length') or 2000
    thermo = kwargs.get('thermo') or 1000
    thermo_style = kwargs.get('thermo_style')
    nonbond_mixing = kwargs.get('nonbond_mixing')
    kspace_style = kwargs.get('kspace_style') or 'pppm 1e-4'

    nanohub = kwargs.get('nanohub') or {}

    pbs = kwargs.get('pbs')
    np = kwargs.get('np')

    dump = kwargs.get('dump') or False
    dump_name = kwargs.get('dump_name')
    dump_append = kwargs.get('dump_append')

    if temp is None:
        t_start = kwargs.get('t_start')
        t_stop = kwargs.get('t_stop')
        if t_start is None:
            t_start = 1000.
        if t_stop is None:
            t_stop = t_start
    else:
        t_start = temp
        t_stop = temp

    command = write_init(s, nb_cut=cutoff, special_bonds=special_bonds,
                         nonbond_mixing=nonbond_mixing, kspace_style=kspace_style)
    if log:
        command += 'log %s append\n' % log
    elif isinstance(name, basestring):
        command += 'log %s.log append\n' % '_'.join(name.split())
    else:
        command += 'log log.lammps append\n'

    if thermo:
        command += 'thermo %s\n' % int(thermo)
    if thermo_style:
        command += 'thermo_style %s\n' % thermo_style

    if dump:
        if dump_name:
            command += ('dump pysimm_dump all atom %s %s.lammpstrj\n'
                        % (dump, dump_name))
        elif name:
            command += ('dump pysimm_dump all atom %s %s.lammpstrj\n'
                        % (dump, name))
        else:
            command += ('dump pysimm_dump all atom %s pysimm_dump.lammpstrj\n'
                        % dump)
        if dump_append:
            command += 'dump_modify pysimm_dump append yes\n'

    if new_v:
        command += 'velocity all create %s %s\n' % (t_start, seed)
    elif scale_v:
        command += 'velocity all scale %s\n' % t_start

    command += 'fix 1 all nve/limit %s\n' % xmax
    command += 'run %s\n' % length
    command += 'unfix 1\n'

    if write:
        command += 'write_data %s\n' % write
    else:
        command += 'write_data pysimm_relax.lmps\n'

    with open('temp.in', 'w') as f:
        f.write(command)

    if name:
        print('%s: starting %s simulation using LAMMPS'
              % (strftime('%H:%M:%S'), name))
    else:
        print('%s: starting nve/limit relaxation using LAMMPS'
              % strftime('%H:%M:%S'))

    if nanohub:
        if name:
            print('%s: sending %s simulation to computer cluster' % (strftime('%H:%M:%S'), name))
        sys.stdout.flush()
        cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
               'lammps-09Dec14-parallel -e both -l none -i temp.in'
               % (nanohub.get('cores'), nanohub.get('walltime')))
        cmd = shlex.split(cmd)
        exit_status, stdo, stde = RapptureExec(cmd)
    if pbs:
        call('mpiexec %s -e both -l log' % LAMMPS_EXEC, shell=True,
             stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
    else:
        if np:
            p = Popen(['mpiexec', '-np', str(np),
                       LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)
        else:
            p = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'],
                      stdin=open('temp.in'), stdout=PIPE, stderr=PIPE)

        while True:
            out = p.stdout.read(1)
            if out == '' and p.poll() is not None:
                break
            if out != '' and print_to_screen:
                sys.stdout.write(out)
                sys.stdout.flush()

    if write:
        n = read_lammps(write, quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    else:
        n = read_lammps('pysimm_relax.lmps', quiet=True,
                        pair_style=s.pair_style,
                        bond_style=s.bond_style,
                        angle_style=s.angle_style,
                        dihedral_style=s.dihedral_style,
                        improper_style=s.improper_style)
    for p in n.particles:
        s.particles[p.tag].x = p.x
        s.particles[p.tag].y = p.y
        s.particles[p.tag].z = p.z
    os.remove('temp.in')

    try:
        os.remove('temp.lmps')
    except OSError as e:
        print e

    if not write:
        try:
            os.remove('pysimm_relax.lmps')
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: nve/limit relaxation using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        except OSError:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: nve/limit relaxation using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))
            return False

    else:
        if os.path.isfile(write):
            if name:
                print('%s: %s simulation using LAMMPS successful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: nve/limit relaxation using LAMMPS successful'
                      % (strftime('%H:%M:%S')))
            return True
        else:
            if name:
                print('%s: %s simulation using LAMMPS UNsuccessful'
                      % (strftime('%H:%M:%S'), name))
            else:
                print('%s: nve/limit relaxation using LAMMPS UNsuccessful'
                      % strftime('%H:%M:%S'))


def write_init(l, **kwargs):
    """pysimm.lmps.write_init

    Create initialization LAMMPS input based on pysimm.system.System data

    Args:
        l: pysimm.system.System object reference
        kwargs:
            atom_style: LAMMPS atom_style default=full
            kspace_style: LAMMPS kspace style default='pppm 1e-4'
            units: LAMMPS set of units to use default=real
            special_bonds: LAMMPS special bonds input
            nonbond_mixing: type of mixing rule for nonbonded interactions default=arithmetic
            nb_cut: cutoff for nonbonded interactions default=12
    """
    atom_style = kwargs.get('atom_style') or 'full'
    kspace_style = kwargs.get('kspace_style') or 'pppm 1e-4'
    units = kwargs.get('units') or 'real'
    nb_cut = kwargs.get('nb_cut') or 12.0
    special_bonds = kwargs.get('special_bonds')
    nonbond_mixing = kwargs.get('nonbond_mixing')
    inner_cutoff = kwargs.get('inner_cutoff', 8.0)

    output = ''

    if type(l) == str and os.path.isfile(l):
        l = read_lammps(l, quiet=True)
    elif type(l) == str:
        return 'init_system failed to read %s' % l
    output += 'units %s\n' % units
    output += 'atom_style %s\n' % atom_style
    pair_style = None
    charge = False
    
    l.write_lammps('temp.lmps')

    if l.charge is None:
        for p in l.particles:
            if p.charge != 0:
                charge = True
                break
    else:
        if l.charge != 0:
            charge = True

    if not l.pair_style:
        if l.particle_types[1].sigma and l.particle_types[1].epsilon:
            if charge:
                if l.ff_class == '2':
                    pair_style = 'lj/class2/coul/long'
                else:
                    pair_style = 'lj/cut/coul/long'
            else:
                if l.ff_class == '2':
                    pair_style = 'lj/class2'
                else:
                    pair_style = 'lj/cut'

        elif (l.particle_types[1].a and l.particle_types[1].rho and
                l.particle_types[1].c):
            if charge:
                pair_style = 'buck/coul/long'
            else:
                pair_style = 'buck'
    else:
        if l.pair_style.startswith('lj') or l.pair_style.startswith('class2'):
            if charge:
                if l.ff_class == '2':
                    pair_style = 'lj/class2/coul/long'
                else:
                    pair_style = 'lj/cut/coul/long'
            else:
                if l.ff_class == '2':
                    pair_style = 'lj/class2'
                else:
                    pair_style = 'lj/cut'
        elif l.pair_style.startswith('buck'):
            if charge:
                pair_style = 'buck/coul/long'
            else:
                pair_style = 'buck'
        elif l.pair_style.startswith('charmm'):
            pair_style = 'lj/charmm/coul/long'

    if pair_style.startswith('lj/charmm'):
        output += 'pair_style %s %s %s\n' % (pair_style, inner_cutoff, nb_cut)
    elif pair_style:
        output += 'pair_style %s %s\n' % (pair_style, nb_cut)
    else:
        error_print('pair style probably not supported')

    if charge:
        output += 'kspace_style %s\n' % kspace_style

    if not pair_style.startswith('buck'):
        if nonbond_mixing == 'arithmetic' or l.mixing_rule.startswith('arith') or (l.forcefield and l.forcefield.lower() in ['amber', 'gaff', 'gaff2', 'dreiding', 'charmm', 'cgenff']):
            output += 'pair_modify mix arithmetic\n'
        elif nonbond_mixing == 'geometric' or l.mixing_rule.startswith('geo') or (l.forcefield and l.forcefield.lower() in ['opls']):
            output += 'pair_modify mix geometric\n'
        elif 'class2' in pair_style or l.mixing_rule.startswith('sixth') or (l.forcefield and l.forcefield.lower() in ['pcff']):
            output += 'pair_modify mix sixthpower\n'

    if l.bond_style:
        output += 'bond_style %s\n' % l.bond_style
    else:
        if l.ff_class == '1':
            output += 'bond_style harmonic\n'
        elif l.ff_class == '2':
            output += 'bond_style class2\n'

    if l.angles.count > 0:
        if l.angle_style:
            output += 'angle_style %s\n' % l.angle_style
        else:
            if l.ff_class == '1':
                output += 'angle_style harmonic\n'
            elif l.ff_class == '2':
                output += 'angle_style class2\n'

    if l.dihedrals.count > 0:
        if l.dihedral_style:
            output += 'dihedral_style %s\n' % l.dihedral_style
        else:
            if l.ff_class == '1':
                output += 'dihedral_style harmonic\n'
            elif l.ff_class == '2':
                output += 'dihedral_style class2\n'

    if l.impropers.count > 0:
        if l.improper_style:
            output += 'improper_style %s\n' % l.improper_style
        else:
            if l.ff_class == '1':
                output += 'improper_style harmonic\n'
            elif l.ff_class == '2':
                output += 'improper_style class2\n'

    if special_bonds:
        output += 'special_bonds %s\n' % special_bonds
    else:
        if 'class2' in pair_style or l.forcefield == 'pcff':
            output += 'special_bonds    lj 0 0 1 coul 0 0 1\n'
        elif 'charmm' in pair_style or l.forcefield == 'charmm':
            output += 'special_bonds charmm\n'
        elif l.forcefield == 'amber' or (l.forcefield and l.forcefield.startswith('gaff')):
            output += 'special_bonds amber\n'
        elif l.forcefield == 'dreiding':
            output += 'special_bonds dreiding\n'
        elif l.forcefield == 'opls':
            output += 'special_bonds lj 0 0 0.5 coul 0 0 0.5\n'

    output += 'read_data temp.lmps\n'

    if pair_style.startswith('buck'):
        for pt1 in l.particle_types:
            for pt2 in l.particle_types:
                if pt1.tag <= pt2.tag:
                    a = pow(pt1.a*pt2.a, 0.5)
                    c = pow(pt1.c*pt2.c, 0.5)
                    rho = 0.5*(pt1.rho+pt2.rho)
                    output += 'pair_coeff %s %s %s %s %s\n' % (pt1.tag, pt2.tag, a, rho, c)

    return output
