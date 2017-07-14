# polymatic module; calls Polymatic perl code and LAMMPS

import os
import sys
import shlex
import shutil
from time import strftime
from subprocess import Popen, PIPE

from pysimm import system, lmps

rappture = True
try:
    import Rappture
except ImportError:
    rappture = False


def pack(script, file_in, nrep, boxl, file_out):
    """pysimm.apps.polymatic.pack

    Calls Polymatic random packing code

    Args:
        script: name of packing script
        file_in: list of file names of reference molecules to pack
        nrep: list of number of monomers for each reference molecule
        boxl: length of one dimension of simulation box for random packing
        file_out: name of output file (packed system)

    Returns:
        output from perl code
    """
    if not isinstance(file_in, list):
        file_in = [file_in]
    if not isinstance(nrep, list):
        nrep = [nrep]
    if len(file_in) != len(nrep) or len(file_in) == 0:
        return False

    cmd = 'perl %s -i ' % script
    cmd += '%s ' % str(len(file_in))
    for i in range(len(file_in)):
        cmd += '%s %s ' % (file_in[i], nrep[i])
    cmd += '-l %s -o %s' % (boxl, file_out)

    o, e = Popen(shlex.split(cmd),
                 stdin=PIPE,
                 stdout=PIPE,
                 stderr=PIPE).communicate()

    if not e and o:
        return o
    else:
        return False


def polymatic(script, file_in, file_out):
    """pysimm.apps.polymatic.polymatic

    Calls Polymatic code. polym.in and types.txt are assumed to exist.

    Args:
        script: name of Polymatic script
        file_in: initial system file name
        file_out: final system file name

    Returns:
        output from perl code
    """
    cmd = ('perl %s -i %s -s polym.in -t types.txt -o %s'
           % (script, file_in, file_out))
    o, e = Popen(shlex.split(cmd),
                 stdin=PIPE,
                 stdout=PIPE,
                 stderr=PIPE).communicate()

    if not e and o and o.split()[0] is not 'Error:':
        return True
    elif not e and o:
        return o
    else:
        return False


def run(settings):
    """pysimm.apps.polymatic.run

    Runs Polymatic algorithm.

    Args:
        settings: object containing Polymatic settings

    Returns:
        (True/False, pysimm.system.System)
    """
    if rappture:
        Rappture.Utils.progress(0, 'Initializing Polymatic...')
    bonds = 0

    os.mkdir('logs')

    polymatic(os.path.join(settings.polymatic_dir, 'polym_init.pl'),
              'data.lmps',
              'step_000.lmps')

    s = system.read_lammps('step_000.lmps', quiet=True)
    s.read_type_names('types.txt')
    s.write_lammps('temp.lmps')

    if rappture:
        Rappture.Utils.progress(0, '%s/%s bonds made: '
                                   'optimizing initial structure...'
                                % (bonds, settings.polym.max_bonds))

    if not lmps_min(s, 'initial optimization', settings):
        s.write_lammps('temp.lmps')
        polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                  'temp.lmps',
                  'final.lmps')
        return False, s

    s.write_lammps('step_000.lmps')
    s.write_lammps('temp.lmps')

    while bonds < settings.polym.max_bonds:
        attempt = 0
        while not polymatic(os.path.join(settings.polymatic_dir, 'polym.pl'),
                            'temp.lmps',
                            'temp.lmps'):
            attempt += 1
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds *
                                            100),
                                        '%s/%s bonds made: attempt #%s to make '
                                        'new bond'
                                        % (bonds, settings.polym.max_bonds,
                                           attempt))
            s = system.read_lammps('temp.lmps', quiet=True)
            s.read_type_names('types.txt')

            if not lmps_step_md(s, bonds, attempt, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

            if attempt >= settings.polym.max_md:
                break

        if attempt >= settings.polym.max_md:
            break

        bonds += 1

        if rappture:
            Rappture.Utils.progress(int(float(bonds)/settings.polym.max_bonds
                                        * 100),
                                    '%s/%s bonds made: '
                                    'optimizing newly formed bond'
                                    % (bonds, settings.polym.max_bonds))

        s = system.read_lammps('temp.lmps', quiet=True)
        s.read_type_names('types.txt')

        print('%s: bond %s made successfully' % (strftime('%H:%M:%S'), bonds))
        sys.stdout.flush()
        if not lmps_min(s, 'bond %s optimization' % bonds, settings):
            s.write_lammps('temp.lmps')
            polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                      'temp.lmps',
                      'final.lmps')
            return False, s
        s.write_lammps('step_%03d.lmps' % bonds)
        s.write_lammps('temp.lmps')

        if (bonds % settings.polym.cycle == 0 and
                (bonds / settings.polym.cycle) % settings.polym.npt_freq == 0):
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing npt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            if not lmps_cycle_npt_md(s, bonds, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

        elif bonds % settings.polym.cycle == 0:
            if rappture:
                Rappture.Utils.progress(int(float(bonds)/settings.
                                            polym.max_bonds
                                            * 100),
                                        '%s/%s bonds made: '
                                        'performing nvt cycle md'
                                        % (bonds, settings.polym.max_bonds))
            if not lmps_cycle_nvt_md(s, bonds, settings):
                s.write_lammps('temp.lmps')
                polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
                          'temp.lmps',
                          'final.lmps')
                return False, s
            s.write_lammps('temp.lmps')

    if rappture:
        Rappture.Utils.progress(99, 'Finalizing Polymatic')

    polymatic(os.path.join(settings.polymatic_dir, 'polym_final.pl'),
              'temp.lmps',
              'final.lmps')

    return True, s


def lmps_min(s, name, settings):
    """pysimm.apps.polymatic.lmps_min

    Runs LAMMPS minimization for the Polymatic algorithm.

    Args:
        s: pysimm.system.System to minimize
        name: name of simulation
        settings: object containing Polymatic settings

    Returns:
        result from lmps.minimize
    """
    if settings.polym.min.cluster:
        nanohub = {'cores': int(settings.polym.min.nanohub_cores),
                   'walltime': int(settings.polym.min.nanohub_walltime)}
        log_name = '%s' % '_'.join(name.split())
    else:
        nanohub = {}
        log_name = 'logs/%s' % '_'.join(name.split())

    if settings.polym.min.user_input:
        lmps.run(s, lammps_in=settings.polym.min.min_in,
                 name='initial optimization',
                 print_to_screen=False, nanohub=nanohub)
    else:
        result = lmps.minimize(s, name=name,
                               cutoff=settings.polym.min.nb_cutoff,
                               sd_etol=settings.polym.min.sd_etol,
                               sd_ftol=settings.polym.min.sd_ftol,
                               sd_maxiter=settings.polym.min.sd_maxiter,
                               sd_maxeval=settings.polym.min.sd_maxeval,
                               cg_etol=settings.polym.min.cg_etol,
                               cg_ftol=settings.polym.min.cg_ftol,
                               cg_maxiter=settings.polym.min.cg_maxiter,
                               cg_maxeval=settings.polym.min.cg_maxeval,
                               log=log_name,
                               np=settings.np,
                               nanohub=nanohub)

    if settings.polym.min.cluster:
        shutil.move(log_name, 'logs')

    return result


def lmps_step_md(s, bonds, attempt, settings):
    """pysimm.apps.polymatic.lmps_step_md

    Runs LAMMPS step md for the Polymatic algorithm.

    Args:
        s: pysimm.system.System to minimize
        bonds: number of bond to be made
        attempt: number of bonding attempt
        settings: object containing Polymatic settings

    Returns:
        result from lmps.md
    """
    if settings.polym.step.cluster:
        nanohub = {'cores': int(settings.polym.step.nanohub_cores),
                   'walltime': int(settings.polym.step.nanohub_walltime)}
        log_name = 'step_%03d_%03d' % (bonds, attempt)
    else:
        nanohub = {}
        log_name = 'logs/step_%03d_%03d' % (bonds, attempt)

    if settings.polym.step.user_input:
        lmps.run(s, lammps_in=settings.polym.step.step_in,
                 name='bond %s attempt #%d' % (bonds + 1, attempt),
                 print_to_screen=False, nanohub=nanohub)
    else:
        result = lmps.md(s, name='bond %s: attempt #%d' % (bonds + 1, attempt),
                         ensemble='nvt',
                         cutoff=settings.polym.step.nb_cutoff,
                         temp=settings.polym.step.temp,
                         new_v=True,
                         length=settings.polym.step.length,
                         log=log_name,
                         np=settings.np,
                         nanohub=nanohub)

    if settings.polym.step.cluster:
        shutil.move(log_name, 'logs')

    return result


def lmps_cycle_nvt_md(s, bonds, settings):
    """pysimm.apps.polymatic.lmps_cycle_nvt_md

    Runs LAMMPS nvt cycle md for the Polymatic algorithm.

    Args:
        s: pysimm.system.System to minimize
        bonds: number of bond to be made
        settings: object containing Polymatic settings

    Returns:
        result from lmps.md
    """
    if settings.polym.cycle_nvt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_nvt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_nvt.nanohub_walltime)}
        log_name = 'cycle_nvt_%03d' % bonds
    else:
        nanohub = {}
        log_name = 'logs/cycle_nvt_%03d' % bonds

    if settings.polym.cycle_nvt.user_input:
        lmps.run(s, lammps_in=settings.polym.cycle_nvt.step_in,
                 name='bond %d cycle nvt' % bonds,
                 print_to_screen=False, nanohub=nanohub)
    else:
        result = lmps.md(s, name='bond %d cycle nvt' % bonds,
                         ensemble='nvt',
                         cutoff=settings.polym.cycle_nvt.nb_cutoff,
                         temp=settings.polym.cycle_nvt.temp,
                         new_v=True,
                         length=settings.polym.cycle_nvt.length,
                         log=log_name,
                         np=settings.np,
                         nanohub=nanohub)

    if settings.polym.cycle_nvt.cluster:
        shutil.move(log_name, 'logs')

    return result


def lmps_cycle_npt_md(s, bonds, settings):
    """pysimm.apps.polymatic.lmps_cycle_npt_md

    Runs LAMMPS npt cycle md for the Polymatic algorithm.

    Args:
        s: pysimm.system.System to minimize
        bonds: number of bond to be made
        settings: object containing Polymatic settings

    Returns:
        result from lmps.md
    """
    if settings.polym.cycle_npt.cluster:
        nanohub = {'cores': int(settings.polym.cycle_npt.nanohub_cores),
                   'walltime': int(settings.polym.cycle_npt.nanohub_walltime)}
        log_name = 'cycle_npt_%03d' % bonds
    else:
        nanohub = {}
        log_name = 'logs/cycle_npt_%03d' % bonds

    if settings.polym.cycle_npt.user_input:
        lmps.run(s, lammps_in=settings.polym.cycle_npt.step_in,
                 name='bond %d cycle npt' % bonds,
                 print_to_screen=False, nanohub=nanohub)
    else:
        result = lmps.md(s, name='bond %d cycle npt' % bonds,
                         ensemble='nvt',
                         cutoff=settings.polym.cycle_npt.nb_cutoff,
                         temp=settings.polym.cycle_npt.temp,
                         new_v=True,
                         pressure=settings.polym.cycle_npt.pressure,
                         length=settings.polym.cycle_npt.length,
                         log=log_name,
                         np=settings.np,
                         nanohub=nanohub)

    if settings.polym.cycle_npt.cluster:
        shutil.move(log_name, 'logs')

    return result


def make_linker_types(s):
    """pysimm.apps.polymatic.make_linker_types

    Identifies linker particles and creates duplicate ParticleType objects with new names.
    Identification is performed by Particle.linker attribute.
    New ParticleType name is prepended with [H or T]L@ to designate head or tail linker

    Args:
        s: system to modify

    Returns:
        None
    """
    for p in s.particles:
        if p.linker == 'head':
            head_linker = s.particle_types.get('HL@%s' % p.type.name)
            if head_linker:
                p.type = head_linker[0]
            else:
                p.type = p.type.copy()
                p.type.name = 'HL@%s' % p.type.name
                s.particle_types.add(p.type)
        elif p.linker == 'tail':
            tail_linker = s.particle_types.get('TL@%s' % p.type.name)
            if tail_linker:
                p.type = tail_linker[0]
            else:
                p.type = p.type.copy()
                p.type.name = 'TL@%s' % p.type.name
                s.particle_types.add(p.type)
        elif p.linker:
            linker = s.particle_types.get('L@%s' % p.type.name)
            if linker:
                p.type = linker[0]
            else:
                p.type = p.type.copy()
                p.type.name = 'L@%s' % p.type.name
                s.particle_types.add(p.type)

def remove_linker_types(s):
    """pysimm.apps.polymatic.remove_linker_types

    Reassigns Particle.type references to original ParticleType objects without linker prepend

    Args:
        s: system to modify

    Returns:
        None
    """
    for p in s.particles:
        if p.type.name.find('@') >= 0:
            pt = s.particle_types.get(p.type.name.split('@')[-1])
            if pt:
                p.type = pt[0]
            else:
                warning_print('cannot find regular type for linker %s'
                      % p.type.name)