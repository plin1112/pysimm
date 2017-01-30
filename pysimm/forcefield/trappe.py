# ******************************************************************************
# pysimm.forcefield.gaff module
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

import os
from itertools import permutations, combinations

from pysimm import error_print
from pysimm.system import Angle, Dihedral, Improper
from pysimm.forcefield import Forcefield


class Trappe(Forcefield):
    """pysimm.forcefield.Trappe

    Forcefield object with typing rules for Gaff model.
    By default reads data file in forcefields subdirectory.

    Attributes:
        ff_name: trappe
        pair_style: lj
        ff_class: 1
    """
    def __init__(self, db_file=None):
        if not db_file and db_file is not False:
            db_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   os.pardir, os.pardir, 'dat', 'forcefields', 'trappe.xml')
        Forcefield.__init__(self, db_file)
        self.ff_name = 'trappe'
        self.pair_style = 'lj'
        self.ff_class = '1'

    def assign_ptypes(self, s):
        """pysimm.forcefield.Trappe.assign_ptypes

        Trappe specific particle typing rules.
        Requires System object Particle objects have Particle.bonds defined.
        *** use System.add_particle_bonding() to ensure this ***

        *** Not entirely inclusive - some atom types not used ***

        Args:
            s: pysimm.system.System

        Returns:
            None
        """
        all_types = set()
        s.pair_style = self.pair_style
        for p in s.particles:
            p.bonded_to = [x.a if p is x.b else x.b for x in p.bonds]
            p.bond_orders = [x.order for x in p.bonds]
            if None in p.bond_orders:
                error_print('error: bond orders are not set')
            p.bond_elements = [x.a.elem if p is x.b else x.b.elem for x in
                               p.bonds]
            p.nbonds = len(p.bond_elements)
            if p.linker:
                p.nbonds += 1
        for p in s.particles:
            if p.elem == 'C':
                if p.nbonds == 3 and 'O' in p.bond_elements:
                    p.type_name = 'c_trappe'
                elif (p.nbonds + p.implicit_h) == 4:
                    if p.implicit_h == 0:
                        p.type_name = 'ch0_trappe'
                    elif p.implicit_h == 1:
                        p.type_name = 'ch1_trappe'
                    elif p.implicit_h == 2:
                        p.type_name = 'ch2_trappe'
                    elif p.implicit_h == 3:
                        p.type_name = 'ch3_trappe'
                    elif p.implicit_h == 4:
                        p.type_name = 'ch4_trappe'
            elif p.elem == 'O':
                if p.nbonds == 1:
                    p.type_name = 'o_trappe'
                elif p.bond_elements.count('H') == 2:
                    p.type_name = 'ow'
                elif p.bond_elements.count('H') == 1:
                    p.type_name = 'oh'
                else:
                    p.type_name = 'os_trappe'
            else:
                print 'cant type particle %s' % p.tag
                return p

            type_ = self.particle_types.get(p.type_name)
            if not type_:
                print(p.tag, p.elem, p.type_name)

            all_types.add(self.particle_types.get(p.type_name)[0])

        for pt in all_types:
            pt = pt.copy()
            s.particle_types.add(pt)

        for p in s.particles:
            pt = s.particle_types.get(p.type_name)
            if pt:
                pt = pt[0]
                p.trappe_name = pt.name
                pt.name = p.type.name
                p.type = pt
                
    def assign_btypes(self, s):
        pass
    
    def assign_atypes(self, s):
        pass
                
    def assign_dtypes(self, s):
        pass
    
    def assign_itypes(self, s):
        pass
    
    def assign_charges(self, s, charges=None):
        pass