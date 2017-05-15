# ******************************************************************************
# pysimm.forcefield module
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
import json
from xml.dom import minidom
from xml.etree import ElementTree as Et
from itertools import permutations, combinations

import gasteiger
from pysimm import error_print, warning_print, debug_print
from pysimm.utils import PysimmError, Item, ItemContainer, compare
from pysimm.system import (ParticleType, BondType, AngleType,
    DihedralType, ImproperType, ParticleTypeContainer, BondTypeContainer, 
    AngleTypeContainer, DihedralTypeContainer, ImproperTypeContainer,
    Angle, Dihedral, Improper)


element_names_by_mass = {1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
                         14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na',
                         24: 'Mg', 27: 'Al', 28: 'Si', 31: 'P', 32: 'S',
                         35: 'Cl', 39: 'K', 40: 'Ca', 80: 'Br', 127: 'I'}


class Forcefield(object):
    """pysimm.forcefield.Forcefield

    Base Forcefield class definition.
    Initialize with force field xml file.

    Attributes:
        ff_class: force field class (1 or 2)
        ff_name: force field name
        particle_types: pysimm.utils.ItemContainer for particle_types
        bond_types: pysimm.utils.ItemContainer for bond_types
        angle_types: pysimm.utils.ItemContainer for angle_types
        dihedral_types: pysimm.utils.ItemContainer for dihedral_types
        improper_types: pysimm.utils.ItemContainer for improper_types
    """
    def __init__(self, file_=None, format=None):
        self.ff_class = 0
        self.ff_name = ''
        self.pair_style = None
        self.particle_types = ParticleTypeContainer()
        self.bond_style = None
        self.bond_types = BondTypeContainer()
        self.angle_style = None
        self.angle_types = AngleTypeContainer()
        self.dihedral_style = None
        self.dihedral_types = DihedralTypeContainer()
        self.improper_style = None
        self.improper_types = ImproperTypeContainer()

        if file_:
            if format == 'json' or file_.split('.')[-1] == 'json':
                self.from_json(file_)
            elif format == 'xml' or file_.split('.')[-1] == 'xml':
                self.from_xml(file_)
            
    
    def from_xml(self, file_):
        tree = Et.parse(file_)
        root = tree.getroot()
        ptypes = root.find('ParticleTypes')
        btypes = root.find('BondTypes')
        atypes = root.find('AngleTypes')
        dtypes = root.find('DihedralTypes')
        itypes = root.find('ImproperTypes')

        for ptype in ptypes:
            pt = ParticleType()
            for k, v in ptype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(pt, k, v)
            pt.tag = pt.name
            self.particle_types.add(pt)

        for btype in btypes:
            bt = BondType()
            for k, v in btype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(bt, k, v)
            bt.tag = bt.name
            bt.rname = ','.join(reversed(bt.name.split(',')))
            self.bond_types.add(bt)

        for atype in atypes:
            at = AngleType()
            for k, v in atype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(at, k, v)
            at.tag = at.name
            at.rname = ','.join(reversed(at.name.split(',')))
            self.angle_types.add(at)

        for dtype in dtypes:
            dt = DihedralType()
            for k, v in dtype.attrib.items():
                if k != 'd' and k != 'n':
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                elif k == 'd' or k == 'n':
                    try:
                        v = int(v)
                    except ValueError:
                        pass
                setattr(dt, k, v)
            dt.tag = dt.name
            dt.rname = ','.join(reversed(dt.name.split(',')))
            self.dihedral_types.add(dt)

        for itype in itypes:
            it = ImproperType()
            for k, v in itype.attrib.items():
                try:
                    v = float(v)
                except ValueError:
                    pass
                setattr(it, k, v)
            it.tag = it.name
            it.rname = ','.join(reversed(it.name.split(',')))
            self.improper_types.add(it)
    
    
    def from_json(self, json_file):
        with file(json_file) as f:
            j = json.loads(f.read())
        self.ff_name = j.get('ff_name')
        self.ff_class = j.get('ff_class')
        self.pair_style = j.get('pair_style')
        self.bond_style = j.get('bond_style')
        self.angle_style = j.get('angle_style')
        self.dihedral_style = j.get('dihedral_style')
        self.improper_style = j.get('improper_style')
        
        for pt in j.get('particle_types'):
            self.particle_types.add(ParticleType(**pt))
        
        for bt in j.get('bond_types'):
            self.bond_types.add(BondType(**bt))
        
        for at in j.get('angle_types'):
            self.angle_types.add(AngleType(**at))
        
        for dt in j.get('dihedral_types'):
            self.dihedral_types.add(DihedralType(**dt))
        
        for it in j.get('improper_types'):
            self.improper_types.add(ImproperType(**it))
        
            
    def write_json(self, out):
        f = {}
        f['ff_name'] = self.ff_name
        f['ff_class'] = self.ff_class
        f['pair_style'] = self.pair_style
        f['bond_style'] = self.bond_style
        f['angle_style'] = self.angle_style
        f['dihedral_style'] = self.dihedral_style
        f['improper_style'] = self.improper_style
        f['particle_types'] = [vars(pt) for pt in self.particle_types]
        f['bond_types'] = [vars(bt) for bt in self.bond_types]
        f['angle_types'] = [vars(at) for at in self.angle_types]
        f['dihedral_types'] = [vars(dt) for dt in self.dihedral_types]
        f['improper_types'] = [vars(it) for it in self.improper_types]
        with file(out, 'w') as _file:
            _file.write(json.dumps(f, indent=4))
            

    def write_xml(self, out):
        """pysimm.forcefield.Forcefield.write

        Write Forcefield object to xml format.

        Args:
            out: file name to write

        Returns:
            None
        """

        ff_elem = Et.Element('Forcefield')
        ff_elem.set('name', self.ff_name)
        ff_elem.set('class', self.ff_class)

        ptypes = Et.SubElement(ff_elem, 'ParticleTypes')
        for pt in self.particle_types:
            ptype = Et.SubElement(ptypes, 'Type')
            for k, v in vars(pt).items():
                ptype.set(k, str(v))

        btypes = Et.SubElement(ff_elem, 'BondTypes')
        for bt in self.bond_types:
            btype = Et.SubElement(btypes, 'Type')
            for k, v in vars(bt).items():
                btype.set(k, str(v))

        atypes = Et.SubElement(ff_elem, 'AngleTypes')
        for at in self.angle_types:
            atype = Et.SubElement(atypes, 'Type')
            for k, v in vars(at).items():
                atype.set(k, str(v))

        dtypes = Et.SubElement(ff_elem, 'DihedralTypes')
        for dt in self.dihedral_types:
            dtype = Et.SubElement(dtypes, 'Type')
            for k, v in vars(dt).items():
                dtype.set(k, str(v))

        itypes = Et.SubElement(ff_elem, 'ImproperTypes')
        for it in self.improper_types:
            itype = Et.SubElement(itypes, 'Type')
            for k, v in vars(it).items():
                itype.set(k, str(v))

        with open(out, 'w') as f:
            f.write(
                minidom.parseString(Et.tostring(ff_elem, 'utf-8')).toprettyxml(
                    indent="  "))
                    
    def assign_bond_types(self, s, order=False):
        untypable = ItemContainer()
        s.bond_style = self.bond_style
        for b in s.bonds:
            bt = self.bond_types.by_name(b.ptype_name(), order=order)
            if bt:
                s_bt = s.bond_types.by_name(bt[0].name, exact=True, order=order)
                if s_bt:
                    b.type = s_bt
                else:
                    b.type = s.bond_types.add(bt[0].copy())
            else:
                untypable.add(b)
        return untypable
                
    def assign_angle_types(self, s, order=False):
        untypable = ItemContainer()
        s.angle_style = self.angle_style
        for p in s.particles:
            try:
                p.bonded_to.count
            except:
                raise PysimmError('bonded_to not defined for particle {}: run System.add_particle_bonding() first'.format(p.tag))
            for p1 in p.bonded_to:
                for p2 in p.bonded_to:
                    if p2.tag > p1.tag:
                        at = self.angle_types.by_name('{},{},{}'.format(
                            p1.type.eq_angle or p1.type.name,
                            p.type.eq_angle or p.type.name,
                            p2.type.eq_angle or p2.type.name
                        ), order=order)
                        if at:
                            a = s.angles.add(Angle(a=p1, b=p, c=p2))
                            s_at = s.angle_types.by_name(at[0].name, exact=True, order=order)
                            if s_at:
                                a.type = s_at
                            else:
                                a.type = s.angle_types.add(at[0].copy())
                        else:
                            untypable.add(Angle(a=p1, b=p, c=p2))
        return untypable
        
    def assign_dihedral_types(self, s, order=False):
        untypable = ItemContainer()
        s.dihedral_style = self.dihedral_style
        for b in s.bonds:
            try:
                b.a.bonded_to.count
            except:
                raise PysimmError('bonded_to not defined for particle {}: run System.add_particle_bonding() first'.format(b.a.tag))
            
            try:
                b.b.bonded_to.count
            except:
                raise PysimmError('bonded_to not defined for particle {}: run System.add_particle_bonding() first'.format(b.b.tag))
            for p1 in b.a.bonded_to:
                for p2 in b.b.bonded_to:
                    if (p1 is not b.b and p2 is not b.a):
                        p1_name = p1.type.eq_dihedral or p1.type.name
                        a_name = b.a.type.eq_dihedral or b.a.type.name
                        b_name = b.b.type.eq_dihedral or b.b.type.name
                        p2_name = p2.type.eq_dihedral or p2.type.name
                        dt = self.dihedral_types.by_name(
                            '{},{},{},{}'.format(
                                p1_name, a_name, b_name, p2_name
                            ), order=order
                        )
                        if dt:
                            d = s.dihedrals.add(
                                Dihedral(
                                    a=p1, b=b.a, c=b.b, d=p2
                                )
                            )
                            s_dt = s.dihedral_types.by_name(dt[0].name, exact=True, order=order)
                            if s_dt:
                                d.type = s_dt
                            else:
                                d.type = s.dihedral_types.add(dt[0].copy())
                        else:
                            untypable.add(Dihedral(a=p1, b=b.a, c=b.b, d=p2))
        return untypable
        
    def assign_improper_types(self, s, order=True):
        s.improper_style = self.improper_style
        for p in s.particles:
            try:
                p.bonded_to.count
            except:
                raise PysimmError('bonded_to not defined for particle {}: run System.add_particle_bonding() first'.format(p.tag))
            if len(p.bonded_to) == 3:
                for perm in permutations(p.bonded_to, 3):
                    p_name = p.type.eq_improper or p.type.name
                    p1_name = perm[0].type.eq_improper or perm[0].type.name
                    p2_name = perm[1].type.eq_improper or perm[1].type.name
                    p3_name = perm[2].type.eq_improper or perm[2].type.name
                    it = self.improper_types.by_name(
                        ','.join([
                            p_name, p1_name, p2_name, p3_name
                        ]), order=order
                    )
                    if it:
                        i = s.impropers.add(
                            Improper(
                                a=p, b=perm[0], c=perm[1], d=perm[2]
                            )
                        )
                        s_it = s.improper_types.by_name(it[0].name, exact=True, order=order)
                        if s_it:
                            i.type = s_it
                        else:
                            i.type = s.improper_types.add(it[0].copy())
                        break