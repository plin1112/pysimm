"""
The pysimm.calc module contains functions for performing calculations on molecular systems.
This includes geometric calculations like distances, angles, and dihedrals, as well as
force field energy calculations.
"""
# ******************************************************************************
# pysimm.calc module
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

from random import random
from math import sin, cos, pi, acos
from typing import List, Tuple, Optional, TYPE_CHECKING, Union

if TYPE_CHECKING:
    from pysimm.system import Particle, System, ParticleType, BondType, AngleType, DihedralType, ImproperType

try:
    import numpy as np
    NDArray = np.ndarray
except ImportError:
    np = None
    NDArray = 'typing.Any'

from pysimm import error_print
from pysimm import warning_print
from pysimm import verbose_print
from pysimm import debug_print
from pysimm import PysimmError


def intersection(line1: List[List[float]], line2: List[List[float]]) -> Tuple[float, float]:
    """
    Finds the intersection between two 2D lines.

    Args:
        line1: A list of two points defining the first line, e.g., [[x1, y1], [x2, y2]].
        line2: A list of two points defining the second line.

    Returns:
        A tuple containing the x and y coordinates of the intersection point.

    Raises:
        PysimmError: If the lines do not intersect.
    """
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a: List[float], b: List[float]) -> float:
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        raise PysimmError('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def find_rotation(a: List[float], b: List[float]) -> NDArray[np.float64]:
    """
    Finds the rotation matrix required to align vector a with vector b.

    Args:
        a: A 3D vector [x, y, z].
        b: A 3D vector [x, y, z].

    Returns:
        The rotation matrix.
    """
    if not np:
        raise PysimmError('pysimm.calc.find_rotation function requires numpy')
    a_np = np.array(a)
    b_np = np.array(b)

    a_x_b = np.cross(a_np, b_np)
    axis = a_x_b / np.linalg.norm(a_x_b)
    theta = acos(np.dot(a_np, b_np) / (np.linalg.norm(a_np) * np.linalg.norm(b_np)))

    skew = np.matrix([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])

    rot_matrix = np.identity(3) + sin(theta) * skew + (1 - cos(theta)) * skew * skew
    return rot_matrix


def rotate_vector(x: float, y: float, z: float, theta_x: Optional[float] = None, theta_y: Optional[float] = None, theta_z: Optional[float] = None) -> List[float]:
    """
    Rotates a 3D vector around the x, y, and z axes by given angles.

    Args:
        x: The x component of the vector.
        y: The y component of the vector.
        z: The z component of the vector.
        theta_x: The angle to rotate around the x-axis.
        theta_y: The angle to rotate around the y-axis.
        theta_z: The angle to rotate around the z-axis.

    Returns:
        The new rotated vector [x, y, z].
    """
    if not np:
        raise PysimmError('pysimm.calc.rotate_vector function requires numpy')
    xt = random() * 2 * pi if theta_x is None else theta_x
    yt = random() * 2 * pi if theta_y is None else theta_y
    zt = random() * 2 * pi if theta_z is None else theta_z

    c = np.array((x, y, z))
    rot_mat_x = np.array([[1, 0, 0],
                          [0, cos(xt), -sin(xt)],
                          [0, sin(xt), cos(xt)]])
    rot_mat_y = np.array([[cos(yt), 0, sin(yt)],
                          [0, 1, 0],
                          [-sin(yt), 0, cos(yt)]])
    rot_mat_z = np.array([[cos(zt), -sin(zt), 0],
                          [sin(zt), cos(zt), 0],
                          [0, 0, 1]])

    return list(np.dot(rot_mat_z, np.dot(rot_mat_y, np.dot(rot_mat_x, c))))


def distance(p1: 'Particle', p2: 'Particle') -> float:
    """
    Calculates the distance between two particles.

    This function calculates the Euclidean distance between the coordinates of two
    :class:`~pysimm.system.Particle` objects, without considering periodic
    boundary conditions.

    Args:
        p1: The first particle.
        p2: The second particle.

    Returns:
        The distance between the two particles.
    """
    if not np:
        raise PysimmError('pysimm.calc.distance function requires numpy')
    return np.linalg.norm(p1.pos - p2.pos)


def angle(p1: 'Particle', p2: 'Particle', p3: 'Particle', radians: bool = False) -> float:
    """
    Calculates the angle formed by three particles.

    Args:
        p1: The first particle.
        p2: The central particle.
        p3: The third particle.
        radians: If True, returns the angle in radians. Otherwise, in degrees.

    Returns:
        The angle between the particles.
    """
    if not np:
        raise PysimmError('pysimm.calc.angle function requires numpy')
    v1 = p1.pos - p2.pos
    v2 = p3.pos - p2.pos
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    if v1_norm == 0 or v2_norm == 0:
        return 0.0
    cos_theta = np.dot(v1, v2) / (v1_norm * v2_norm)
    theta = acos(np.clip(cos_theta, -1.0, 1.0))
    if not radians:
        theta = theta * 180 / pi
    return theta


def dihedral(p1: 'Particle', p2: 'Particle', p3: 'Particle', p4: 'Particle', radians: bool = False) -> float:
    """
    Calculates the dihedral angle formed by four particles.

    Args:
        p1: The first particle.
        p2: The second particle.
        p3: The third particle.
        p4: The fourth particle.
        radians: If True, returns the angle in radians. Otherwise, in degrees.

    Returns:
        The dihedral angle between the particles.
    """
    if not np:
        raise PysimmError('pysimm.calc.dihedral function requires numpy')
    b1 = p2.pos - p1.pos
    b2 = p3.pos - p2.pos
    b3 = p4.pos - p3.pos

    n1 = np.cross(b1, b2)
    n1_norm = np.linalg.norm(n1)
    if n1_norm > 1e-6:
        n1 /= n1_norm

    n2 = np.cross(b2, b3)
    n2_norm = np.linalg.norm(n2)
    if n2_norm > 1e-6:
        n2 /= n2_norm

    b2_norm = np.linalg.norm(b2)
    if b2_norm < 1e-6:
        return 0.0 # Cannot compute dihedral angle

    m1 = np.cross(n1, b2 / b2_norm)

    x = np.clip(np.dot(n1, n2), -1.0, 1.0)
    y = np.dot(m1, n2)

    theta = np.arctan2(y, x)

    if not radians:
        theta = theta * 180 / pi
    return theta


def chiral_angle(a: 'Particle', b: 'Particle', c: 'Particle', d: 'Particle') -> float:
    """
    Calculates the chiral angle for a group of four particles.

    The chiral angle is defined as the angle between the vector resulting from
    vec(a->c) X vec(a->d) and vec(a->b). This is used to help define tacticity
    where the backbone follows b'--a--b and c and d are side groups.

    Args:
        a: The central particle.
        b: A backbone particle bonded to a.
        c: A side group particle bonded to a.
        d: Another side group particle bonded to a.

    Returns:
        The chiral angle in degrees.
    """
    if not np:
        raise PysimmError('pysimm.calc.chiral_angle function requires numpy')

    ht = a.pos - b.pos
    ht_norm = np.linalg.norm(ht)
    if ht_norm > 0: ht /= ht_norm

    hmethyl = a.pos - d.pos
    hmethyl_norm = np.linalg.norm(hmethyl)
    if hmethyl_norm > 0: hmethyl /= hmethyl_norm

    hside = a.pos - c.pos
    hside_norm = np.linalg.norm(hside)
    if hside_norm > 0: hside /= hside_norm

    side_x_methyl = np.cross(hside, hmethyl)

    side_dot_methyl = np.dot(hside, hmethyl)
    side_theta_methyl = acos(np.clip(side_dot_methyl, -1.0, 1.0))

    if sin(side_theta_methyl) > 1e-6:
        side_x_methyl /= sin(side_theta_methyl)

    cos_theta = np.dot(side_x_methyl, ht)

    return acos(np.clip(cos_theta, -1.0, 1.0))/pi*180


def tacticity(s: 'System', a_tag: int, b_tag: int, c_tag: int, d_tag: int, offset: int, return_angles: bool = True, unwrap: bool = True,
              rewrap: bool = True, skip_first: bool = False) -> Union[Tuple[str, List[float]], str]:
    """
    Determines the tacticity for a polymer chain.

    This function iterates through groups of four particles given by their tags
    and an offset, assuming equivalent atoms in each monomer are perfectly offset.

    Args:
        s: The molecular system.
        a_tag: The tag of the first 'a' particle (chiral center).
        b_tag: The tag of the first 'b' particle.
        c_tag: The tag of the first 'c' particle.
        d_tag: The tag of the first 'd' particle.
        offset: The offset of particle tags between monomers.
        return_angles: If True, returns the chiral angles of all monomers.
        unwrap: If True, unwraps the system before calculation.
        rewrap: If True, rewraps the system after calculation.
        skip_first: If True, skips the first monomer.

    Returns:
        The tacticity ('isotactic', 'syndiotactic', 'atactic', or 'undetermined'),
        and optionally a list of chiral angles.
    """
    if not np:
        raise PysimmError('pysimm.calc.tacticity function requires numpy')
    if a_tag is None or b_tag is None or c_tag is None or d_tag is None:
        raise PysimmError('particle tags for chiral center are required (a, b, c, and d)')

    if offset is None:
        raise PysimmError('offset for tags in each monomer is required, i.e. - number of particles in each monomer')

    if unwrap:
        s.unwrap()

    a_ = [s.particles[i] for i in range(a_tag, s.particles.count, offset)]
    b_ = [s.particles[i] for i in range(b_tag, s.particles.count, offset)]
    c_ = [s.particles[i] for i in range(c_tag, s.particles.count, offset)]
    d_ = [s.particles[i] for i in range(d_tag, s.particles.count, offset)]

    stereochem_angles = []

    if skip_first:
        for a, b, c, d in zip(a_[1:], b_[1:], c_[1:], d_[1:]):
            stereochem_angles.append(chiral_angle(a, b, c, d))
    else:
        for a, b, c, d in zip(a_, b_, c_, d_):
            stereochem_angles.append(chiral_angle(a, b, c, d))

    last = None
    iso_diads = 0
    syn_diads = 0
    for a in stereochem_angles:
        if last is not None:
            if (a < 90 and last < 90) or (a > 90 and last > 90):
                iso_diads += 1
            else:
                syn_diads += 1
        last = a

    if not stereochem_angles:
        t = 'undetermined'
    elif iso_diads == (len(stereochem_angles) - 1):
        t = 'isotactic'
    elif syn_diads == (len(stereochem_angles) - 1):
        t = 'syndiotactic'
    else:
        t = 'atactic'

    if rewrap:
        s.wrap()

    if syn_diads + iso_diads > 0:
        print('{:.1f}% syn diads'.format(100*float(syn_diads)/(syn_diads+iso_diads)))
        print('{:.1f}% iso diads'.format(100*float(iso_diads)/(syn_diads+iso_diads)))
    if return_angles:
        return t, stereochem_angles
    else:
        return t


def frac_free_volume(v_sp: float, v_void: float) -> float:
    """
    Determines the fractional free volume for a porous system.

    Args:
        v_sp: The specific volume.
        v_void: The void volume.

    Returns:
        The fractional free volume.
    """
    return (-0.3 * v_sp + 1.3 * v_void)/v_sp


def pbc_distance(s: 'System', p1: 'Particle', p2: 'Particle') -> float:
    """
    Calculates the distance between two particles considering periodic boundary conditions.

    Args:
        s: The molecular system containing the particles.
        p1: The first particle.
        p2: The second particle.

    Returns:
        The distance between the particles, considering PBC.
    """
    if not np:
        raise PysimmError('pysimm.calc.pbc_distance function requires numpy')
    box_dims = np.array([s.dim.dx, s.dim.dy, s.dim.dz])
    frac_pos1 = p1.pos / box_dims
    frac_pos2 = p2.pos / box_dims

    frac_d = frac_pos1 - frac_pos2
    frac_d = frac_d - np.round(frac_d)

    d_vec = frac_d * box_dims

    return np.linalg.norm(d_vec)


def LJ_12_6(pt: 'ParticleType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the 12-6 Lennard-Jones potential."""
    inv_d = pt.sigma / d
    return 4 * pt.epsilon * (np.power(inv_d, 12) - np.power(inv_d, 6))


def LJ_9_6(pt: 'ParticleType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the 9-6 Lennard-Jones potential."""
    inv_d = pt.sigma / d
    return pt.epsilon * (2 * np.power(inv_d, 12) - 3 * np.power(inv_d, 6))


def buckingham(pt: 'ParticleType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the Buckingham potential."""
    return pt.a * np.exp(-d / pt.rho) - (pt.c / np.power(d, 6))


def harmonic_bond(bt: 'BondType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the harmonic bond potential."""
    return bt.k * np.power(d - bt.r0, 2)


def class2_bond(bt: 'BondType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the class 2 bond potential."""
    return bt.k2 * np.power(d - bt.r0, 2) + bt.k3 * np.power(d - bt.r0, 3) + bt.k4 * np.power(d - bt.r0, 4)


def harmonic_angle(at: 'AngleType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the harmonic angle potential."""
    return at.k * np.power(d - at.theta0, 2)


def class2_angle(at: 'AngleType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the class 2 angle potential."""
    return at.k2 * np.power(d - at.theta0, 2) + at.k3 * np.power(d - at.theta0, 3) + at.k4 * np.power(d - at.theta0, 4)


def harmonic_dihedral(dt: 'DihedralType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the harmonic dihedral potential."""
    return dt.k * (1 + dt.d * np.cos(dt.n * np.radians(d)))


def class2_dihedral(dt: 'DihedralType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the class 2 dihedral potential."""
    return (
            dt.k1 * (1 - np.cos(np.radians(d) - np.radians(dt.phi1))) +
            dt.k2 * (1 - np.cos(2*np.radians(d) - np.radians(dt.phi2))) +
            dt.k3 * (1 - np.cos(3*np.radians(d) - np.radians(dt.phi3)))
    )


def opls_dihedral(dt: 'DihedralType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the OPLS dihedral potential."""
    return (
            0.5 * dt.k1 * (1 + np.cos(np.radians(d))) +
            0.5 * dt.k2 * (1 - np.cos(2 * np.radians(d))) +
            0.5 * dt.k3 * (1 + np.cos(3 * np.radians(d))) +
            0.5 * dt.k4 * (1 - np.cos(4 * np.radians(d)))
    )


def fourier_dihedral(dt: 'DihedralType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the Fourier dihedral potential."""
    d_rad = np.radians(d)
    energy = 0
    for k, n, d_ in zip(dt.k, dt.n, dt.d):
        energy += k * (1 + np.cos(n * d_rad - np.radians(d_)))
    return energy


def harmonic_improper(it: 'ImproperType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the harmonic improper potential."""
    return it.k * np.power(d - it.x0, 2)


def cvff_improper(it: 'ImproperType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the CVFF improper potential."""
    return it.k * (1 + it.d * np.cos(it.n * np.radians(d)))


def umbrella_improper(it: 'ImproperType', d: Union[float, NDArray[np.float64]]) -> Union[float, NDArray[np.float64]]:
    """Calculates the umbrella improper potential."""
    return it.k * (1 - np.cos(np.radians(d)))
