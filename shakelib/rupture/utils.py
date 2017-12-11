#!/usr/bin/env python

# stdlib modules
import copy

# third party imports
import numpy as np
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection

from impactutils.vectorutils.ecef import latlon2ecef
from impactutils.vectorutils.vector import Vector
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture import constants


def is_quad(q):
    """
    Checks that an individual quad is coplanar.

    Args:
        q (list): A quadrilateral; list of four OQ Points.

    Returns:
        tuple: Bool for whether or not the points are planar within tolerance;
            and also the corrected quad where p2 is adjusted to be on the same
            plane as the other points.
    """
    P0, P1, P2, P3 = q

    # Convert points to ECEF
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Unit vector along top edge
    v0 = (p1 - p0).norm()

    # Distance along bottom edge
    d = (p3 - p2).mag()

    # New location for p2 by extending from p3 the same distance and
    # direction that p1 is from p0:
    new_p2 = p3 + v0 * d

    # How far off of the plane is the origin p2?
    planepoints = [p0, p1, p2]
    dist = get_distance_to_plane(planepoints, p2)

    # Is it close enough?
    if dist / d > constants.OFFPLANE_TOLERANCE:
        on_plane = False
    else:
        on_plane = True

    # Fixed quad
    fquad = [p0.toPoint(),
             p1.toPoint(),
             new_p2.toPoint(),
             p3.toPoint()]

    return (on_plane, fquad)


def get_quad_width(q):
    """
    Return width of an individual planar trapezoid, where the p0-p1 distance
    represents the long side.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Width of planar trapezoid.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    AB = p0 - p1
    AC = p0 - p3
    t1 = (AB.cross(AC).cross(AB)).norm()
    width = t1.dot(AC)

    return width


def get_quad_mesh(q, dx):
    """
    Create a mesh from a quadrilateal.

    Args:
        q (list): A quadrilateral; list of four points.
        dx (float):  Target dx in km; used to get nx and ny of mesh, but mesh
            snaps to edges of rupture so actual dx/dy will not actually equal
            this value in general.

    Returns:
        dict: Mesh dictionary, which includes numpy arrays:

        - llx: lower left x coordinate in ECEF coords.
        - lly: lower left y coordinate in ECEF coords.
        - llz: lower left z coordinate in ECEF coords.
        - ulx: upper left x coordinate in ECEF coords.
        - etc.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Get nx based on length of top edge, minimum allowed is 2
    toplen_km = get_quad_length(q)
    nx = int(np.max([round(toplen_km / dx, 0) + 1, 2]))

    # Get array of points along top and bottom edges
    xfac = np.linspace(0, 1, nx)
    topp = [p0 + (p1 - p0) * a for a in xfac]
    botp = [p3 + (p2 - p3) * a for a in xfac]

    # Get ny based on mean length of vectors connecting top and bottom points
    ylen_km = np.ones(nx)
    for i in range(nx):
        ylen_km[i] = (topp[i] - botp[i]).mag() / 1000
    ny = int(np.max([round(np.mean(ylen_km) / dx, 0) + 1, 2]))
    yfac = np.linspace(0, 1, ny)

    # Build mesh: dict of ny by nx arrays (x, y, z):
    mesh = {'x': np.zeros([ny, nx]), 'y': np.zeros(
        [ny, nx]), 'z': np.zeros([ny, nx])}
    for i in range(nx):
        mpts = [topp[i] + (botp[i] - topp[i]) * a for a in yfac]
        mesh['x'][:, i] = [a.x for a in mpts]
        mesh['y'][:, i] = [a.y for a in mpts]
        mesh['z'][:, i] = [a.z for a in mpts]

    # Make arrays of pixel corners
    mesh['llx'] = mesh['x'][1:, 0:-1]
    mesh['lrx'] = mesh['x'][1:, 1:]
    mesh['ulx'] = mesh['x'][0:-1, 0:-1]
    mesh['urx'] = mesh['x'][0:-1, 1:]
    mesh['lly'] = mesh['y'][1:, 0:-1]
    mesh['lry'] = mesh['y'][1:, 1:]
    mesh['uly'] = mesh['y'][0:-1, 0:-1]
    mesh['ury'] = mesh['y'][0:-1, 1:]
    mesh['llz'] = mesh['z'][1:, 0:-1]
    mesh['lrz'] = mesh['z'][1:, 1:]
    mesh['ulz'] = mesh['z'][0:-1, 0:-1]
    mesh['urz'] = mesh['z'][0:-1, 1:]
    mesh['cpx'] = np.zeros_like(mesh['llx'])
    mesh['cpy'] = np.zeros_like(mesh['llx'])
    mesh['cpz'] = np.zeros_like(mesh['llx'])

    # i and j are indices over subruptures
    ni, nj = mesh['llx'].shape
    for i in range(0, ni):
        for j in range(0, nj):
            # Rupture corner points
            pp0 = Vector(
                mesh['ulx'][i, j], mesh['uly'][i, j], mesh['ulz'][i, j])
            pp1 = Vector(
                mesh['urx'][i, j], mesh['ury'][i, j], mesh['urz'][i, j])
            pp2 = Vector(
                mesh['lrx'][i, j], mesh['lry'][i, j], mesh['lrz'][i, j])
            pp3 = Vector(
                mesh['llx'][i, j], mesh['lly'][i, j], mesh['llz'][i, j])
            # Find center of quad
            mp0 = pp0 + (pp1 - pp0) * 0.5
            mp1 = pp3 + (pp2 - pp3) * 0.5
            cp = mp0 + (mp1 - mp0) * 0.5
            mesh['cpx'][i, j] = cp.x
            mesh['cpy'][i, j] = cp.y
            mesh['cpz'][i, j] = cp.z
    return mesh


def get_local_unit_slip_vector(strike, dip, rake):
    """
    Compute the components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (deg) between rupture plane and the horizontal plane
            normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike) + \
        np.cos(rake) * np.sin(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike) + \
        np.cos(rake) * np.cos(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_DS(strike, dip, rake):
    """
    Compute the DIP SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_SS(strike, dip, rake):
    """
    Compute the STRIKE SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.cos(rake) * np.sin(strike)
    sy = np.cos(rake) * np.cos(strike)
    sz = 0.0
    return Vector(sx, sy, sz)


def reverse_quad(q):
    """
    Reverse the verticies of a quad in the sense that the strike direction
    is flipped.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        list: Reversed quadrilateral.

    """
    return [q[1], q[0], q[3], q[2]]


def get_quad_slip(q, rake):
    """
    Compute the unit slip vector in ECEF space for a quad and rake angle.

    Args:
        q (list): A quadrilateral; list of four points.
        rake (float): Direction of motion of the hanging wall relative to
        the foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in ECEF space.

    """
    P0, P1, P2 = q[0:3]
    strike = P0.azimuth(P1)
    dip = get_quad_dip(q)
    s1_local = get_local_unit_slip_vector(strike, dip, rake)
    s0_local = Vector(0, 0, 0)
    qlats = [a.latitude for a in q]
    qlons = [a.longitude for a in q]
    proj = get_orthographic_projection(
        np.min(qlons), np.max(qlons), np.min(qlats), np.max(qlats))
    s1_ll = proj(np.array([s1_local.x]), np.array([s1_local.y]), reverse=True)
    s0_ll = proj(np.array([s0_local.x]), np.array([s0_local.y]), reverse=True)
    s1_ecef = Vector.fromTuple(latlon2ecef(s1_ll[1], s1_ll[0], s1_local.z))
    s0_ecef = Vector.fromTuple(latlon2ecef(s0_ll[1], s0_ll[0], s0_local.z))
    slp_ecef = (s1_ecef - s0_ecef).norm()
    return slp_ecef


def get_quad_length(q):
    """
    Length of top eduge of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Length of quadrilateral in km.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    qlength = (p1 - p0).mag() / 1000
    return qlength


def get_quad_dip(q):
    """
    Dip of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Dip in degrees.

    """
    N = get_quad_normal(q)
    V = get_vertical_vector(q)
    dip = np.degrees(np.arccos(Vector.dot(N, V)))
    return dip


def get_quad_normal(q):
    """
    Compute the unit normal vector for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized normal vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    v1 = p1 - p0
    v2 = p3 - p0
    vn = Vector.cross(v2, v1).norm()
    return vn


def get_quad_strike_vector(q):
    """
    Compute the unit vector pointing in the direction of strike for a
    quadrilateral in ECEF coordinates. Top edge assumed to be horizontal.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing in strike direction in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    v1 = (p1 - p0).norm()
    return v1


def get_quad_down_dip_vector(q):
    """
    Compute the unit vector pointing down dip for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing down dip in ECEF coords.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p0p1 = p1 - p0
    qnv = get_quad_normal(q)
    ddv = Vector.cross(p0p1, qnv).norm()
    return ddv


def get_vertical_vector(q):
    """
    Compute the vertical unit vector for a quadrilateral
    in ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized vertical vector for the quadrilateral in ECEF
                coords.
    """
    P0, P1, P2, P3 = q
    P0_up = copy.deepcopy(P0)
    P0_up.depth = P0_up.depth - 1.0
    p0 = Vector.fromPoint(P0)   # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P0_up)
    v1 = (p1 - p0).norm()
    return v1


def _quad_distance(q, points, horizontal=False):
    """
    Calculate the shortest distance from a set of points to a rupture surface.

    Args:
        q (list): A quadrilateral; list of four points.
        points (array): Numpy array Nx3 of points (ECEF) to calculate distance
            from.
        horizontal:  Boolean indicating whether to treat points inside quad as
            0 distance.

    Returns:
        float: Array of size N of distances (in km) from input points to
            rupture surface.
    """
    P0, P1, P2, P3 = q

    if horizontal:
        # project this quad to the surface
        P0 = Point(P0.x, P0.y, 0)
        P1 = Point(P1.x, P1.y, 0)
        P2 = Point(P2.x, P2.y, 0)
        P3 = Point(P3.x, P3.y, 0)

    # Convert to ecef
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Make a unit vector normal to the plane
    normalVector = (p1 - p0).cross(p2 - p0).norm()

    dist = np.ones_like(points[:, 0]) * np.nan

    p0d = p0.getArray() - points
    p1d = p1.getArray() - points
    p2d = p2.getArray() - points
    p3d = p3.getArray() - points

    # Create 4 planes with normals pointing outside rectangle
    n0 = (p1 - p0).cross(normalVector).getArray()
    n1 = (p2 - p1).cross(normalVector).getArray()
    n2 = (p3 - p2).cross(normalVector).getArray()
    n3 = (p0 - p3).cross(normalVector).getArray()

    sgn0 = np.signbit(np.sum(n0 * p0d, axis=1))
    sgn1 = np.signbit(np.sum(n1 * p1d, axis=1))
    sgn2 = np.signbit(np.sum(n2 * p2d, axis=1))
    sgn3 = np.signbit(np.sum(n3 * p3d, axis=1))

    inside_idx = (sgn0 == sgn1) & (sgn1 == sgn2) & (sgn2 == sgn3)
    if horizontal:
        dist[inside_idx] = 0.0
    else:
        dist[inside_idx] = np.power(np.abs(
            np.sum(p0d[inside_idx, :] * normalVector.getArray(), axis=1)), 2)

    outside_idx = np.logical_not(inside_idx)
    s0 = _distance_sq_to_segment(p0d, p1d)
    s1 = _distance_sq_to_segment(p1d, p2d)
    s2 = _distance_sq_to_segment(p2d, p3d)
    s3 = _distance_sq_to_segment(p3d, p0d)

    smin = np.minimum(np.minimum(s0, s1), np.minimum(s2, s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist) / 1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0], 1)
    if np.any(np.isnan(dist)):
        raise ShakeLibException("Could not calculate some distances!")
    dist = np.fliplr(dist)
    return dist


def get_distance_to_plane(planepoints, otherpoint):
    """
    Calculate a point's distance to a plane.  Used to figure out if a
    quadrilateral points are all co-planar.

    Args:
        planepoints (list): List of three points (from Vector class) defining a
            plane.
        otherpoint (Vector): 4th Vector to compare to points defining the
            plane.

    Returns:
        float: Distance (in meters) from otherpoint to plane.

    """
    # from
    # https://en.wikipedia.org/wiki/Plane_(geometry)#Describing_a_plane_through_three_points
    p0, p1, p2 = planepoints
    x1, y1, z1 = p0.getArray()
    x2, y2, z2 = p1.getArray()
    x3, y3, z3 = p2.getArray()
    D = np.linalg.det(np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]))
    if D != 0:
        d = -1
        at = np.linalg.det(np.array([[1, y1, z1], [1, y2, z2], [1, y3, z3]]))
        bt = np.linalg.det(np.array([[x1, 1, z1], [x2, 1, z2], [x3, 1, z3]]))
        ct = np.linalg.det(np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
        a = (-d / D) * at
        b = (-d / D) * bt
        c = (-d / D) * ct

        numer = np.abs(a * otherpoint.x +
                       b * otherpoint.y +
                       c * otherpoint.z + d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer / denom
    else:
        dist = 0
    return dist


def _distance_sq_to_segment(p0, p1):
    """
    Calculate the distance^2 from the origin to a segment defined by two
    vectors

    Args:
        p0 (array): Numpy array Nx3 (ECEF).
        p1 (array): Numpy array Nx3 (ECEF).

    Returns:
        float: The squared distance from the origin to a segment.
    """
    # /*
    #  * This algorithm is from (Vince's) CS1 class.
    #  * It returns the distance^2 from the origin to a segment defined
    #  * by two vectors
    #  */

    dist = np.zeros_like(p1[:, 0])
    # /* Are the two points equal? */
    idx_equal = (p0[:, 0] == p1[:, 0]) & (
        p0[:, 1] == p1[:, 1]) & (p0[:, 2] == p1[:, 2])
    dist[idx_equal] = np.sqrt(p0[idx_equal, 0]**2 +
                              p0[idx_equal, 1]**2 + p0[idx_equal, 2]**2)

    v = p1 - p0

    # /*
    #  * C1 = c1/|v| is the projection of the origin O on line (P0,P1).
    #  * If C1 is negative, then O is outside the segment and
    #  * closer to the P0 side.
    #  * If C1 is positive and >V then O is on the other side.
    #  * If C1 is positive and <V then O is inside.
    #  */

    c1 = -1 * np.sum(p0 * v, axis=1)
    idx_neg = c1 <= 0
    dist[idx_neg] = p0[idx_neg, 0]**2 + p0[idx_neg, 1]**2 + p0[idx_neg, 2]**2

    c2 = np.sum(v * v, axis=1)
    idx_less_c1 = c2 <= c1
    dist[idx_less_c1] = p1[idx_less_c1, 0]**2 + \
        p1[idx_less_c1, 1]**2 + p1[idx_less_c1, 2]**2

    idx_other = np.logical_not(idx_neg | idx_equal | idx_less_c1)

    nr, nc = p0.shape
    t1 = c1 / c2
    t1.shape = (nr, 1)
    tmp = p0 + (v * t1)
    dist[idx_other] = tmp[idx_other, 0]**2 + \
        tmp[idx_other, 1]**2 + tmp[idx_other, 2]**2

    return dist


def rake_to_mech(rake):
    """
    Convert rake to mechanism.

    Args:
        rake (float): Rake angle in degrees.

    Returns:
        str: Mechanism.
    """
    mech = 'ALL'
    if rake is not None:
        if (rake >= -180 and rake <= -150) or \
           (rake >= -30 and rake <= 30) or \
           (rake >= 150 and rake <= 180):
            mech = 'SS'
        if rake >= -120 and rake <= -60:
            mech = 'NM'
        if rake >= 60 and rake <= 120:
            mech = 'RS'
    return mech
