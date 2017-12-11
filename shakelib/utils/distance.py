import numexpr as ne

EARTH_RADIUS = 6371.0

def geodetic_distance_fast(lons1, lats1, lons2, lats2):
    """
    Calculate the geodetic distance between two points or two collections
    of points using a formula that is substantially faster than the
    Haversine formula, but is nearly as accurate for distances up to a
    few hundred km.

    Parameters are coordinates in RADIANS. They could be scalar
    float numbers or numpy arrays, in which case they should be either
    the same dimensions (in which case the returned distances will be
    an array of the same shape with point for point distances), or they
    should broadcast together (in which case the returned distances will
    be a matrix of each point in the first set to every point in the
    second set.

    Args:
        lons1 (numpy array): An array of longitudes.
        lats1 (numpy array): An array of latitudes the same shape as lons1.
        lons2 (numpy array): An array of longitudes.
        lats2 (numpy array): An array of latitudes the same shape as lons2.

    Returns:
        (numpy array): Distances in km.
    """
    d = ne.evaluate("EARTH_RADIUS * sqrt(((lons1 - lons2) * cos(0.5 * "
                    "(lats1 + lats2)))**2.0 + (lats1 - lats2)**2.0)")
    return d


def geodetic_distance_haversine(lons1, lats1, lons2, lats2):
    """
    Calculate the geodetic distance between two points or two collections
    of points using the Haversine formula.

    Parameters are coordinates in RADIANS. They could be scalar
    float numbers or numpy arrays, in which case they should be either
    the same dimensions (in which case the returned distances will be
    an array of the same shape with point for point distances), or they
    should broadcast together (in which case the returned distances will
    be a matrix of each point in the first set to every point in the
    second set.

    Args:
        lons1 (numpy array): An array of longitudes.
        lats1 (numpy array): An array of latitudes the same shape as lons1.
        lons2 (numpy array): An array of longitudes.
        lats2 (numpy array): An array of latitudes the same shape as lons2.

    Returns:
        (numpy array): Distances in km.
    """
    diameter = 2 * EARTH_RADIUS
    distance = ne.evaluate(
            "diameter * arcsin(sqrt(sin((lats1 - lats2) / 2.0)**2.0 "
            "+ cos(lats1) * cos(lats2) * sin((lons1 - lons2) / 2.0)**2.0))")
    return distance

