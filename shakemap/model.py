import re

import numpy as np
import numpy.ma as ma
"""
Functions to support sm_model.
"""


def get_period_index_from_imt_str(imtstr, imt_per_ix):
        """
        Get the index for the period of the specified IMT.

        Args:
            imtstr (str): The (OQ-style) IMT string.
            imt_per_ix (dict): Dictionary relating periods to their
                indices.

        Returns:
            int: The index corresponding to the period of the IMT.
        """
        if imtstr == 'PGA':
            return imt_per_ix['0.01']
        elif imtstr in ('PGV', 'MMI'):
            return imt_per_ix['1.0']
        else:
            return imt_per_ix[imtstr.replace('SA(', '').replace(')', '')]


def get_period_array(*args):
    """
    Return an array of the periods represented by the IMT list(s) in
    the input.

    Args:
        *args (list): One or more lists of IMTs.

    Returns:
        array: Numpy array of the (sorted) periods represented by the
        IMTs in the input list(s).
    """
    imt_per = set()
    for imt_list in args:
        if imt_list is None:
            continue
        for imtstr in imt_list:
            if imtstr == 'PGA':
                imt_per.add(0.01)
            elif imtstr == 'PGV' or imtstr == 'MMI':
                imt_per.add(1.0)
            else:
                imt_per.add(float(imtstr.replace('SA(', '').replace(')', '')))
    return np.array(sorted(imt_per))


def get_imts(imtstr, imtset):
    """
    Return the input imt or its closest surrogarte (or bracket) found
    in imtset.

    Args:
        imtstr (str): An (OQ-style) IMT string.
        imtsset (set): A set of IMTs to search for imtstr or its closest
            surrogate (or bracket).

    Returns:
        tuple: The IMT, it's closest surrogate, or a bracket of periods on
        either side of the IMT's period, from the IMTs in intset.
    """

    if imtstr in imtset:
        return (imtstr, )

    salist = [x for x in imtset if x.startswith('SA(')]
    periodlist = [float(x.replace('SA(', '').replace(')', '')) for x in salist]
    periodlist = sorted(periodlist)
    periodlist_str = [str(x) for x in periodlist]

    #
    # If we're here, then we know that IMT isn't in the inputs. Try
    # some alternatives.
    #
    if imtstr == 'PGA':
        #
        # Use the highest frequency in the inputs, otherwise use PGV
        #
        if len(salist):
            return ('SA(' + periodlist_str[0] + ')', )
        elif 'PGV' in imtset:
            return ('PGV', )
        else:
            return ()
    elif imtstr == 'PGV':
        #
        # Use 1.0 sec SA (or its bracket) if it's there, otherwise
        # use PGA
        #
        sa_tuple = get_sa_bracket(1.0, periodlist, periodlist_str)
        if sa_tuple != ():
            return sa_tuple
        if 'PGA' in imtset:
            return ('PGA', )
        else:
            return ()
    elif imtstr == 'MMI':
        #
        # Use PGV if it's there, otherwise use 1.0 sec SA (or its
        # bracket)
        #
        if 'PGV' in imtset:
            return ('PGV', )
        return get_sa_bracket(1.0, periodlist, periodlist_str)
    elif imtstr.startswith('SA('):
        myper = float(imtstr.replace('SA(', '').replace(')', ''))
        return get_sa_bracket(myper, periodlist, periodlist_str)
    else:
        raise ValueError('Unknown IMT %s in get_imt_bracket' % imtstr)


def get_sa_bracket(myper, plist, plist_str):

    """
    For a given period, look through the input periods and return a tuple of
    a) the single IMT string representing the period itself if it is found
    in the input; b) a pair of IMT strings representing the periods
    bracketing the given period; or c) the single IMT representing the first
    or last period in the input list if the given period is off the end of
    the list.

    Args:
        myper (float): The period to search for in the input lists.
        plist (array): A sorted list of periods as floats.
        plist_str (array): The periods in 'plist' (above) as strings.

    Returns:
        tuple: One or two strings representing the IMTs closest to or
        bracketing the input period.

    """
    if not len(plist):
        return ()
    try:
        return ('SA(' + plist_str[plist.index(myper)] + ')', )
    except ValueError:
        pass
    for i, v in enumerate(plist):
        if v > myper:
            break
    if i == 0 or v < myper:
        return ('SA(' + plist_str[i] + ')', )
    else:
        return ('SA(' + plist_str[i-1] + ')', 'SA(' + plist_str[i] + ')')


def get_sta_imts(imtstr, sdf, ix, imtset):
    """
    Get the IMT, its closest surrogate, or its bracket for a stataion
    in a particular station dataframe, accounting for missing or
    flagged data.

    Args:
        imtstr (str): The desired IMT as an OQ-style string.

        sdf (dict): The dataframe containing the station.

        ix (int): The index of the station in the dataframe.

        imtset (set): The list of IMTs (as OQ-style strings) in the
            dataframe.

    Returns:
        tuple: The IMT, its closest surrogate, or the pair of IMTs
        bracketing it in period, gathered from the valid data for the
        station.
    """
    myimts = set()
    for this_imt in imtset:
        if not np.isnan(sdf[this_imt][ix]) and \
           not sdf[this_imt + '_outliers'][ix]:
            myimts.add(this_imt)
    return get_imts(imtstr, myimts)


def get_map_grade(do_grid, outsd, psd, moutgrid):
    """
    Computes a 'grade' for the map. Essentially looks at the ratio of
    the computed PGA uncertainty to the predicted PGA uncertainty for
    the area inside the MMI 6 contour. If the maximum MMI is less than
    6, or the map is not a grid, the grade and mean ratio are set to '-'.

    Args:
        do_grid (bool): Is the map a grid (True) or a list of points
            (False)?

        outsd (dict): A dictionary of computed uncertainty arrays.

        psd (dict): A dictionary of predicted uncertainty arrays.

        moutgrid (dict): A dictionary of landmasked output ground
            motion arrays.

    Returns:
        tuple: The mean uncertainty ratio and the letter grade
    """
    mean_rat = '-'
    mygrade = '-'
    if not do_grid or 'PGA' not in outsd or 'PGA' not in psd:
        return mean_rat, mygrade
    sd_rat = outsd['PGA'] / psd['PGA']
    mmimasked = ma.masked_less(moutgrid['MMI'], 6.0)
    mpgasd_rat = ma.masked_array(sd_rat, mask=mmimasked.mask)
    if not np.all(mpgasd_rat.mask):
        gvals = [0.96, 0.98, 1.05, 1.25]
        grades = ['A', 'B', 'C', 'D', 'F']
        mean_rat = mpgasd_rat.mean()
        for ix, val in enumerate(gvals):
            if mean_rat <= val:
                mygrade = grades[ix]
                break
        if mygrade == '-':
            mygrade = 'F'
    return mean_rat, mygrade


# we need a way to get units information about intensity measure types
# and translate between openquake naming convention and ShakeMap grid naming
# convention.
def get_layer_info(layer):
    layer_out = layer
    layer_units = ''
    layer_digits = 4  # number of significant digits
    if layer.startswith('SA'):
        match = re.match(r'^SA\(([^)]+?)\)', layer)
        period = float(match.group(1))
        layer_out = ('PSA%g' % period).replace('.', 'p')
        if layer.endswith('_sd'):
            layer_out = layer_out.replace('$', '_sd')
        layer_units = 'ln(g)'
    elif layer.startswith('PGA'):
        layer_units = 'ln(g)'
    elif layer.startswith('PGV'):
        layer_units = 'ln(cm/s)'
    elif layer.startswith('MMI'):
        layer_units = 'intensity'
        layer_digits = 2
    elif layer.startswith('vs30'):
        layer_units = 'm/s'
    else:
        raise ValueError('Unknown layer type: %s' % layer)

    return (layer_out, layer_units, layer_digits)


#
# Helper function to call get_mean_and_stddevs for the
# appropriate object given the IMT
#
def gmas(ipe, gmpe, sx, rx, dx, oqimt, stddev_types):
    if 'MMI' in oqimt:
        pe = ipe
    else:
        pe = gmpe
    return pe.get_mean_and_stddevs(sx, rx, dx, oqimt, stddev_types)


#
# This is a stand-in for tau when the gmpe set doesn't provide it
# It is an educated guess based on the NGA-east work and BC Hydro
# gmpe. It's not perfect, but probably isn't too far off. It is
# only used when the GMPEs don't provide a breakdown of the
# uncertainty terms.
#
def get_default_std_inter():
    return 0.35
