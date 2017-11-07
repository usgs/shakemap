def get_period_index_from_imt_str(imtstr, imt_per_ix):
        if imtstr == 'PGA':
            return imt_per_ix['0.01']
        elif imtstr in ('PGV', 'MMI'):
            return imt_per_ix['1.0']
        else:
            return imt_per_ix[imtstr.replace('SA(', '').replace(')', '')]

# TODO Need to doc all of these functions
def get_period_array(*args):
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
    myimts = set()
    for this_imt in imtset:
        if not np.isnan(sdf[this_imt][ix]) and \
           not sdf[this_imt + '_outliers'][ix]:
            myimts.add(this_imt)
    return get_imts(imtstr, myimts)

def get_map_grade(do_grid, outsd, psd, moutgrid):
    mean_rat = '-'
    mygrade = '-'
    if not do_grid or 'PGA' not in outsd:
        return mean_rat, mygrade, []
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
    return mean_rat, mygrade, sd_rat

