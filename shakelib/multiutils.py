import numpy as np

from openquake.hazardlib import const
from openquake.hazardlib.contexts import ContextMaker


def gmpe_gmas(gmpe, ctx, imt, stddev_types):
    """ """
    N = len(ctx)
    mean = np.zeros((1, N))
    sig = np.zeros((1, N))
    tau = np.zeros((1, N))
    phi = np.zeros((1, N))
    if gmpe.compute.__annotations__.get("ctx") is np.recarray:
        if isinstance(ctx.mag, np.ndarray):
            magstr = "%.2f" % ctx.mag[0]
        else:
            magstr = "%.2f" % ctx.mag
        param = dict(
            imtls={imt.string: [0]},
            maximum_distance=4000,
            truncation_level=3,
            investigation_time=1.0,
            mags=[magstr],
        )
        cmaker = ContextMaker("*", [gmpe], param)
        if not isinstance(ctx, np.ndarray):
            ctx = cmaker.recarray([ctx])
    try:
        gmpe.compute(ctx, [imt], mean, sig, tau, phi)
    except NotImplementedError:
        mean, stddevs = gmpe.get_mean_and_stddevs(ctx, ctx, ctx, imt, stddev_types)
        return mean, stddevs
    except:
        raise
    else:
        stddevs = []
        for stddev_type in stddev_types:
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(sig[0])
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(tau[0])
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(phi[0])
        return mean[0], stddevs
