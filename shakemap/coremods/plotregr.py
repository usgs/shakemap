# stdlib
import os.path
import concurrent.futures as cf
import json

# third party
from configobj import ConfigObj
import matplotlib.pyplot as plt
import numpy as np

# neic imports
from impactutils.io.smcontainers import ShakeMapOutputContainer
from openquake.hazardlib import imt

# local imports
from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error)
from .base import CoreModule, Contents
from shakelib.utils.imt_string import oq_to_file


class PlotRegr(CoreModule):
    """
    plotregr -- Plot the regression curves from a data file
    """

    command_name = 'plotregr'
    targets = [r'products/.*_regr\.png']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(PlotRegr, self).__init__(eventid)
        self.contents = Contents('Regression Plots', 'regression', eventid)

    def execute(self):
        """
        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        oc = ShakeMapOutputContainer.load(datafile)
        if oc.getDataType() != 'grid':
            raise NotImplementedError('plotregr module can only operate on '
                                      'gridded data not sets of points')

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        config = ConfigObj(config_file, configspec=spec_file)
        results = config.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        # If mapping runs in parallel, then we want this module too, as well.
        # Otherwise we get weird errors from matplotlib
        max_workers = config['products']['mapping']['max_workers']

        components = oc.getComponents()
        for component in components:
            #
            # Cheating here a bit by assuming that the IMTs are the same
            # as the regression IMTs
            #
            rockgrid = {}
            soilgrid = {}
            rocksd = {}
            soilsd = {}
            imtlist = oc.getIMTs(component)
            for myimt in imtlist:
                rockgrid[myimt], _ = oc.getArray(['attenuation', 'rock',
                                                  myimt], 'mean')
                soilgrid[myimt], _ = oc.getArray(['attenuation', 'soil',
                                                  myimt], 'mean')
                rocksd[myimt], _ = oc.getArray(['attenuation', 'rock',
                                                myimt], 'std')
                soilsd[myimt], _ = oc.getArray(['attenuation', 'soil',
                                                myimt], 'std')
            distances, _ = oc.getArray(['attenuation', 'distances'], 'rrup')

            stations = oc.getStationDict()

            if component == 'GREATER_OF_TWO_HORIZONTAL':
                ctype = ''
                cfile = ''
            else:
                ctype = ' %s' % component
                cfile = '_%s' % component
            #
            # Make plots
            #
            alist = []
            for myimt in imtlist:
                a = {'myimt': myimt,
                     'rockgrid': rockgrid,
                     'soilgrid': soilgrid,
                     'rocksd': rocksd,
                     'soilsd': soilsd,
                     'stations': stations,
                     'distances': distances,
                     'eventid': self._eventid,
                     'datadir': datadir,
                     'component': cfile
                     }
                alist.append(a)
                if myimt == 'MMI':
                    self.contents.addFile('miRegr', 'Intensity Regression',
                                          'Regression plot of macroseismic '
                                          'intensity%s.' % ctype,
                                          'mmi_regr%s.png' % cfile,
                                          'image/png')
                elif myimt == 'PGA':
                    self.contents.addFile(
                        'pgaRegr', 'PGA Regression', 'Regression plot of peak '
                        'ground acceleration (%%g)%s.' % ctype,
                        'pga_regr%s.png' % cfile,
                        'image/png')
                elif myimt == 'PGV':
                    self.contents.addFile('pgvRegr', 'PGV Regression',
                                          'Regression plot of peak ground '
                                          'velocity (cm/s)%s.' % ctype,
                                          'pgv_regr%s.png' % cfile,
                                          'image/png')
                else:
                    oqimt = imt.from_string(myimt)
                    period = str(oqimt.period)
                    filebase = oq_to_file(myimt)
                    psacap = ('Regression plot of ' + period +
                              ' sec 5%% damped pseudo-spectral '
                              'acceleration(%%g)%s.' % ctype)
                    self.contents.addFile(filebase + 'Regr',
                                          'PSA ' + period + ' sec Regression',
                                          psacap,
                                          filebase + '_regr%s.png' % cfile,
                                          'image/png')

            if max_workers > 0:
                with cf.ProcessPoolExecutor(max_workers=max_workers) as ex:
                    results = ex.map(make_plots, alist)
                    list(results)
            else:
                for adict in alist:
                    make_plots(adict)

            #
            # Make attenuation_curves.json
            #
            jdict = {'eventid': self._eventid}
            jdict['gmpe'] = {}
            for site in ['soil', 'rock']:
                jdict['gmpe'][site] = {}
                for myimt in imtlist:
                    jdict['gmpe'][site][myimt] = {}
                    jdict['gmpe'][site][myimt]['mean'] = oc.getArray(
                        ['attenuation', site, myimt], 'mean')[0].round(
                            decimals=5).tolist()
                    jdict['gmpe'][site][myimt]['stddev'] = oc.getArray(
                        ['attenuation', site, myimt], 'std')[0].round(
                            decimals=5).tolist()
            jdict['distances'] = {}
            for dtype in ['repi', 'rhypo', 'rjb', 'rrup']:
                jdict['distances'][dtype] = oc.getArray(
                    ['attenuation', 'distances'], dtype)[0].round(
                        decimals=5).tolist()
            jdict['mean_bias'] = {}
            info = oc.getMetadata()
            for myimt in imtlist:
                jdict['mean_bias'][myimt] = info['output']['ground_motions'][
                                                 myimt]['bias']
            jstring = json.dumps(jdict, allow_nan=False)
            jfile = os.path.join(datadir, 'attenuation_curves%s.json' % cfile)
            f = open(jfile, 'wt')
            f.write(jstring)
            f.close()
            oc.close()
            cap = "Nominal attenuation curves"
            self.contents.addFile('attenuationCurves', 'Attenuation Curves',
                                  cap, jfile, 'application/json')


def make_plots(adict):
    myimt = adict['myimt']
    eventid = adict['eventid']
    datadir = adict['datadir']
    rockgrid = adict['rockgrid']
    soilgrid = adict['soilgrid']
    rocksd = adict['rocksd']
    soilsd = adict['soilsd']
    stations = adict['stations']
    distances = adict['distances']
    cfile = adict['component']

    plt.figure(figsize=(10, 10))

    plt.semilogx(distances, rockgrid[myimt], 'r', label='rock')
    plt.semilogx(distances, soilgrid[myimt], 'g', label='soil')
    plt.semilogx(distances, rockgrid[myimt] + rocksd[myimt], 'r--',
                 label='rock +/- stddev')
    plt.semilogx(distances, rockgrid[myimt] - rocksd[myimt], 'r--')
    plt.semilogx(distances, soilgrid[myimt] + soilsd[myimt], 'g--',
                 label='soil +/- stddev')
    plt.semilogx(distances, soilgrid[myimt] - soilsd[myimt], 'g--')

    for station in stations['features']:
        dist = station['properties']['distances']['rrup']
        if dist > distances[-1]:
            continue
        if station['properties']['station_type'] == 'seismic':
            symbol = '^'
            if myimt == 'MMI':
                value = station['properties']['intensity']
                if value != 'null':
                    plt.semilogx(dist, value, symbol + 'k', mfc='none')
            else:
                imtstr = myimt.lower()
                value = np.nan
                for chan in station['properties']['channels']:
                    if chan['name'].endswith('Z') or \
                       chan['name'].endswith('U'):
                        continue
                    for amp in chan['amplitudes']:
                        if amp['name'] != imtstr:
                            continue
                        if amp['flag'] != '' and amp['flag'] != '0':
                            break
                        if amp['value'] is None or \
                                amp['value'] == 'null':
                            break
                        if isinstance(amp['value'], str):
                            thisamp = float(amp['value'])
                        else:
                            thisamp = amp['value']
                        if thisamp <= 0:
                            break
                        if myimt == 'PGV':
                            tmpval = np.log(thisamp)
                        else:
                            tmpval = np.log(thisamp / 100.)
                        if np.isnan(value) or tmpval > value:
                            value = tmpval
                        break
                if not np.isnan(value):
                    plt.semilogx(dist, value, symbol + 'k', mfc='none')
        else:
            symbol = 'o'
            if myimt == 'MMI':
                amp = station['properties']['intensity']
                flag = station['properties']['intensity_flag']
                if flag == '' or flag == '0':
                    if amp is not None and amp != 'null':
                        if isinstance(amp, str):
                            value = float(amp)
                        else:
                            value = amp
                        plt.semilogx(dist, value, symbol + 'k',
                                     mfc='none')
            else:
                imtstr = myimt.lower()
                for thing in station['properties']['pgm_from_mmi']:
                    if thing['name'] != imtstr:
                        continue
                    amp = thing['value']
                    if amp is not None and amp != 'null' and amp != 0:
                        if myimt == 'PGV':
                            amp = np.log(amp)
                        else:
                            amp = np.log(amp / 100.)
                        plt.semilogx(dist, amp, symbol + 'k',
                                     mfc='none')
                    break

    plt.title(eventid + ': ' + myimt + ' mean')
    plt.xlabel('Rrup (km)')
    if myimt == 'MMI':
        plt.ylabel('MMI')
    elif myimt == 'PGV':
        plt.ylabel('PGV ln(cm/s)')
    else:
        plt.ylabel(myimt + ' ln(g)')
    plt.legend()

    fileimt = oq_to_file(myimt)
    pfile = os.path.join(datadir, fileimt + '_regr%s.png' % cfile)
    plt.savefig(pfile)
    plt.close()
