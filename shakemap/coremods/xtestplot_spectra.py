# stdlib
import os.path

# third party

import matplotlib.pyplot as plt

import numpy as np

# neic imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from shakemap.utils.config import get_config_paths
from .base import CoreModule


class XTestPlotSpectra(CoreModule):
    """
    xtestplot_spectra -- Plot spectra of test events.
    """

    command_name = 'xtestplot_spectra'

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
        container = ShakeMapOutputContainer.load(datafile)
        if container.getDataType() != 'points':
            raise NotImplementedError('xtestplot module can only operate on '
                                      'sets of points, not gridded data')

        datalist = []
        stddevlist = []
        periodlist = []
        imtlist = container.getIMTs('GREATER_OF_TWO_HORIZONTAL')
        for myimt in imtlist:
            if not myimt.startswith('SA('):
                continue
            ddict = container.getIMTArrays(myimt, 'GREATER_OF_TWO_HORIZONTAL')
            datalist.append(ddict['mean'][0])
            stddevlist.append(ddict['std'][0])
            periodlist.append(float(myimt.replace('SA(', '').replace(')', '')))
            self.logger.debug(myimt, datalist[-1])
        datalist = np.array(datalist)
        stddevlist = np.array(stddevlist)
        periodlist = np.array(periodlist)
        indxx = np.argsort(periodlist)

        #
        # Make plots
        #
        fig = plt.figure(figsize=(10, 8))
        plt.semilogx(periodlist[indxx],
                     datalist[indxx],
                     color='k', label='mean')
        plt.semilogx(periodlist[indxx],
                     datalist[indxx] + stddevlist[indxx],
                     '--b', label='mean +/- stddev')
        plt.semilogx(periodlist[indxx],
                     datalist[indxx] - stddevlist[indxx],
                     '--b')
        plt.semilogx(periodlist[indxx],
                     stddevlist[indxx],
                     '-.r', label='stddev')
        plt.xlabel('Period (s)')
        plt.ylabel('ln(SA) (g)')
        plt.legend(loc='best')
        plt.title(self._eventid)
        plt.grid()
        pfile = os.path.join(datadir, self._eventid + '_spectra_plot.pdf')
        plt.savefig(pfile)
        plt.close()
