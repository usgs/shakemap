# stdlib
import os.path

# third party

import matplotlib.pyplot as plt

import numpy as np

# neic imports
from impactutils.io.smcontainers import ShakeMapOutputContainer

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

        container.close()

        #
        # Make plots
        #
        fig, axa = plt.subplots(2, sharex=True, figsize=(10, 8))
        plt.subplots_adjust(hspace=0.1)
        axa[0].semilogx(periodlist[indxx],
                        datalist[indxx],
                        color='k', label='mean')
        axa[0].semilogx(periodlist[indxx],
                        datalist[indxx] + stddevlist[indxx],
                        '--b', label='mean +/- stddev')
        axa[0].semilogx(periodlist[indxx],
                        datalist[indxx] - stddevlist[indxx],
                        '--b')
        axa[1].semilogx(periodlist[indxx],
                        stddevlist[indxx],
                        '-.r', label='stddev')
        axa[1].set_xlabel('Period (s)')
        axa[0].set_ylabel('Mean ln(SA) (g)')
        axa[1].set_ylabel('Stddev ln(SA) (g)')
        axa[0].legend(loc='best')
        axa[1].legend(loc='best')
        axa[0].set_title(self._eventid)
        axa[0].grid()
        axa[1].grid()
        axa[1].set_ylim(bottom=0)
        pfile = os.path.join(datadir, self._eventid + '_spectra_plot.pdf')
        plt.savefig(pfile)
        pfile = os.path.join(datadir, self._eventid + '_spectra_plot.png')
        plt.savefig(pfile)
        plt.close()
