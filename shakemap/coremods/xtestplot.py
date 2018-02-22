# stdlib
import os.path

# third party

import matplotlib.pyplot as plt

# neic imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from shakemap.utils.config import get_config_paths
from .base import CoreModule
from shakelib.utils.imt_string import oq_to_file


class XTestPlot(CoreModule):
    """
    xtestplot -- Plot 1D sections of test events.
    """

    command_name = 'xtestplot'

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

        datadict = {}
        imtlist = container.getIMTs('GREATER_OF_TWO_HORIZONTAL')
        for myimt in imtlist:
            datadict[myimt] = container.getIMTArrays(myimt,
                                                     'GREATER_OF_TWO_HORIZONTAL')

        #
        # Make plots
        #
        for myimt in imtlist:
            data = datadict[myimt]
            fig = plt.figure(figsize=(10, 8))
            plt.plot(data['lons'],
                     data['mean'],
                     color='k', label='mean')
            plt.plot(data['lons'],
                     data['mean'] + data['std'],
                     '--b', label='mean +/- stddev')
            plt.plot(data['lons'],
                     data['mean'] - data['std'],
                     '--b')
            plt.plot(data['lons'],
                     data['std'],
                     '-.r', label='stddev')
            plt.xlabel('Longitude')
            plt.ylabel('ln(%s) (g)' % myimt)
            plt.legend(loc='best')
            plt.title(self._eventid)
            plt.grid()
            filleimt = oq_to_file(myimt)
            pfile = os.path.join(datadir, self._eventid + '_' + fileimt + '.pdf')
            plt.savefig(pfile)
            plt.close()
