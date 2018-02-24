# stdlib
import os.path

# third party

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

# neic imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from shakemap.utils.config import get_config_paths
from .base import CoreModule
from shakelib.utils.imt_string import oq_to_file


class XTestImage(CoreModule):
    """
    xtestimage -- Plot 2D images of ShakeMap arrays
    """

    command_name = 'xtestimage'

    def execute(self):
        """
        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)
        if container.getDataType() != 'grid':
            raise NotImplementedError('xtestimage module can only operate on '
                                      'gridded data not sets of points')

        datadict = {}
        imtlist = container.getIMTs('GREATER_OF_TWO_HORIZONTAL')
        for myimt in imtlist:
            datadict[myimt] = container.getIMTGrids(myimt,
                                                    'GREATER_OF_TWO_HORIZONTAL')

        #
        # Make plots
        #
        for myimt in imtlist:
            fileimt = oq_to_file(myimt)

            data = datadict[myimt]
            gridobj = data['mean']
            grddata = gridobj.getData()
            metadata = gridobj.getGeoDict().asDict()

            fig = plt.figure(figsize=(10, 10))
            gs = plt.GridSpec(4, 4, hspace=0.1, wspace=0.1)
            ax0 = fig.add_subplot(gs[:-1, 1:])
            plt.title(self._eventid + ': ' + myimt + ' mean')
            ax0.imshow(grddata,
                       extent=(metadata['xmin'], metadata['xmax'],
                               metadata['ymin'], metadata['ymax']))
            ycut = fig.add_subplot(gs[:-1, 0], sharey=ax0)
            xcut = fig.add_subplot(gs[-1, 1:], sharex=ax0)
            rows, cols = grddata.shape
            midrow = int(rows / 2)
            midcol = int(cols / 2)
            xvals = np.linspace(metadata['xmin'], metadata['xmax'], cols)
            yvals = np.linspace(metadata['ymin'], metadata['ymax'], rows)
            ycut.plot(grddata[:, midcol], yvals)
            xcut.plot(xvals, grddata[midrow, :])
            ycut.set(ylabel='Latitude')
            xcut.set(xlabel='Longitude')
            ax0.label_outer()

            pfile = os.path.join(datadir, self._eventid + '_' + fileimt + '.pdf')
            plt.savefig(pfile)
            plt.close()

            gridobj = data['std']
            grddata = gridobj.getData()

            fig = plt.figure(figsize=(10, 10))
            gs = plt.GridSpec(4, 4, hspace=0.1, wspace=0.1)
            ax0 = fig.add_subplot(gs[:-1, 1:])
            plt.title(self._eventid + ': ' + myimt + ' stddev')
            ax0.imshow(grddata,
                       extent=(metadata['xmin'], metadata['xmax'],
                               metadata['ymin'], metadata['ymax']))
            ycut = fig.add_subplot(gs[:-1, 0], sharey=ax0)
            xcut = fig.add_subplot(gs[-1, 1:], sharex=ax0)
            rows, cols = grddata.shape
            midrow = int(rows / 2)
            midcol = int(cols / 2)
            xvals = np.linspace(metadata['xmin'], metadata['xmax'], cols)
            yvals = np.linspace(metadata['ymin'], metadata['ymax'], rows)
            ycut.plot(grddata[:, midcol], yvals)
            xcut.plot(xvals, grddata[midrow, :])
            ycut.set(ylabel='Latitude')
            xcut.set(xlabel='Longitude')
            ax0.label_outer()

            pfile = os.path.join(datadir, self._eventid + '_' + fileimt + '_sd.pdf')
            plt.savefig(pfile)
            plt.close()
