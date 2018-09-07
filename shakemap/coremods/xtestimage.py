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
            datadict[myimt] = container.getIMTGrids(
                myimt, 'GREATER_OF_TWO_HORIZONTAL')

        container.close()

        #
        # Make plots
        #
        for myimt in imtlist:
            if myimt == 'MMI':
                yunits = '(MMI)'
            elif myimt == 'PGV':
                yunits = '[ln(cm/s)]'
            else:
                yunits = '[ln(g)]'

            fileimt = oq_to_file(myimt)

            #
            # Do the ground motion plots
            #
            data = datadict[myimt]
            grddata = data['mean']
            metadata = data['mean_metadata']

            fig = plt.figure(figsize=(10, 10))
            gs = plt.GridSpec(4, 4, hspace=0.2, wspace=0.1)
            ax0 = fig.add_subplot(gs[:-1, 1:])
            plt.title(self._eventid + ': ' + myimt + ' mean')
            im1 = ax0.imshow(grddata,
                             extent=(metadata['xmin'], metadata['xmax'],
                                     metadata['ymin'], metadata['ymax']))
            cbax = fig.add_axes([0.915, .34, .02, .5])
            plt.colorbar(im1, ax=ax0, cax=cbax)
            ycut = fig.add_subplot(gs[:-1, 0], sharey=ax0)
            xcut = fig.add_subplot(gs[-1, 1:], sharex=ax0)
            rows, cols = grddata.shape
            midrow = int(rows / 2)
            midcol = int(cols / 2)
            xvals = np.linspace(metadata['xmin'], metadata['xmax'], cols)
            yvals = np.linspace(metadata['ymin'], metadata['ymax'], rows)
            ycut.plot(grddata[:, midcol], yvals)
            xcut.plot(xvals, grddata[midrow, :])
            ycut.set(xlabel=myimt + ' ' + yunits, ylabel='Latitude')
            xcut.set(xlabel='Longitude', ylabel=myimt + ' ' + yunits)
            ycut.set_ylim((metadata['ymin'], metadata['ymax']))
            xcut.set_xlim((metadata['xmin'], metadata['xmax']))
            ax0.label_outer()

            pfile = os.path.join(datadir,
                                 self._eventid + '_' + fileimt + '.pdf')
            plt.savefig(pfile, bbox_inches='tight')
            plt.close()

            #
            # Do the stddev plots
            #
            grddata = data['std']

            fig = plt.figure(figsize=(10, 10))
            gs = plt.GridSpec(4, 4, hspace=0.2, wspace=0.1)
            ax0 = fig.add_subplot(gs[:-1, 1:])
            plt.title(self._eventid + ': ' + myimt + ' stddev')
            im1 = ax0.imshow(grddata,
                             extent=(metadata['xmin'], metadata['xmax'],
                                     metadata['ymin'], metadata['ymax']))
            cbax = fig.add_axes([0.915, .34, .02, .5])
            plt.colorbar(im1, ax=ax0, cax=cbax)
            ycut = fig.add_subplot(gs[:-1, 0], sharey=ax0)
            xcut = fig.add_subplot(gs[-1, 1:], sharex=ax0)
            rows, cols = grddata.shape
            midrow = int(rows / 2)
            midcol = int(cols / 2)
            xvals = np.linspace(metadata['xmin'], metadata['xmax'], cols)
            yvals = np.linspace(metadata['ymin'], metadata['ymax'], rows)
            ycut.plot(grddata[:, midcol], yvals)
            xcut.plot(xvals, grddata[midrow, :])
            ycut.set(xlabel='stddev ' + yunits, ylabel='Latitude')
            xcut.set(xlabel='Longitude', ylabel='stddev ' + yunits)
            xcut.set_xlim((metadata['xmin'], metadata['xmax']))
            xcut.set_ylim(bottom=0, top=np.max(grddata[midrow, :]) * 1.1)
            ycut.set_xlim(left=0, right=np.max(grddata[:, midcol] * 1.1))
            ycut.set_ylim((metadata['ymin'], metadata['ymax']))
            ax0.label_outer()

            pfile = os.path.join(datadir,
                                 self._eventid + '_' + fileimt + '_sd.pdf')
            plt.savefig(pfile, bbox_inches='tight')
            plt.close()
