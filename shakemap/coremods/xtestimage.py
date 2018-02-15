# stdlib
import os.path

# third party

import matplotlib.pyplot as plt

# neic imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from shakemap.utils.config import get_config_paths
from .base import CoreModule


class TestPlot(CoreModule):
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
            data = datadict[myimt]
            gridobj = data['mean']
            metadata = gridobj.getGeoDict().asDict()
            plt.figure(figsize=(10, 10))
            plt.imshow(gridobj.getData(),
                       extent=(metadata['xmin'], metadata['xmax'],
                               metadata['ymin'], metadata['ymax']))
            plt.colorbar(shrink=0.6)
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.title(self._eventid + ': ' + myimt + ' mean')
            pfile = os.path.join(datadir, self._eventid + '_' + myimt + '.pdf')
            plt.savefig(pfile)

            gridobj = data['std']
            plt.figure(figsize=(10, 10))
            plt.imshow(gridobj.getData(),
                       extent=(metadata['xmin'], metadata['xmax'],
                               metadata['ymin'], metadata['ymax']))
            plt.colorbar(shrink=0.6)
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.title(self._eventid + ': ' + myimt + ' stddev')
            pfile = os.path.join(datadir, self._eventid + '_' + myimt + '_sd.pdf')
            plt.savefig(pfile)
