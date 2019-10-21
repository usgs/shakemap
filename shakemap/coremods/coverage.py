# stdlib imports
import os.path
import json

# third party imports
import numpy as np
from scipy.ndimage import gaussian_filter
from impactutils.io.smcontainers import ShakeMapOutputContainer
from openquake.hazardlib import imt

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths
from shakelib.utils.imt_string import oq_to_file

# Not really relevant, but seemingly necessary
COMPONENT = 'GREATER_OF_TWO_HORIZONTAL'


class CoverageModule(CoreModule):
    """
    coverage -- Create JSON coverage(s) of the ground motion layers.
    """

    command_name = 'coverage'
    targets = [r'products/coverage_h\.json', r'products/coverage_m\.json',
               r'products/coverage_l\.json']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(CoverageModule, self).__init__(eventid)
        self.contents = Contents('JSON Coverages', 'coverages', eventid)

    def execute(self):
        """Create high, medium, and low resolution coverage of the mapped
        parameters.

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

        # get all of the grid layers and the geodict
        if container.getDataType() != 'grid':
            raise NotImplementedError('coverage module can only function on '
                                      'gridded data, not sets of points')

        imtlist = container.getIMTs()
        for imtype in imtlist:
            component, imtype = imtype.split('/')
            fileimt = oq_to_file(imtype)
            oqimt = imt.from_string(imtype)

            imtdict = container.getIMTGrids(imtype, component)
            grid_data = imtdict['mean']
            metadata = imtdict['mean_metadata']

            if imtype == 'MMI':
                description = 'Modified Mercalli Intensity',
                property_id = "https://earthquake.usgs.gov/learn/topics/mercalli.php",  # noqa
                decimals = 1
            elif imtype == 'PGA':
                description = 'Peak Ground Acceleration',
                units = 'natural logarithm of "g"'
                symbol = 'ln(g)'
                decimals = 2
            elif imtype == 'PGV':
                description = 'Peak Ground Velocity',
                units = 'natural logarithm of centimeters per second'
                symbol = 'ln(cm/s)'
                decimals = 2
            elif imtype.startswith('SA'):
                description = str(oqimt.period) + \
                    '-second Spectral Acceleration',
                units = 'natural logarithm of "g"'
                symbol = 'ln(g)'
                decimals = 2
            else:
                raise TypeError("Unknown IMT in coverage module")

            for i in range(3):
                if i == 0:
                    resolution = 'high'
                    fgrid = grid_data
                    decimation = 1
                elif i == 1:
                    resolution = 'medium'
                    fgrid = gaussian_filter(grid_data, sigma=1)
                    decimation = 2
                elif i == 2:
                    resolution = 'low'
                    fgrid = gaussian_filter(grid_data, sigma=2)
                    decimation = 4

                rgrid = fgrid[::decimation, ::decimation]
                ny, nx = rgrid.shape
                rnd_grd = np.flipud(np.around(
                    rgrid, decimals=decimals)).flatten()
                if imtype == 'MMI':
                    rnd_grd = np.clip(rnd_grd, 1.0, 10.0)
                xstart = metadata["xmin"]
                xstop = metadata["xmin"] + \
                    (nx - 1) * decimation * metadata["dx"]
                ystart = metadata["ymin"]
                ystop = metadata["ymin"] + \
                    (ny - 1) * decimation * metadata["dy"]

                coverage = {"type": "Coverage",
                            "domain": {
                                "type": "Domain",
                                "domainType": "Grid",
                                "axes": {
                                    "x": {"start": xstart,
                                          "stop": xstop,
                                          "num": nx},
                                    "y": {"start": ystart,
                                          "stop": ystop,
                                          "num": ny}
                                },
                                "referencing": [{
                                    "coordinates": ["x", "y"],
                                    "system": {
                                        "type": "GeographicCRS",
                                        "id": "http://www.opengis.net/def/crs/OGC/1.3/CRS84"  # noqa
                                    }
                                }]
                            },
                            "parameters": {
                                imtype: {
                                    "type": "Parameter",
                                    "description": {
                                        "en": description
                                    },
                                    "observedProperty": {
                                        "id": property_id,
                                        "label": {
                                            "en": imtype
                                        }
                                    },
                                }
                            },
                            "ranges": {
                                imtype: {
                                    "type": "NdArray",
                                    "dataType": "float",
                                    "axisNames": ["y", "x"],
                                    "shape": [ny, nx],
                                    "values": rnd_grd.tolist()
                                }
                            }
                }
                if imtype == 'MMI':
                    coverage["parameters"]["MMI"]["preferredPalette"] = {
                                        "colors": [
                                            "rgb(255, 255, 255)",
                                            "rgb(255, 255, 255)",
                                            "rgb(191, 204, 255)",
                                            "rgb(160, 230, 255)",
                                            "rgb(128, 255, 255)",
                                            "rgb(122, 255, 147)",
                                            "rgb(255, 255, 0)",
                                            "rgb(255, 200, 0)",
                                            "rgb(255, 145, 0)",
                                            "rgb(255, 0, 0)",
                                            "rgb(200, 0, 0)"
                                        ],
                                        "extent": [0, 10],
                                        "interpolation": "linear"
                                    }
                else:
                    coverage["parameters"][imtype]["unit"] = {
                        'label': {
                            "en": units
                        },
                        "symbol": {
                            'value': symbol,
                            'type': "http://www.opengis.net/def/uom/UCUM/"
                        }
                    }

                if component == 'GREATER_OF_TWO_HORIZONTAL':
                    fname = 'coverage_%s_%s_res.covjson' % (fileimt,
                                                            resolution)
                else:
                    fname = 'coverage_%s_%s_%s_res.covjson' % (fileimt,
                                                               resolution,
                                                               component)
                filepath = os.path.join(datadir, fname)
                with open(filepath, 'w') as outfile:
                    json.dump(coverage, outfile, separators=(',', ':'))
                self.contents.addFile(
                    imtype + "_" + resolution + '_res_coverage',
                    resolution + '-res ' + imtype.upper() + ' Coverage',
                    'Coverage of ' + resolution + ' resolution ' + imtype,
                    fname, 'application/json')
        container.close()
