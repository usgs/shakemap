# stdlib imports
import os.path
import logging
from datetime import datetime
import re
from collections import OrderedDict

# third party imports
import numpy as np
from impactutils.io.smcontainers import ShakeMapOutputContainer
from mapio.shake import ShakeGrid
from mapio.geodict import GeoDict

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths
import shakemap
from shakelib.rupture import constants


def _oq_to_gridxml(oqimt):
    """
    Convert openquake IMT nomenclature to grid.xml friendly form.

    Note: The grid.xml form only handles periods up to 9.9, after
    that there is no way to tell the difference between 10.0 and 1.0.

    Examples:
    SA(1.0) (Spectral Acceleration at 1 second) -> PSA10
    SA(0.3) (Spectral Acceleration at 0.3 second) -> PSA03
    SA(15.0) (Spectral Acceleration at 15 seconds) -> NOT SUPPORTED
    SA(3) (Spectral Acceleration at 3 seconds) -> PSA30
    SA(.5) (Spectral Acceleration at 0.5 seconds) -> PSA05


    Args:
        oqimt (str): Openquake IMT nomenclature string.
    Returns:
        str: grid.xml friendly IMT string.
    Raises:
        ValueError: when there is no corresponding filename-friendly
            IMT representation, or when frequency exceeds 9.9.
    """
    if oqimt in ['PGA', 'PGV', 'MMI']:
        return oqimt
    float_pattern = r"[-+]?\d*\.\d+|\d+"
    periods = re.findall(float_pattern, oqimt)
    if not len(periods):
        fmt = 'IMT string "%s" has no file-name friendly representation.'
        raise ValueError(fmt % oqimt)
    period = periods[0]
    if period.find('.') < 0:
        integer = period
        fraction = '0'
    else:
        integer, fraction = period.split('.')
        if not len(integer):
            integer = '0'
    if int(integer) >= 10:
        raise ValueError('Periods >= than 10 seconds not supported.')
    fileimt = 'PSA%s%s' % (integer, fraction)
    return fileimt


class GridXMLModule(CoreModule):
    """
    gridxml -- Create grid.xml and uncertainty.xml files from shake_result.hdf.
    """

    command_name = 'gridxml'
    targets = [r'products/grid\.xml', r'products/uncertainty\.xml']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(GridXMLModule, self).__init__(eventid)
        self.contents = Contents('XML Grids', 'gridxml', eventid)

    def execute(self):
        """Create grid.xml and uncertainty.xml files.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        logger = logging.getLogger(__name__)
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
            raise NotImplementedError('gridxml module can only function on '
                                      'gridded data, not sets of points')

        components = container.getComponents()
        for component in components:
            xml_types = ['grid', 'uncertainty']
            for xml_type in xml_types:
                layers = OrderedDict()
                field_keys = OrderedDict()
                gridnames = container.getIMTs(component)
                for gridname in gridnames:
                    imt_field = _oq_to_gridxml(gridname)
                    imtdict = container.getIMTGrids(gridname, component)
                    if xml_type == 'grid':
                        grid_data = imtdict['mean']
                        metadata = imtdict['mean_metadata']
                    elif xml_type == 'uncertainty':
                        grid_data = imtdict['std']
                        metadata = imtdict['std_metadata']

                    units = metadata['units']
                    digits = metadata['digits']
                    # convert from HDF units to legacy grid.xml units
                    if xml_type == 'grid':
                        if units == 'ln(cm/s)':
                            grid_data = np.exp(grid_data)
                            units = 'cm/s'
                        elif units == 'ln(g)':
                            grid_data = np.exp(grid_data) * 100
                            units = '%g'
                        else:
                            pass

                    if xml_type == 'grid':
                        layers[imt_field] = grid_data
                        field_keys[imt_field] = (units, digits)
                    else:
                        layers['STD' + imt_field] = grid_data
                        field_keys['STD' + imt_field] = (units, digits)

                if xml_type == 'grid':
                    grid_data, _ = container.getArray([], 'vs30')
                    units = 'm/s'
                    digits = metadata['digits']
                    layers['SVEL'] = grid_data
                    field_keys['SVEL'] = (units, digits)

                geodict = GeoDict(metadata)

                config = container.getConfig()

                # event dictionary
                info = container.getMetadata()
                event_info = info['input']['event_information']
                event_dict = {}
                event_dict['event_id'] = event_info['event_id']
                event_dict['magnitude'] = float(event_info['magnitude'])
                event_dict['depth'] = float(event_info['depth'])
                event_dict['lat'] = float(event_info['latitude'])
                event_dict['lon'] = float(event_info['longitude'])
                try:
                    event_dict['event_timestamp'] = datetime.strptime(
                        event_info['origin_time'], constants.TIMEFMT)
                except ValueError:
                    event_dict['event_timestamp'] = datetime.strptime(
                        event_info['origin_time'], constants.ALT_TIMEFMT)
                event_dict['event_description'] = event_info['location']
                event_dict['event_network'] = \
                    info['input']['event_information']['eventsource']
                event_dict['intensity_observations'] =\
                    info['input']['event_information']['intensity_observations']
                event_dict['seismic_stations'] =\
                    info['input']['event_information']['seismic_stations']
                if info['input']['event_information']['fault_ref'] == 'Origin':
                    event_dict['point_source'] = 'True'
                else:
                    event_dict['point_source'] = 'False'

                # shake dictionary
                shake_dict = {}
                shake_dict['event_id'] = event_dict['event_id']
                shake_dict['shakemap_id'] = event_dict['event_id']
                shake_dict['shakemap_version'] = \
                    info['processing']['shakemap_versions']['map_version']
                shake_dict['code_version'] = shakemap.__version__
                ptime = info['processing']['shakemap_versions']['process_time']
                try:
                    shake_dict['process_timestamp'] = datetime.strptime(
                        ptime, constants.TIMEFMT)
                except ValueError:
                    shake_dict['process_timestamp'] = datetime.strptime(
                        ptime, constants.ALT_TIMEFMT)

                shake_dict['shakemap_originator'] = \
                    config['system']['source_network']
                shake_dict['map_status'] = config['system']['map_status']
                shake_dict['shakemap_event_type'] = 'ACTUAL'
                if event_dict['event_id'].endswith('_se'):
                    shake_dict['shakemap_event_type'] = 'SCENARIO'

                shake_grid = ShakeGrid(
                    layers, geodict, event_dict,
                    shake_dict, {}, field_keys=field_keys)
                if component == 'GREATER_OF_TWO_HORIZONTAL':
                    fname = os.path.join(datadir, '%s.xml' % xml_type)
                else:
                    fname = os.path.join(datadir, '%s_%s.xml' % (xml_type,
                                                                 component))
                logger.debug('Saving IMT grids to %s' % fname)
                shake_grid.save(fname)  # TODO - set grid version number
                cname = os.path.split(fname)[1]

                if xml_type == 'grid':
                    self.contents.addFile(
                        'xmlGrids',
                        'XML Grid',
                        'XML grid of %s ground motions' % component,
                        cname, 'text/xml')
                else:
                    self.contents.addFile(
                        'uncertaintyGrids',
                        'Uncertainty Grid',
                        'XML grid of %s uncertainties' % component,
                        cname, 'text/xml')

        container.close()
