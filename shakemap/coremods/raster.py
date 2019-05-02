# stdlib imports
import os.path
import zipfile

# third party imports
from impactutils.io.smcontainers import ShakeMapOutputContainer
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths
from shakelib.utils.imt_string import oq_to_file

FORMATS = {
    'shapefile': ('ESRI Shapefile', 'shp'),
    'geojson': ('GeoJSON', 'json')
}

DEFAULT_FILTER_SIZE = 10


class RasterModule(CoreModule):
    """
    raster -- Generate GIS raster files of all IMT values from
                    shake_result.hdf.
    """

    command_name = 'raster'
    targets = [r'products/raster\.zip']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(RasterModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)

    def execute(self):
        """
        Write raster.zip file containing ESRI Raster files of all the IMTs
        in shake_result.hdf.

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
        if container.getDataType() != 'grid':
            raise NotImplementedError('raster module can only operate on '
                                      'gridded data, not sets of points')

        # create GIS-readable .flt files of imt and uncertainty
        self.logger.debug('Creating GIS grids...')
        layers = container.getIMTs()

        # Package up all of these files into one zip file.
        zfilename = os.path.join(datadir, 'raster.zip')
        zfile = zipfile.ZipFile(zfilename, mode='w',
                                compression=zipfile.ZIP_DEFLATED)

        files_written = []
        for layer in layers:
            _, layer = layer.split('/')
            fileimt = oq_to_file(layer)
            # This is a bit hacky -- we only produce the raster for the
            # first IMC returned. It should work as long as we only have
            # one IMC produced per ShakeMap run.
            imclist = container.getComponents(layer)
            imtdict = container.getIMTGrids(layer, imclist[0])
            mean_grid = Grid2D(imtdict['mean'],
                               GeoDict(imtdict['mean_metadata']))
            std_grid = Grid2D(imtdict['std'],
                              GeoDict(imtdict['std_metadata']))
            mean_gdal = GDALGrid.copyFromGrid(mean_grid)
            std_gdal = GDALGrid.copyFromGrid(std_grid)
            mean_fname = os.path.join(datadir, '%s_mean.flt' % fileimt)
            mean_hdr = os.path.join(datadir, '%s_mean.hdr' % fileimt)
            std_fname = os.path.join(datadir, '%s_std.flt' % fileimt)
            std_hdr = os.path.join(datadir, '%s_std.hdr' % fileimt)
            self.logger.debug('Saving %s...' % mean_fname)
            mean_gdal.save(mean_fname)
            files_written.append(mean_fname)
            files_written.append(mean_hdr)
            self.logger.debug('Saving %s...' % std_fname)
            std_gdal.save(std_fname)
            files_written.append(std_fname)
            files_written.append(std_hdr)
            zfile.write(mean_fname, '%s_mean.flt' % fileimt)
            zfile.write(mean_hdr, '%s_mean.hdr' % fileimt)
            zfile.write(std_fname, '%s_std.flt' % fileimt)
            zfile.write(std_hdr, '%s_std.hdr' % fileimt)

        zfile.close()

        # nuke all of the copies of the files we just put in the zipfile
        for file_written in files_written:
            os.remove(file_written)

        self.contents.addFile('rasterData', 'ESRI Raster Files',
                              'Data and uncertainty grids in ESRI raster '
                              'format',
                              'raster.zip', 'application/zip')
        container.close()
