import matplotlib.pyplot as plt
from osgeo import osr
from osgeo import ogr
import numpy as np
import pandas as pd
import os
# import ast
from shutil import copyfile
from osgeo import gdal
gdal.UseExceptions()


class SourceProcessing(object):
    '''
    Geospatial functions for use with rectangular grid finite-difference models. The purpose
    of this class is to populate a raster grid (in GeoTiff format) with new information, 
    although it may also be useful for reading, extracting, and displaying grids from existing 
    GeoTiffs. To properly initialize this class, you must use either read_raster or create_raster
    to create the initial raster grid ("old_array"). The "new_array" is optional, and it can be
    created by scaling and sampling an existing raster, by rasterizing an existing vector data
    layer, or by simply assigning an existing numpy grid (with the same dimensions as the 
    starting raster) to either old_array or new_array.

    Parameters
    ----------
    nodata: (float, int, or np.nan)
        Set the value to be as no data for both input and output

    Attributes
    ----------
    old_array : array (nrow, ncol)
        The numpy array that is either read as part of the input raster (read_raster), or, 
        if a new raster is created rather being read in, an array of zeros (create_raster).
    new_array : array (nrow, ncol)
        The array that is produced by either rasterizing vector data or resampling and 
        reprojecting raster data.  new_array has the same spatial attributes as old_array. 
    ncol, nrow : ints
        number of rows and columns. These may coincide with a MODFLOW grid.
    gt : list of floats
        6-element geotransform list [C, A, B, F, E, D]. 
    output_raster_prj: WKT format
        projected coordinate system for output in Well-Known Text format (e.g. shapefile *.prj file)
    nodata : int or float or np.nan
        value to use as missing data in the output raster
        
    Methods
    -------
    read_raster(src_pth)
        Reads an existing raster.
    create_raster(theta, origin, LX, LY, nrow, ncol, output_raster_proj)
        Creates a blank model grid in memory.
    prj_coords_to_array_coords(x, y)
        Transforms coordinates.
    array_coords_to_prj_coords(r, c)
        Transforms coordinates.
    process_raster_data(src, method, conversion=1.0)
        Interpolates raster data on to model grid.
    process_vector_data(src, method, conversion=1.0)
        Interpolate vector data on to model grid.
    write_raster(dst_file)
        Writes an array to a GeoTiff file.
    plot_raster(which_raster='old')
        Plots a raster image in real-world coordinates
    make_clockwise(coords)
        Detects a counterclockwise polygon and reverses it.
    dbf2df(dbf_path, index=None, cols=False, incl_index=False)
        Reads a dbf file into a Pandas DataFrame
        
    Notes
    -----
    Geotransform list gives the coordinates of each pixel from the upper left pixel.
    If there is no rotation, B=D=0. If cells are square, A=-E.   
    Letter designations come from the original GDAL documentation.

    C = x coordinate in map units of the upper left corner of the upper left pixel
    A = distance from C along x axis to upper right pixel corner of the upper left pixel
    B = distance from C along x axis to lower left pixel corner of the upper left pixel,
    F = y coordinate in map units of the upper left corner of the upper left pixel
    E = distance from C along y axis to lower left pixel corner of the upper left pixel
    D = distance from C along y axis to upper right pixel corner of the upper left pixel
    
    Returns
    -------
    SourceProcessing instance
    '''

    def __init__(self, nodata=-9999):
        self.nodata = nodata

    def read_raster(self, src_pth):
        '''
        Reads an existing raster that can be used 
        1. to simply plot the image in real world coordinates,
        2. as a model-grid template for populating the grid with new information, and 
        3. to supply the information needed to transform both ways between a set of 
           model (row, col) coordinates and real-world coordinates.

        Parameters
        ----------
        src_pth: str
            Relative or absolute path, including file name, to a valid data source.
            Valid sources include shapefiles, GeoTiff, and ESRI Grid formats, and many more'''
        assert os.path.exists(src_pth), 'raster source does not exist'
        self.src = gdal.Open(src_pth)
        band = self.src.GetRasterBand(1)

        self.nrow = int(self.src.RasterYSize)
        self.ncol = int(self.src.RasterXSize)
        self.gt = self.src.GetGeoTransform()
        self.output_raster_prj = self.src.GetProjection()
        self.old_array = band.ReadAsArray()

        self.src = None
        
        self._get_coordinates()
        self._make_transforms()

    def create_raster(self, theta, origin, LX, LY, nrow, ncol, output_raster_proj):
        '''
        Creates a new blank raster that can be used 1. as a model-grid template for populating 
        the grid with new information, and 2. to supply the information needed to transform 
        between a set of model (row, col) coordinates and real-world coordinates.

        Parameters
        ----------
        theta: float
            counterclockwise rotation from positive x axis in radians
        origin: tuple of floats
            projected coordinates of the upper left corner
        LX, LY: floats
            grid cell dimensions (pixel size) in x and y 
        nrow, ncol: floats
            number of cells (pixels in each dirction)
        output_raster_prj: WKT format
            projected coordinate system for output in Well-Known Text format (e.g. shapefile *.prj file)'''
        self.theta = theta
        self.origin = origin
        A = LX * np.cos(theta)
        B = LY * np.sin(theta)
        D = LX * np.sin(theta)
        E = LY * -np.cos(theta)

        self.nrow = int(nrow)
        self.ncol = int(ncol)
        self.gt = (origin[0], A, B, origin[1], D, E)
        self.output_raster_prj = output_raster_proj
        self.old_array = np.zeros((self.nrow, self.ncol))
        
        self._get_coordinates()
        self._make_transforms()

    def _make_transforms(self):
        # make sure the gt list has been created
        assert isinstance(
            self.gt, tuple), 'Make sure either read_raster or create_model_grid has been run first'
        # format the geotransformation list into an affine transformation matrix))
        forward_transform = np.array(self.gt).reshape(2, -1)
        # add a row to get homogeneous coodinates (offsets are in the first column)
        self.forward_transform = np.vstack((forward_transform, [1, 0, 0]))
        # invert the forward transform
        self.reverse_transform = np.linalg.inv(self.forward_transform)

    def prj_coords_to_array_coords(self, x, y):
        '''
        Transform real-world coordinates to model grid coordinates (row, column). 

        Parameters
        ----------
        x, y: arrays of floats, one each for x and y (n, 1)
            real-world coordinates in whatever system the model grid is in

        Returns
        -------
        array : (n, 2)
            row and col coordinates. Integer portion is the zero-based row or column;
            decimal portion is the relative coordinate within a cell (0 to 1). 
            Row coordinate increases from top (row 0) to bottom (last row). 
        '''

        self._make_transforms()
        assert x.shape[0] == y.shape[0], 'x and y have to have the same dimensions'
        ones = np.ones(x.shape[0])
        wpts = np.column_stack((x, y, ones))
        wpp = self.reverse_transform.dot(wpts.T)
        return wpp[1:, ].T

    def array_coords_to_prj_coords(self, r, c):
        '''
        Transform model grid (row, col) coordinates to real-world coordinates

        Parameters
        ----------
        r, c: arrays of floats, one each for row and col (n, 1)
            real-world coordinates in whatever system the model grid is in

        Returns
        -------
        array : (n, 2) 
            x and y real-world coordinates
        '''
        self._make_transforms()
        assert r.shape[0] == c.shape[0], 'r and c have to have the same dimensions'
        ones = np.ones(r.shape[0])
        wpts = np.column_stack((ones, c, r))
        dat = self.forward_transform.dot(wpts.T).T
        return dat[:, :-1]

    def process_raster_data(self, src, method, conversion=1.0):
        '''
        Takes a raster data source (ESRI grid, GeoTiff, .IMG and many other formats)
        and returns a numpy array. Arrangement of pixels is given as input and may 
        correspond to a MODFLOW grid.

        Parameters
        ----------
        src : str
            complete path to raster data source
        method : str
            gdal method for interpolation. See notes.
        conversion : float
            factor to be applied to raw data values, for eaxmple to change units

        Returns
        -------
        array : (nrow, ncol)
            Raster data source projected onto model grid. Returns a zero array with the 
            correct shape if the source does not exist.
            
        Notes
        -----
        Choices for method are (not all may be available depending on version of gdal)      
            gdal.GRA_NearestNeighbour 
                Nearest neighbour (select on one input pixel)
            gdal.GRA_Bilinear
                Bilinear (2x2 kernel)
            gdal.GRA_Cubic
                Cubic Convolution Approximation (4x4 kernel)
            gdal.GRA_CubicSpline
                Cubic B-Spline Approximation (4x4 kernel)
            gdal.GRA_Lanczos
                Lanczos windowed sinc interpolation (6x6 kernel)
            gdal.GRA_Average
                Average (computes the average of all non-NODATA contributing pixels)
            gdal.GRA_Mode
                Mode (selects the value which appears most often of all the sampled points)
            gdal.GRA_Max
                Max (selects maximum of all non-NODATA contributing pixels)
            gdal.GRA_Min
                Min (selects minimum of all non-NODATA contributing pixels)
            gdal.GRA_Med
                Med (selects median of all non-NODATA contributing pixels)
            gdal.GRA_Q1
                Q1 (selects first quartile of all non-NODATA contributing pixels)
            gdal.GRA_Q3
                Q3 (selects third quartile of all non-NODATA contributing pixels)

        '''
        if os.path.exists(src):
            rast = gdal.Open(src)
            dest = self._make_grid()

            gdal.ReprojectImage(rast, dest, rast.GetProjection(),
                                self.output_raster_prj, method)

            grid = dest.GetRasterBand(1).ReadAsArray()
            grid = grid * conversion

            dest = None
            rast = None

        else:
            grid = np.ones((self.nrow, self.ncol)) * self.nodata
            print(
                'Data not processed for\n{}\n Check that the file exists and path is correct'.format(src))

        self.new_array = grid

    def process_vector_data(self, src, attribute, layer=0, all_touched=False):
        '''
        Takes a vector data source (e.g. ESRI shapefile) and returns a numpy array.
        Arrangement of pixels is given as input and may correspond to a MODFLOW grid.

        Parameters
        ----------
        src : str
            complete path to vector data source
        attribute : str
            field in data table to assign to rasterized pixels

        Returns
        -------
        array : (nrow, ncol)
            Vector data source projected onto model grid. Returns a zero array 
            with the correct shape if the source does not exist.
        '''
        if os.path.exists(src):
            datasource = ogr.Open(src)
            layer = datasource.GetLayer(iLayer=layer)

            dest = self._make_grid()
            args = 'ATTRIBUTE={}'.format(attribute)
            args2 = 'ALL_TOUCHED={}'.format(all_touched)
            gdal.RasterizeLayer(dest, [1], layer, options=[args, args2])

            grid = dest.GetRasterBand(1).ReadAsArray()

            src = None
            dst = None

        else:
            grid = np.ones((self.nrow, self.ncol)) * self.nodata
            print(
                'Data not processed for\n{}\n Check that the file exists and path is correct'.format(src))

        self.new_array = grid

    def write_raster(self, dst_file, which_raster='old'):
        '''
        Writes a numpy array to a GeoTiff file.

        Parameters
        ----------
        dst_file : str
            path name of file to write
        '''
        driver = gdal.GetDriverByName("GTiff")
        dst = driver.Create(dst_file, self.ncol,
                            self.nrow, 1, gdal.GDT_Float32)
        dst.SetGeoTransform(self.gt)
        dst.SetProjection(self.output_raster_prj)
        band = dst.GetRasterBand(1)
        band.SetNoDataValue(self.nodata)
        if which_raster == 'old':
            band.WriteArray(self.old_array)
        elif which_raster == 'new':
            band.WriteArray(self.new_array)
        else:
            raise Exception(
                "which_raster has to be 'old' or 'new'. \n You entered {}".format(which_raster))
        dst = None

    def _make_grid(self):
        '''
        Creates a blank raster image in memory.
        '''
        mem_drv = gdal.GetDriverByName('MEM')
        grid_ras = mem_drv.Create(
            '', self.ncol, self.nrow, 1, gdal.GDT_Float32)
        grid_ras.SetGeoTransform(self.gt)
        grid_ras.SetProjection(self.output_raster_prj)
        band = grid_ras.GetRasterBand(1)
        band.SetNoDataValue(self.nodata)
        self.new_array = np.zeros((self.nrow, self.ncol))
        band.WriteArray(self.new_array)
        return grid_ras
        
    def _get_coordinates(self):
        '''
        Computes the x and y coordinates for pixels at the edges of
        cells (for pcolormesh and other plotting routines) and at the
        center of cells (for contourf). The row and column attributes
        (r and c) are for cell centers where the upper left cell is 
        at (0.5, 0.5) and the bottom right cell is at (nrow - 0.5, ncol - 0.5)
        '''
        self.r, self.c = np.indices((self.nrow, self.ncol)) + 0.5
        r, c = self.r.ravel(), self.c.ravel()
        x, y = self.array_coords_to_prj_coords(r, c).T

        self.x_center = x.reshape(self.nrow, self.ncol)
        self.y_center = y.reshape(self.nrow, self.ncol)
        
        re, ce = np.indices((self.nrow + 1, self.ncol + 1))
        re, ce = re.ravel(), ce.ravel()
        x, y = self.array_coords_to_prj_coords(re, ce).T

        self.x_edge = x.reshape(self.nrow + 1, self.ncol + 1)
        self.y_edge = y.reshape(self.nrow + 1, self.ncol + 1)

    def plot_raster(self, which_raster='old', plot=True, sk={'figsize':(5,5)}, pk={'cmap':plt.cm.coolwarm}, ck={'shrink':0.5}):
        '''
        Parameters
        ----------
        which_raster: str
           valid values are 'old' to plot a raster that was read in
           with read_raster; 'new' will plot the newly created raster
        sk : dictionary
           keywords to add to subplots
        pk : dictionary
           keywords to add to pcolormesh
        ck : dictionary
           keywords to add to colorbar

        Returns
        -------
        matplotlib figure instance
            fig, ax
        '''
        if plot:
            fig, ax = plt.subplots(1, 1, **sk)
            if which_raster == 'old':
                im = ax.pcolormesh(self.x_edge, self.y_edge, self.old_array, **pk)
                fig.colorbar(im, **ck)
            elif which_raster == 'new':
                im = ax.pcolormesh(self.x_edge, self.y_edge, self.new_array, **pk)
                fig.colorbar(im, **ck)
            else:
                raise Exception(
                    "which_raster has to be 'old' or 'new'. \n You entered {}".format(which_raster))
            ax.set_aspect(1)
            return fig, ax

    def make_clockwise(self, coords):
        '''
        Function to determine direction of vertices of a polygon (clockwise or CCW).
        Probably not needed, but here just in case. 

        Parameters
        ----------
        coords : array (n, 2)
            n is number of vertices in the polygon. The last vertex is the same 
            as the first to close the polygon. The first column is x and the second is y.
        '''
        # if the points are counterclockwise, reverse them
        x1 = coords[:-1, 0]
        x2 = coords[1:, 0]
        y1 = coords[:-1, 1]
        y2 = coords[1:, 1]
        ccw = np.sum((x2 - x1) * (y2 + y1)) < 0
        if ccw:
            coords = np.flipud(coords)
            print('yup, coordinates are ccw')
            print("let's change them to CW")
        return coords

    # test data for make_clockwise

    # print('clockwise')
    # x = np.array([1, 1, 2, 2, 1])
    # y = np.array([1, 2, 2, 1, 1])
    # coords = np.array(zip(x, y))
    # c = make_clockwise(coords)
    # print( c)
    # print('\n')
    # print('CCW')
    # x = np.array([1, 2, 2, 1, 1])
    # y = np.array([1, 1, 2, 2, 1])
    # coords = np.array(zip(x, y))
    # c = make_clockwise(coords)
    # print( c)

    def dbf2df(self, dbf_path, pandas=True):
        '''
        Read a dbf file as a geopandas.GeoDataFrame, optionally converting to Pandas.

        ...

        Parameters
        ----------
        dbf_path : str
            Path to the DBF file to be read
        pandas : bool
            Flag to convert geodataframe to datatframe
 
        Returns
        -------
        df : DataFrame or GeoDataFrame
        '''
        db = gpd.read_file(dbf_path)
        if pandas:
            db = pd.DataFrame(db)

def download_and_extract(url, id=None, destination='.'):
    # importing the requests module
    import requests, zipfile, py7zr
    from io import BytesIO

    print('Downloading started')
    # Downloading the file by sending the request to the URL
    req = requests.get(url)
    
    if id is None:
        # extracting the zip file contents
        if url.endswith('zip'):
            zippy = zipfile.ZipFile(BytesIO(req.content))
            zippy.extractall(path=destination)
            
        if url.endswith('7z'):
            zippy = py7zr.SevenZipFile(BytesIO(req.content))
            zippy.extractall(path=destination)
            
    else:
        # extracting the zip file contents where its name contains 'id'
        
        if url.endswith('zip'):
            zippy = zipfile.ZipFile(BytesIO(req.content))
            file_list = [item for item in zippy.namelist() if id in item]
            for file in file_list:
                zippy.extract(file, path=destination)
                
        if url.endswith('7z'):
           zippy = py7zr.SevenZipFile(BytesIO(req.content))
           file_list = [item for item in zippy.namelist() if id in item]
           for file in file_list:
                zippy.extract(file, path=destination)        
        return file_list
    print('Downloading Completed')

