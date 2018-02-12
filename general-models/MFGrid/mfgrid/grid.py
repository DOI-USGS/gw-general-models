import numpy as np

class ModflowGrid(object):
    """
    ModflowGrid is a MODFLOW grid class.  This is a structured grid of shape 
    (nlay, nrow, ncol).  The structure assumes that the first layer is the top 
    layer, and that the first row corresponds to the north side of the grid, 
    so that index [1, 1, 1] corresponds to the top layer, and the northwest 
    corner of the grid.
    
    """

    def __init__(self, nlay, nrow, ncol, delr, delc, botm, **kwargs):
        """
        This method is used to instantiate a ModflowGrid object.  The 
        following arguments are required:
            nlay: number of layers
            nrow: number of rows
            ncol: number of columns
            delr: a ModflowArray object of size [ncol]
            delc: a ModflowArray object of size [nrow]
            botm: a ModflowArray object of size [nlay + 1, nrow, ncol]
            
        Optional arguments can be specified using **kwargs.  These are:
            xoffset: the offset of the grid in the x direction
            yoffset: the offset of the grid in the y direction
            rotation: the grid rotation angle in degrees relative to the lower  
                      left corner.
                      
        The local coordinate system has a (0, 0) x,y location that corresponds 
        to the lower left corner of the grid.
        """
        
        from gridutil import ModflowArray

        #assign the arguments to the object.
        self.nlay = nlay
        self.nrow = nrow
        self.ncol = ncol
        self.nodes = nrow * ncol
        self.delr = ModflowArray(delr, arrayshape=(ncol,)).as_ndarray()
        self.delc = ModflowArray(delc, arrayshape=(nrow,)).as_ndarray()

        #set defaults
        self.xoffset = 0.
        self.yoffset = 0.
        self.rotation = 0.
        if 'xoffset' in kwargs: self.xoffset = kwargs['xoffset']
        if 'yoffset' in kwargs: self.yoffset = kwargs['yoffset']
        if 'rotation' in kwargs: self.rotation = kwargs['rotation']

        #assign the arguments to the object.
        self.nodes = nlay * nrow * ncol
        self.botm = ModflowArray(botm, arrayshape=(nlay + 1, nrow, ncol), 
                                dtype=float).as_ndarray()

        #assign centroid and edge coordinate arrays
        self.X = self.get_local_x_array()
        self.Y = self.get_local_y_array()
        self.Xe = self.get_local_Xe_array()
        self.Ye = self.get_local_Ye_array()

        return

    def get_local_x_array(self):
        """
        Return a numpy one-dimensional float array that has the cell center x
        coordinate for every cell in the grid.

        """
        x = np.add.accumulate(self.delr) - 0.5 * self.delr
        x += self.xoffset
        return x

    def get_local_y_array(self):
        """
        Return a numpy one-dimensional float array that has the cell center x
        coordinate for every cell in the grid.

        """
        Ly = np.add.reduce(self.delc)
        y = Ly - ( np.add.accumulate(self.delc) - 0.5 * self.delc )
        y += self.yoffset
        return y

    def get_local_Xe_array(self):
        """
        Return a numpy one-dimensional float array that has the cell edge x
        coordinates for every cell in the grid.  Array is of size (ncol+1)

        """
        Xe = np.concatenate(
             ([0.], np.add.accumulate(self.delr) )
             )
        Xe += self.xoffset
        return Xe

    def get_local_Ye_array(self):
        Ly = np.add.reduce(self.delc)
        Ye = np.concatenate(
            ([Ly], Ly - np.add.accumulate(self.delc) )
             )
        Ye += self.yoffset
        return Ye

    def get_extent(self):
        '''
        Return a tuple containing two points (the lower left and the upper
        right) in global grid coordinates.
        '''
        xmin = self.Xe[0]
        xmax = self.Xe[-1]
        ymax = self.Ye[0]
        ymin = self.Ye[-1]
        return ((xmin, ymin), (xmax, ymax))

    def get_nodeid(self, k, i, j):
        '''
        Return the nodeid for the modflow grid.
        '''
        return k * self.nrow * self.ncol + i * self.ncol + j

    def get_indices(self, nodeid):
        """
        Return the indices for the specified nodeid.
        """
        if nodeid > self.nlay * self.nrow * self.ncol: 
            raise Exception('Nodeid is not valid...')
        k = int( nodeid / self.nrow / self.ncol )
        i = int( (nodeid - k * self.nrow * self.ncol) / self.ncol )
        j = nodeid - k * self.nrow * self.ncol - i * self.ncol
        return (k, i, j)

    def get_vertices(self, nodeid):
        """
        Return the four vertices for the specified nodeid.  This currently
        returns two-dimensional vertices
        3-------2   7-------6
        |  top  |   |  bot  |
        |       |   |       |
        0-------1   4-------5
        """
        (k, i, j) = self.get_indices(nodeid)
        x0 = self.Xe[j]
        y0 = self.Ye[i+1]
        x1 = self.Xe[j + 1]
        y1 = self.Ye[i + 1]
        x2 = self.Xe[j + 1]
        y2 = self.Ye[i]
        x3 = self.Xe[j]
        y3 = self.Ye[i]
        return ((x0, y0), (x1, y1), (x2, y2), (x3, y3), (x0, y0))

    def get_3dvertices(self, nodeid):
        """
        Return the eight vertices for the specified nodeid.  This currently
        returns two-dimensional vertices
        3-------2   7-------6
        |  top  |   |  bot  |
        |       |   |       |
        0-------1   4-------5
        """
        (k, i, j) = self.get_indices(nodeid)
        p0 = (self.Xe[j], self.Ye[i + 1], self.botm[k, i, j])
        p1 = (self.Xe[j + 1], self.Ye[i + 1], self.botm[k, i, j])
        p2 = (self.Xe[j + 1], self.Ye[i], self.botm[k, i, j])
        p3 = (self.Xe[j], self.Ye[i], self.botm[k, i, j])

        p4 = (self.Xe[j], self.Ye[i + 1], self.botm[k, i, j])
        p5 = (self.Xe[j + 1], self.Ye[i + 1], self.botm[k, i, j])
        p6 = (self.Xe[j + 1], self.Ye[i], self.botm[k, i, j])
        p7 = (self.Xe[j], self.Ye[i], self.botm[k, i, j])

        return (p0, p1, p2, p3, p4, p5, p6, p7)

    def draw(self, **kwargs):
        import gridplot
        gridplot.draw(self, **kwargs)
        return

    def label(self, arr=None, **kwargs):
        """
        Label the grid.  If arr=None, then (i,j) will be labeled for each cell

        """
        import gridplot
        gridplot.label(self, arr, **kwargs)
        return

    def colorflood(self, arr, **kwargs):
        """
        Color flood the array

        """
        import gridplot
        gridplot.colorflood(self, arr, **kwargs)
        return

    def intersection(self, feature, featuretype, **kwargs):
        """
        Method for intersecting the model grid with a point, line, rectangle,
        or polygon.  A rectangle intersection is faster than a polygon
        intersection.

        """
        #determine type of feature
        import gridintersect as intersect
        featuretype = featuretype.lower()
        if featuretype == 'point':
            return intersect.PointGridIntersection(self, feature, **kwargs).nodelist
        elif featuretype == 'line':
            lg = intersect.LineGridIntersection(self, feature, **kwargs)
            nodes = lg.nodelist
            lengths = lg.lengths
            return nodes, lengths
        elif featuretype == 'rectangle':
            return intersect.PointGridIntersection(self, feature, **kwargs).nodelist
        elif featuretype == 'polygon':
            # make a quick check that the polygon is closed and, if not, add the
            # first point to the end of the list
            if not isinstance(feature, np.ndarray):
                feature = np.array(feature)
            if (feature[0] != feature[-1]).any():
                feature = np.vstack((feature, feature[0]))
            pg = intersect.PolygonGridIntersection(self, feature, **kwargs)
            nodes = pg.nodelist
            areas = pg.areas
            return nodes, areas
        else:
            msg = ('Unknown feature type: ' + featuretype + '\n' +
            'feature must be point, line or polygon.')
            grid.logger.error(msg)
            raise Exception()
        return


if __name__ == '__main__':

    import grid
    import matplotlib.pyplot as plt

    #set grid defaults
    nrow = 50
    ncol = 25
    nlay = 2
    delr = 15.
    delc = 15.
    delv = 1.
    botm = np.empty((nlay + 1), dtype=float)
    for k in range(1, nlay):
        botm[k] = botm[k - 1] - delv
    mfg = grid.ModflowGrid(nlay, nrow, ncol, delr, delc, botm)

    #set arrays containing cell coordinate information
    X = mfg.X
    Y = mfg.Y
    Xe = mfg.Xe
    Ye = mfg.Ye

    pltnum = 1
    #create a figure showing points at cell centers
    try:
        plt.close('all')
    except:
        pass
    plt.subplot(1, 1, 1, aspect='equal')
    mfg.draw()
    fname = 'mfg' + str(pltnum) + '.png'
    print ('saving figure: ', fname)
    plt.savefig('mfg' + str(pltnum) + '.png')
    pltnum += 1

    #create a figure showing grid and nodenumbers
    try:
        plt.close('all')
    except:
        pass
    plt.subplot(1, 1, 1, aspect='equal')
    mfg.draw()
    mfg.label()
    plt.xlim(Xe[0], Xe[-1])
    plt.ylim(Ye[-1], Ye[0])
    fname = 'mfg' + str(pltnum) + '.png'
    print ('saving figure: ', fname)
    plt.savefig(fname)
    pltnum += 1

    #test the grid-point intersection
    feature = (150., 300.)
    featuretype = 'point'
    i0, j0 = mfg.intersection(feature, featuretype)
    intersection_array = np.zeros( (nrow, ncol), dtype=np.int)
    intersection_array[i0, j0] = 1
    print ('point intersection result: ', i0, j0)
        #create a figure showing grid and nodenumbers
    try:
        plt.close('all')
    except:
        pass
    plt.subplot(1, 1, 1, aspect='equal')
    imshowextent = (Xe[0], Xe[-1], Ye[-1], Ye[0])
    plt.imshow(intersection_array, interpolation='nearest', extent=imshowextent)
    plt.plot(feature[0], feature[1], 'bo')
    plt.xlim(Xe[0], Xe[-1])
    plt.ylim(Ye[-1], Ye[0])
    fname = 'mfg' + str(pltnum) + '.png'
    print ('saving figure: ', fname)
    plt.savefig(fname)
    pltnum += 1

    #test the grid-line intersection
    feature = [(10, 10), (50, 50), (500, 300), (100, 700)]
    featuretype = 'line'
    nodelist, lengths = mfg.intersection(feature, featuretype)
    print ('line intersection result: ', nodelist, lengths)
    intersection_array = np.zeros( (nrow, ncol), dtype=np.float)
    for (i0, j0), length in zip(nodelist, lengths):
        intersection_array[i0, j0] = length
    try:
        plt.close('all')
    except:
        pass
    plt.subplot(1, 1, 1, aspect='equal')
    plt.imshow(intersection_array, interpolation='nearest', extent=imshowextent)
    line = np.array(feature)
    plt.plot(line[:, 0], line[:, 1], 'k-')
    plt.xlim(Xe[0], Xe[-1])
    plt.ylim(Ye[-1], Ye[0])
    fname = 'mfg' + str(pltnum) + '.png'
    print ('saving figure: ', fname)
    plt.savefig(fname)
    pltnum += 1

    #test the grid-polygon intersection
    feature = [(10, 10), (50, 50), (500, 300), (100, 700), (10, 10)]
    featuretype = 'polygon'
    nodelist, areas = mfg.intersection(feature, featuretype)
    print ('polygon intersection result: ', nodelist, areas)
    intersection_array = np.zeros( (nrow, ncol), dtype=np.float)
    for (i0, j0), area in zip(nodelist, areas):
        intersection_array[i0, j0] = area
    try:
        plt.close('all')
    except:
        pass
    plt.subplot(1, 1, 1, aspect='equal')
    plt.imshow(intersection_array, interpolation='nearest', extent=imshowextent)
    line = np.array(feature)
    plt.plot(line[:, 0], line[:, 1], 'k-')
    plt.xlim(Xe[0], Xe[-1])
    plt.ylim(Ye[-1], Ye[0])
    fname = 'mfg' + str(pltnum) + '.png'
    print ('saving figure: ', fname)
    plt.savefig(fname)
    pltnum += 1


    
