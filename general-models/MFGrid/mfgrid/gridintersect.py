'''
Classes and methods for performing spatial analyses with a ModflowGrid
object or a ModflowGrid2D object.
'''

class PointGridIntersection(object):
    '''
    
    '''
    def __init__(self, grid, point, **kwargs):
        '''
        Find the point in the grid.  Return None if not found in range, 
        otherwise return (row, column) in zero-based indexing.  If x, y, fall 
        directly on a cell edge, return the smaller cell indices.
        
        Arguments:
        
            *point*: A tuple containing the x and y coordinates of the point,
                (x, y) or a tuple containing the x, y, and z coordinates of 
                the point (x, y, z).
        
        Store the node (ipos, jpos) or (kpos, ipos, jpos) in self.nodelist.
        Set self.nodelist = None if point is not in grid.        
        '''
        from gridutil import ModflowGridIndices
        
        self.nodelist = None
        
        #two dimensional point
        Xe = grid.Xe
        Ye = grid.Ye
        x = point[0]
        jpos = ModflowGridIndices.find_position_in_array(Xe, x)
        y = point[1]
        ipos = ModflowGridIndices.find_position_in_array(Ye, y)
        self.nodelist = (ipos, jpos)

        #three dimensional point
        if len(point) == 3:
            #find k
            z = point[2]
            kpos = ModflowGridIndices.find_position_in_array(
                    grid.botm[:, ipos, jpos], z)
            if kpos is None:
                self.nodelist = None
            else:
                self.nodelist = (kpos, ipos, jpos)
        return

class LineGridIntersection(object):
    '''
    This class contains the methods for intersecting a polyline with a
    MODFLOW grid object.  The nodes intersecting the line and the lengths of
    each line are stored in the nodes and lengths properties of the class
    object.
    '''

    def __init__(self, grid, line, keepzerolengths=False):
        '''
        Create the line intersection object.  Store the list of nodes that
        intersect the line in self.nodelist.  Store the corresponding lengths in
        self.lengths.
        
        Nodes are represented as (row, col).
        
        Arguments:
        
            *grid*: A ModflowGrid or ModflowGrid2D object.
            
            *line*: A tuple or list of points defining a line.
            
            *keepzerolengths*: A true or false flag indicating whether line
                segments that have zero lengths should be included in the
                list of nodes and list of lengths.  Zero length line segments
                can occur when a line touches a cell edge.
            
        '''
        from shapely.geometry import LineString, Polygon, MultiLineString, box
        self.grid = grid
        self.line = line
        (xmin, ymin), (xmax, ymax) = self.grid.get_extent()
        pl = box(xmin, ymin, xmax, ymax)
        ls = LineString(self.line)
        lineclip = ls.intersection(pl)
        self.nodelist = []
        self.lengths = []
        if lineclip.length == 0.:            
            return
        if lineclip.geom_type is 'MultiLineString': #there are multiple lines
            for l in lineclip:
                self.get_nodes_intersecting_linestring(l)
        else:
            self.get_nodes_intersecting_linestring(lineclip)

        #eliminate any nodes that have a zero length
        if not keepzerolengths:
            tempnodes = []
            templengths = []
            for i in range(len(self.nodelist)):
                if self.lengths[i] > 0:
                    tempnodes.append(self.nodelist[i])
                    templengths.append(self.lengths[i])
            self.nodelist = tempnodes
            self.lengths = templengths
        return
        
    def get_nodes_intersecting_linestring(self, linestring):
        '''
        Intersect the linestring with the model grid and return a list of 
        node indices and the length of the line in that node.
        '''
        from shapely.geometry import LineString, Polygon, MultiLineString, box

        #start at the beginning of the line
        x, y = linestring.xy
        x0 = x[0]
        y0 = y[0]
        i, j = self.grid.intersection((x0, y0), 'point')
        xmin = self.grid.Xe[j]
        xmax = self.grid.Xe[j + 1]
        ymax = self.grid.Ye[i]
        ymin = self.grid.Ye[i + 1]
        pl = box(xmin, ymin, xmax, ymax)
        length = linestring.intersection(pl).length
        self.lengths.append(length)
        self.nodelist.append( (i, j) )
        n = 0
        while True:
            (i, j) = self.nodelist[n]
            self.check_adjacent_cells_intersecting_line(linestring, i, j)
            if n == len(self.nodelist) - 1:
                break
            n += 1
        return
        
    def check_adjacent_cells_intersecting_line(self, linestring, i, j):
        from shapely.geometry import LineString, Polygon, box
        
        #check to left
        if j > 0:
            ii = i
            jj = j - 1
            if (ii, jj) not in self.nodelist:
                xmin = self.grid.Xe[jj]
                xmax = self.grid.Xe[jj + 1]
                ymax = self.grid.Ye[ii]
                ymin = self.grid.Ye[ii + 1]
                pl = box(xmin, ymin, xmax, ymax)
                if linestring.intersects(pl):
                    length = linestring.intersection(pl).length
                    self.nodelist.append( (ii,jj) )
                    self.lengths.append(length)

        #check to right
        if j < self.grid.ncol - 1:
            ii = i
            jj = j + 1
            if (ii, jj) not in self.nodelist:
                xmin = self.grid.Xe[jj]
                xmax = self.grid.Xe[jj + 1]
                ymax = self.grid.Ye[ii]
                ymin = self.grid.Ye[ii + 1]
                pl = box(xmin, ymin, xmax, ymax)
                if linestring.intersects(pl):
                    length = linestring.intersection(pl).length
                    self.nodelist.append( (ii,jj) )
                    self.lengths.append(length)
        
        #check to back
        if i > 0:
            ii = i - 1
            jj = j
            if (ii, jj) not in self.nodelist:
                xmin = self.grid.Xe[jj]
                xmax = self.grid.Xe[jj + 1]
                ymax = self.grid.Ye[ii]
                ymin = self.grid.Ye[ii + 1]
                pl = box(xmin, ymin, xmax, ymax)
                if linestring.intersects(pl):
                    length = linestring.intersection(pl).length
                    self.nodelist.append( (ii,jj) )
                    self.lengths.append(length)

        #check to front
        if i < self.grid.nrow - 1:
            ii = i + 1
            jj = j
            if (ii, jj) not in self.nodelist:
                xmin = self.grid.Xe[jj]
                xmax = self.grid.Xe[jj + 1]
                ymax = self.grid.Ye[ii]
                ymin = self.grid.Ye[ii + 1]
                pl = box(xmin, ymin, xmax, ymax)
                if linestring.intersects(pl):
                    length = linestring.intersection(pl).length
                    self.nodelist.append( (ii,jj) )
                    self.lengths.append(length)

        return


class RectangleGridIntersection(object):
    '''
    
    '''
    def __init__(self, grid, rectangle, **kwargs):
        '''
        Given a rectangle defined as [(xmin, ymin), (xmax, ymax)]
        return the cells (k, i, j) that are within or touching
        the rectangle.  This is faster than using the more generic
        PolygonGridIntersect approach.
        
        Arguments:
        
            *grid*: The ModflowGrid or ModflowGrid2D object.
            
            *rectangle*: A tuple containing ((xmin, ymin), (xmax, ymax))
        
        
        '''
        from shapely.geometry import Point, Polygon, box
        from gridutil import ModflowGridIndices
        
        self.nodelist = []

        #return if rectangle does not contain any cells
        (xmin, ymin), (xmax, ymax) = grid.get_extent()
        bgrid = box(xmin, ymin, xmax, ymax)
        (xmin, ymin), (xmax, ymax) = rectangle
        b = box(xmin, ymin, xmax, ymax)
        if not b.intersects(bgrid):
            #return with nodelist as an empty list
            return
        
        Xe = grid.Xe
        Ye = grid.Ye
        
        jmin = ModflowGridIndices.find_position_in_array(Xe, xmin)
        if jmin is None:
            if xmin <= Xe[0]:
                jmin = 0
            elif xmin >= Xe[-1]:
                jmin = grid.ncol - 1
                
        jmax = ModflowGridIndices.find_position_in_array(Xe, xmax)
        if jmax is None:
            if xmax <= Xe[0]:
                jmax = 0
            elif xmax >= Xe[-1]:
                jmax = grid.ncol - 1

        imin = ModflowGridIndices.find_position_in_array(Ye, ymax)
        if imin is None:
            if ymax >= Ye[0]:
                imin = 0
            elif ymax <= Ye[-1]:
                imin = grid.nrow - 1
                
        imax = ModflowGridIndices.find_position_in_array(Ye, ymin)
        if imax is None:
            if ymin >= Ye[0]:
                imax = 0
            elif ymin <= Ye[-1]:
                imax = grid.nrow - 1

        for i in range(imin, imax + 1):
            for j in range(jmin, jmax + 1):
                self.nodelist.append( (i, j) )
        return


class PolygonGridIntersection(object):
    '''
    
    '''
    def __init__(self, grid, polygon, **kwargs):
        '''
        Find the nodes within the polygon.
        
        '''
        from gridutil import ModflowGridIndices
        from shapely.geometry import Point, Polygon

        holes = None
        if kwargs.has_key('holes'): holes = kwargs['holes']

        #initialize the result arrays
        self.nodelist = []
        self.areas = []
        self.containscentroid = []
        pg = Polygon(polygon, holes=holes)

        #set X and Y
        X = grid.X
        Y = grid.Y

        #use the bounds of the polygon to restrict the cell search
        minx, miny, maxx, maxy = pg.bounds
        rectangle = ((minx, miny), (maxx, maxy))
        nodelist = RectangleGridIntersection(grid, rectangle).nodelist
        
        for (i, j) in nodelist:
            nodenumber = ModflowGridIndices.nn0_from_kij(0, i, j,
                            grid.nrow, grid.ncol)
            #node = grid.get_nodeobj(nodenumber)
            #node_polygon = Polygon(node.vertices)
            node_polygon = Polygon(grid.get_vertices(nodenumber))
            if pg.intersects(node_polygon):
                area = pg.intersection(node_polygon).area
                if area > 0.:
                    self.nodelist.append( (i, j) )
                    self.areas.append(area)
                    #pt = Point(node.position)
                    #if pg.contains(pt):
                    #    self.containscentroid = True
                    #else:
                    #    self.containscentroid = False
        return


if __name__ == '__main__':
    from pylab import *
    try:
        close('all')
    except:
        pass

    from grid import *
    from plotting import *
    nlay = 3
    nrow = 3
    ncol = 3
    delr = 1.
    delc = 1.
    
    npr = 2
    npc = 2
    iplotnum = 1

    #2d
    mfg2d = ModflowGrid2D(nrow, ncol, delr, delc)
    subplot(npr, npc, iplotnum, aspect='equal')
    iplotnum += 1
    drawgrid(mfg2d)
    points = [
        (0., 3.), (.5, 3.), (1.5, 3.), (2.5, 3.), (3.0, 3.),
        (0., 2.5), (.5, 2.5), (1.5, 2.5), (2.5, 2.5), (3.0, 2.5),
        (0., 1.5), (.5, 1.5), (1.5, 1.5), (2.5, 1.5), (3.0, 1.5),
        (0., 0.5), (.5, 0.5), (1.5, 0.5), (2.5, 0.5), (3.0, 0.5),
        (0., 0.0), (.5, 0.0), (1.5, 0.0), (2.5, 0.0), (3.0, 0.0),
        ]
    for feature in points:
        featuretype = 'point'
        (i, j) = mfg2d.intersection(feature, featuretype)
        plot(feature[0], feature[1], 'bo')
        text(feature[0], feature[1], str( (i, j) ) )

    #3d
    mfg3d = ModflowGrid(nlay, nrow, ncol, delr, delc, numpy.arange(nlay, -1, -1))
    subplot(npr, npc, iplotnum, aspect='equal')
    iplotnum += 1
    drawgrid(mfg3d)
    z = 0.
    points = [
        (0., 3.0, z), (.5, 3.0, z), (1.5, 3.0, z), (2.5, 3.0, z), (3.0, 3.0, z),
        (0., 2.5, z), (.5, 2.5, z), (1.5, 2.5, z), (2.5, 2.5, z), (3.0, 2.5, z),
        (0., 1.5, z), (.5, 1.5, z), (1.5, 1.5, z), (2.5, 1.5, z), (3.0, 1.5, z),
        (0., 0.5, z), (.5, 0.5, z), (1.5, 0.5, z), (2.5, 0.5, z), (3.0, 0.5, z),
        (0., 0.0, z), (.5, 0.0, z), (1.5, 0.0, z), (2.5, 0.0, z), (3.0, 0.0, z),
        ]
    for feature in points:
        featuretype = 'point'
        (k, i, j) = mfg3d.intersection(feature, featuretype)
        plot(feature[0], feature[1], 'bo')
        text(feature[0], feature[1], str( (k, i, j) ) )

    subplot(npr, npc, iplotnum, aspect='equal')
    iplotnum += 1
    drawgrid(mfg2d)
    featuretype = 'line'
    a = zeros( (mfg2d.nodes), dtype=int)
    lines = [ [(-1, -1), (0.5, 0.5), (4.0, 0.5), (0.0, 3.5)] ]
    for feature in lines:
        nodes, lengths = mfg2d.intersection(feature, featuretype)
        for (i, j) in nodes:
            nodeid = mfg2d.get_nodeid(i, j)
            a[nodeid] = 1
        xl = []; yl = []
        for p in feature:
            xl.append(p[0]); yl.append(p[1])
        plot(xl, yl, 'r-')
    print (nodes, lengths)
    colorflood(mfg2d, a, alpha=0.2)

    show()
    draw()
    