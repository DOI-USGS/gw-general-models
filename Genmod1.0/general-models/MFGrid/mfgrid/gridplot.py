
import matplotlib.pyplot
import numpy
import grid
#import quadpatch

def draw(mfg, **kwargs):
    """
    Method for quickly drawing a 2D grid.  Adapted from:
        http://exnumerus.blogspot.com/2011/02/how-to-quickly-plot-multiple-line.html
    """
    ((xmin, ymin), (xmax, ymax)) = mfg.get_extent()
    xpairs = []
    ypairs = []
    for j in range(mfg.ncol + 1):
        xpairs.append((mfg.Xe[j], mfg.Xe[j]))
        ypairs.append((ymin, ymax))
    for i in range(mfg.nrow + 1):
        xpairs.append((xmin, xmax))
        ypairs.append((mfg.Ye[i], mfg.Ye[i]))
    xlist = []
    ylist = []
    for xends,yends in zip(xpairs,ypairs):
        xlist.extend(xends)
        xlist.append(None)
        ylist.extend(yends)
        ylist.append(None)        
    matplotlib.pyplot.plot(xlist, ylist, **kwargs)
    return

def label(mfg, arr, **kwargs):
    """
    Function for labeling a 2D grid.
    """
    for i in range(mfg.nrow):
        for j in range(mfg.ncol):
            if arr is None:
                slabel = str((i,j))
            else:
                slabel = arr[i, j]
            matplotlib.pyplot.text(mfg.X[j], mfg.Y[i], slabel, **kwargs)
    return

def colorflood(gridobject, arr, **kwargs):
    """
    Method for color flooding an array for a Modflow grid.
    """

    #imports
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    #set defaults from kwargs
    ax = matplotlib.pyplot.gca()
    arrexclude = []
    arrmin = arr.min()
    arrmax = arr.max()
    alpha = 1.0
    cb = False
    if kwargs.has_key('ax'): ax = kwargs['ax']
    if kwargs.has_key('arrexclude'): arrexclude = kwargs['arrexclude']
    if kwargs.has_key('arrmin'): arrmin = kwargs['arrmin']
    if kwargs.has_key('arrmax'): arrmax = kwargs['arrmax']
    if kwargs.has_key('alpha'): alpha = kwargs['alpha']
    if kwargs.has_key('cb'): cb = kwargs['cb'] #plot colorbar?

    norm = colors.normalize(arrmin, arrmax)
    patches = []
    for nodenumber in range(gridobject.nodes):
        nodeobj = gridobject.get_nodeobj(nodenumber)
        (k, i, j) = gridobject.get_indices(nodenumber)
        x, y = nodeobj.position[0], nodeobj.position[1]
        dx, dy = nodeobj.dxdydz[0], nodeobj.dxdydz[1]
        xmin = x - dx / 2.
        xmax = x + dx / 2.
        ymin = y - dy / 2.
        ymax = y + dy / 2.
        vertices = [ (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
        if arr[k, i, j] in arrexclude:
            #cannot get visible and alpha to work, so this is the hack.
            vertices = [(0,0), (0,0), (0,0), (0,0) ]
        poly = Polygon(vertices, linewidth=0., fc = None, 
                       ec='none')
        patches.append(poly)

    #set up the patchcollection for faster coloring
    p = PatchCollection(patches, match_original=True)
    p.set_array(arr.reshape(-1))
    p.set_norm(norm)
    p.set_edgecolor(None)
    p.set_alpha(alpha)
    ax.add_collection(p)
    
    #create the colorbar
    if cb:
        pc = []
        color = cm.jet(norm(arrmin))
        poly = Polygon( ((0,0), (1,0), (1,1), (0,1)), linewidth=0., 
            fc = color, ec='none', alpha=0.4)
        pc.append(poly)
        color = cm.jet(norm(arrmax))
        poly = Polygon( ((0,0), (1,0), (1,1), (0,1)), linewidth=0., 
            fc = color, ec='none', alpha=0.4)
        pc.append(poly)
        p = PatchCollection(pc, match_original=True)
        p.set_array(numpy.array([arrmin, arrmax]))
        matplotlib.pyplot.colorbar(p, shrink=0.5)

    return    

