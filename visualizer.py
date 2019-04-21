from utilities import dump_coords
import numpy as np

def plot_shape(shape, plt, c = 'k', lw=1, alpha = 1.0, fill=False):
    gtype = shape.geom_type.lower()
    if(gtype.startswith("multi") or gtype == "geometrycollection"):
        for part in shape:
            plot_shape(part, plt, fill=fill)
    else:
        if(gtype=="point"):
            plt.plot(shape.x, shape.y,'o'+c)
        if(gtype=='linestring'):
            coords = np.asarray(shape.coords)
            plt.plot(coords[:,0], coords[:,1],"-"+c, linewidth=lw, alpha=alpha)
        if(gtype=='polygon'):
            coords = np.asarray(shape.exterior.coords)
            plt.plot(coords[:,0], coords[:,1],"-"+c, linewidth=lw, alpha=alpha)
            if(fill is True):
                    plt.fill(coords[:,0], coords[:,1], 'gray')
            for interior in shape.interiors:
                coords = np.asarray(interior.coords)
                plt.plot(coords[:,0], coords[:,1],"-"+c, linewidth=lw, alpha = alpha)
                if(fill is True):
                    plt.fill(coords[:,0], coords[:,1], 'gray')

def annote_vertex(shape, plt):
    coords = dump_coords(shape)
    for i in range(len(coords)):
        plt.plot(coords[i][0], coords[i][1], 'ok', ms=2)
        plt.annotate(str(i), coords[i])