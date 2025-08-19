import numpy as np
import matplotlib.pyplot as plt


def plot_raster_simple(x, y, data, plotname, cmap='Spectral', vmin=None, vmax=None, figsize=(10, 8)):
    """
    General function to visualize a raster and save it into PNG.

    :param x: 1d array
    :param y: 1d array
    :param data: 2d array
    :param plotname: string, filename
    :param cmap: default 'Spectral'
    :param vmin: float, optional, None
    :param vmax: float, optional, None
    :param figsize: tuple of (xwidth, yheight) values
    :return:
    """
    vmax = vmax if vmax is not None else np.nanmax(data)
    vmin = vmin if vmin is not None else np.nanmin(data)
    plt.figure(figsize=figsize, dpi=300)
    X, Y = np.meshgrid(x, y)
    plt.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.savefig(plotname)
    plt.close()
    return
