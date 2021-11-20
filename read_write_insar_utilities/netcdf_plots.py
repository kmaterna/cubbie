"""
Netcdf plotting functions
For the InSAR library, only Netcdf3 and Netcdf4 files with PIXEL NODE REGISTRATION are valid.
The assumption for 2D Netcdf files is 3 variables.
"""

from matplotlib import pyplot as plt
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3, read_any_grd


def produce_output_contourf(netcdfname, plottitle, plotname, cblabel):
    # Read in the dataset
    [xread, yread, zread] = read_netcdf3(netcdfname);

    # Make a plot
    _fig = plt.figure(figsize=(7, 10));
    plt.contourf(xread, yread, zread)
    plt.title(plottitle);
    plt.gca().set_xlabel("Range", fontsize=16);
    plt.gca().set_ylabel("Azimuth", fontsize=16);
    cb = plt.colorbar();
    cb.set_label(cblabel, size=16);
    plt.savefig(plotname);
    plt.close();
    return;


def produce_output_plot(netcdfname, plottitle, plotname, cblabel, aspect=1.0, invert_yaxis=True, dot_points=None,
                        vmin=None, vmax=None, cmap='rainbow'):
    # Read in the dataset
    [_, _, zread] = read_any_grd(netcdfname);

    # Make a plot
    fig = plt.figure(figsize=(7, 10));
    _ax1 = fig.add_axes([0.0, 0.1, 0.9, 0.8]);
    if vmin is not None:
        plt.imshow(zread, aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax);
    else:
        plt.imshow(zread, aspect=aspect, cmap=cmap);
    if invert_yaxis:
        plt.gca().invert_yaxis()  # for imshow, rows get labeled in the downward direction
    # plt.gca().get_xaxis().set_ticks([]);
    # plt.gca().get_yaxis().set_ticks([]);
    if dot_points is not None:
        plt.plot(dot_points[0], dot_points[1], color='black', marker='*', markersize=10);
    plt.title(plottitle);
    plt.gca().set_xlabel("Range", fontsize=16);
    plt.gca().set_ylabel("Azimuth", fontsize=16);
    cb = plt.colorbar();
    cb.set_label(cblabel, size=16);
    plt.savefig(plotname);
    plt.close();
    return;
