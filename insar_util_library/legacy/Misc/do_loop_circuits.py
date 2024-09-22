# This script forms loops of 3 interferograms, and then analyzes their phase unwrapping errors. 

# Step 1: Glob the intf_all
# Step 2: Form as many loops as possible
# Step 3: Make images of all the possible loops. 

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from Tectonic_Utils.read_write import netcdf_read_write

# BROKE AT LOOP # 120 OR # 121. NOT SURE WHY.
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3


def identify_all_loops():
    # This function takes the glob intf_all directories and then makes all possible triangles.
    # It populates a list of loops, each containing three images.
    loops, edges, nodes = [], [], []
    directories = glob.glob("intf_all/???????_???????")
    for item in directories:
        edgepair = item.split('/')[-1]
        edges.append(edgepair)
        nodes.append(edgepair.split('_')[0])
        nodes.append(edgepair.split('_')[1])
    nodes = list(set(nodes))

    # Making a graph of nodes and edges
    startnode = []
    endnode = []
    for i in range(len(nodes)):
        startnode.append(nodes[i])
        thisend = []
        for j in range(len(nodes)):
            if nodes[i] + '_' + nodes[j] in edges or nodes[j] + '_' + nodes[i] in edges:
                thisend.append(nodes[j])
        endnode.append(thisend)

    graph = dict([(startnode[i], endnode[i]) for i in range(len(startnode))])

    # Now we find loops.
    # Start in chronological order. Step through the main loop.
    # Within, have a loop of connectors, and another loop of connectors. If the innermost
    # layer comes back to the starter, then you have a loop of 3.
    for node0 in startnode:
        for node1 in graph[node0]:
            possible_list = graph[node1]
            for node2 in possible_list:
                if node2 == node0:
                    continue
                elif node0 in graph[node2]:
                    loop_try = [node0, node1, node2]
                    loop_try = sorted(loop_try)
                    if loop_try not in loops:
                        loops.append(loop_try)

    print("Number of Images (nodes): % s" % (len(startnode)))
    print("Number of Loops (circuits): % s" % (len(loops)))
    return loops


def compute_loops(all_loops, loops_dir, loops_guide, rowref, colref):
    os.makedirs(loops_dir, exist_ok=True)
    ofile = open(loops_dir + loops_guide, 'w')
    for i in range(len(all_loops)):
        ofile.write("Loop %d: %s %s %s\n" % (i, all_loops[i][0], all_loops[i][1], all_loops[i][2]))
    ofile.close()

    unwrapped = 'unwrap.grd'
    wrapped = 'phasefilt.grd'
    filename = os.path.join('intf_all', all_loops[0][0] + '_' + all_loops[0][1], unwrapped)
    z1_sample = read_netcdf3(filename)[2]
    number_of_errors = np.zeros(np.shape(z1_sample))

    for i in range(0, len(all_loops)):
        edge1 = all_loops[i][0] + '_' + all_loops[i][1]
        edge2 = all_loops[i][1] + '_' + all_loops[i][2]
        edge3 = all_loops[i][0] + '_' + all_loops[i][2]
        [_, _, z1] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge1, unwrapped))
        [_, _, z2] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge2, unwrapped))
        [_, _, z3] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge3, unwrapped))

        [xdata, ydata, wr_z1] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge1, wrapped))
        [_, _, wr_z2] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge2, wrapped))
        [_, _, wr_z3] = netcdf_read_write.read_any_grd(os.path.join('intf_all', edge3, wrapped))

        print("Loop " + str(i) + ":")

        rowdim, coldim = np.shape(z1)

        histdata_raw, histdata_fix = [], []
        znew_raw = np.zeros(np.shape(z1))
        znew_fix = np.zeros(np.shape(z1))
        errorcount = 0

        for j in range(rowdim):
            for k in range(coldim):

                wr1 = wr_z1[j][k] - wr_z1[rowref, colref]
                z1_adj = z1[j][k] - z1[rowref, colref]
                wr2 = wr_z2[j][k] - wr_z2[rowref, colref]
                z2_adj = z2[j][k] - z2[rowref, colref]
                wr3 = wr_z3[j][k] - wr_z3[rowref, colref]
                z3_adj = z3[j][k] - z3[rowref, colref]

                # Using equation from Heresh Fattahi's PhD thesis to isolate unwrapping errors.
                wrapped_closure_raw = wr_z1[j][k] + wr_z2[j][k] - wr_z3[j][k]
                wrapped_closure_fix = wr1 + wr2 - wr3
                offset_before_unwrapping = np.mod(wrapped_closure_fix, 2 * np.pi)
                if offset_before_unwrapping > np.pi:
                    offset_before_unwrapping = offset_before_unwrapping - 2 * np.pi  # send it to the -pi to pi realm.

                unwrapped_closure_raw = z1[j][k] + z2[j][k] - z3[j][k]
                unwrapped_closure_fix = z1_adj + z2_adj - z3_adj

                znew_raw[j][k] = unwrapped_closure_raw - offset_before_unwrapping
                znew_fix[j][k] = unwrapped_closure_fix - offset_before_unwrapping

                if ~np.isnan(znew_raw[j][k]):
                    histdata_raw.append(znew_raw[j][k] / np.pi)
                if ~np.isnan(znew_fix[j][k]):
                    histdata_fix.append(znew_fix[j][k] / np.pi)
                if abs(znew_fix[j][k]) > 0.5:  # if this pixel has
                    errorcount = errorcount + 1
                    number_of_errors[j][k] = number_of_errors[j][k] + 1

        errorpixels = round(100 * float(errorcount) / len(histdata_fix), 2)
        print("Most common raw loop sum: ")
        print(np.median(histdata_raw))
        print("Most common fix loop sum: ")
        print(np.median(histdata_fix))
        print("\n")

        make_plot(xdata, ydata, znew_fix, loops_dir + 'phase_closure_' + str(i) + '.eps', errorpixels)
        make_histogram(histdata_fix, loops_dir + 'histogram_' + str(i) + '.eps')

    return [xdata, ydata, number_of_errors]


def make_plot(x, y, z, plotname, errorpixels):
    fig = plt.figure()
    plt.imshow(z, cmap='jet', aspect=0.7)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    plt.title('Phase Closure: %.2f %% error pixels' % errorpixels)
    plt.gca().set_xlabel("Range", fontsize=16)
    plt.gca().set_ylabel("Azimuth", fontsize=16)
    cb = plt.colorbar()
    cb.set_label('phase residual', size=16)
    plt.savefig(plotname)
    plt.close()
    return


def make_histogram(histdata, plotname):
    plt.figure()
    plt.hist(histdata, bins=700)
    plt.yscale('log')
    plt.xlabel('Phase Difference (#*pi radians)')
    plt.ylabel('number of pixels')
    plt.savefig(plotname)
    plt.close()
    return


if __name__ == "__main__":
    loops_dir = "Phase_Circuits/"
    loops_guide = "loops.txt"
    rowref = 237  # doing correction for phase ambiguity
    colref = 172
    # reference_pixel = []  # not doing any correction for phase ambiguity

    all_loops = identify_all_loops()
    [xdata, ydata, number_of_errors] = compute_loops(all_loops, loops_dir, loops_guide, rowref, colref)

    # Print how often phase unwrapping errors affect different pixels.
    outfile = loops_dir + "how_many_errors.grd"
    netcdf_read_write.produce_output_netcdf(xdata, ydata, number_of_errors, 'phase unwrapping errors', outfile)
    netcdf_read_write.flip_if_necessary(outfile)
    make_plot(xdata, ydata, number_of_errors, loops_dir + 'Number of Unwrapping Errors', 0.00)
