

import matplotlib.pyplot as plt
import numpy as np


def before_after_images(before_phase, after_phase, outfilename):
    """
    Produce a comparison plot of before-and-after images of unwrapped phase.
    """
    plt.figure(figsize=(12, 7), dpi=300)
    plt.subplot(1, 2, 1)
    plt.imshow(before_phase, cmap='viridis', vmin=np.nanmin(before_phase), vmax=np.nanmax(after_phase))
    plt.title('Before Correction')
    plt.subplot(1, 2, 2)
    plt.imshow(after_phase, cmap='viridis', vmin=np.nanmin(before_phase), vmax=np.nanmax(after_phase))
    plt.title('After Correction')
    cb = plt.colorbar()
    cb.set_label("Unwrapped Phase (radians)", size=12)
    print("Saving figure %s " % outfilename)
    plt.savefig(outfilename)
    plt.close()


def linear_topo_phase_plot(phase_array_1d, dem_array_1d, corrected_phase_1d, outfilename):
    """
    Plot topography and phase as x-y scatter plot.
    """
    plt.figure(figsize=(8, 6), dpi=300)
    plt.plot(dem_array_1d, phase_array_1d, '.', label='Before correction')
    plt.plot(dem_array_1d, corrected_phase_1d, '.r', alpha=0.15, label='After correction')
    plt.ylabel('phase')
    plt.xlabel('topo')
    plt.title('Initial and Corrected Phase vs. Topography')
    plt.legend()
    print("Saving figure %s " % outfilename)
    plt.savefig(outfilename)
    plt.close()
    return
