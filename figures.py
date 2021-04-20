from matplotlib import pyplot as plt
import numpy as np


def makePlot(stack):
    fig = plt.figure(figsize=(10.12, 9))

    cols = len(stack) // 3
    rows = len(stack) % cols

    gs = fig.add_gridspec(rows, cols, hspace=0, wspace=0)
    axs = gs.subplots(sharex='all', sharey='all')

    i = 0
    for row in axs:
        for ax in row:
            mappable = ax.imshow(
                image, origin='lower',  cmap=plt.cm.hsv, vmin=-10, vmax=10)

            if sample_point:
                ax.plot(sample_point[1], sample_point[0],
                        'x', markersize=4, color='black')

            ax.set_aspect(7)
            ax.grid(color='gray', alpha=0.4)
            ax.text(.045, .95, dates[i].isoformat().split('T')[0], fontweight='bold', horizontalalignment='left',
                    verticalalignment='center', transform=ax.transAxes)
            i += 1

    ax = fig.add_subplot(111, frameon=False)
    # # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)
    plt.xlabel("Pixels (Range)", labelpad=10, fontdict={
               'fontsize': 8, 'weight': 'bold'})
    plt.ylabel("Pixels (Aziumth)", labelpad=10, fontdict={
               'fontsize': 8, 'weight': 'bold'})

    cbar = plt.colorbar(mappable, ax=axs)
    # cbar.ax.get_yaxis().labelpad = 100
    # cbar.ax.get_xaxis().labelpad = 250

    cbar.ax.set_ylabel('Phase Difference (Radians)', rotation=270)

    plt.savefig('./figures/crater_unwrapped.png', dpi=300)
    plt.show()
