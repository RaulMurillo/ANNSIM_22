import sys
import numpy as np
import matplotlib.pyplot as plt


def plot3d_lorenz(file, title=None, color=None):
    data = np.genfromtxt(file)
    assert data.shape[1]==4
    fig = plt.figure(figsize = (10,10))
    ax = plt.axes(projection='3d')
    ax.grid()
    ax.plot3D(data[:,1], data[:,2], data[:,3], color=color)
    if title is not None:
        ax.set_title(title)
        # plt.savefig('lorenz_'+title+'.pdf')
        plt.savefig('lorenz_'+title+'.png', bbox_inches='tight', dpi=400)
    else:
        # plt.savefig('lorenz_3d.pdf')
        plt.savefig('lorenz_3d.png', bbox_inches='tight', dpi=400)
    # plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 2:
        plot3d_lorenz(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        plot3d_lorenz(sys.argv[1])
