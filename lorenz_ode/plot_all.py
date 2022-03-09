import sys
import numpy as np
import matplotlib.pyplot as plt


def plot3d_lorenz(fig, nrows, ncols, index, file, title=None):
    data = np.genfromtxt(file)
    assert data.shape[1]==4
    # fig = plt.figure(figsize = (10,10))
    # ax = plt.axes(projection='3d')
    ax = fig.add_subplot(nrows, ncols, index, projection='3d')
    ax.grid()
    ax.plot3D(data[:,1], data[:,2], data[:,3])
    if title is not None:
        ax.set_title(title, fontdict = {'fontsize' : 20})
    #     plt.savefig('lorenz_'+title+'.pdf')
    # else:
    #     plt.savefig('lorenz_3d.pdf')
    # plt.show()


if __name__ == "__main__":
    # if len(sys.argv) > 2:
    #     plot3d_lorenz(sys.argv[1], sys.argv[2])
    # else:
    #     plot3d_lorenz(sys.argv[1])



    fig, axs = plt.subplots(2,3, figsize=(20,15)) # , sharey=True
    # fig.subplots_adjust(wspace=0, hspace=0)
    fig.subplots_adjust(hspace=0)


    index = 1

    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_double.txt", title="Float64"); index+=1
    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_float.txt", title="Float32"); index+=1
    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_float16.txt", title="Float16"); index+=1
    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_posit_32_2.txt", title="Posit32"); index+=1
    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_posit_20_2.txt", title="Posit20"); index+=1
    plot3d_lorenz(fig, 2, 3, index, "build/lorenz_ode_data_posit_16_2.txt", title="Posit16"); index+=1

    for ax in axs.flat:
        ax.axis('off')


    plt.tight_layout()
    plt.savefig('lorenz_3d.pdf', bbox_inches='tight')
    # plt.show()


