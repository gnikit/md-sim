import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
"""
    To use visualise, first enable VISUALISE = true in the MD object e.g.
    MD run (dir, steps);
    run.VISUALISE = true;
    // Now run MD simulation
    run.Simulation(...)

    Then change dir to the dir of the files
"""

# TODO: Make filename searcher
# TODO: Use a more optimised file format/higher read/writes?

# Read file
dir = "/home/gn/Desktop/test_data/sample"
PARTICLES = 1000
STEPS = 10000
RHO = 0.5
TEMPERATURE = 0.5
A = 0.5
POWER = 6
BOX_LENGTH = (PARTICLES / RHO) ** (1./3.0)


# Read files into vectors of vectors
def load_data(dir):
    # each line corresponds to a single step
    # hence each line has N particles (1000) by default
    x_all = np.loadtxt(f"{dir}/x_data.dat")
    y_all = np.loadtxt(f"{dir}/y_data.dat")
    z_all = np.loadtxt(f"{dir}/z_data.dat")

    return x_all, y_all, z_all


def update_figure(num):
    x = X_ALL[num]
    y = Y_ALL[num]
    z = Z_ALL[num]
    # TODO: NN calculation between all particles & -> assign colormap
    text.set_text(f"Frame: {num}")
    graph._offsets3d = (x, y, z)
    return graph

X_ALL, Y_ALL, Z_ALL = load_data(dir)

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection="3d")
graph = ax.scatter(X_ALL[0], Y_ALL[0], Z_ALL[0], color='orange')
text = fig.text(0, 1, "TEXT", va='top')  # for debugging

ax.set_xlim3d(0, BOX_LENGTH)
ax.set_ylim3d(0, BOX_LENGTH)
ax.set_zlim3d(0, BOX_LENGTH)

# Creating the Animation object
ani = animation.FuncAnimation(fig, update_figure, frames=10000, interval=0, blit=False)
# ani.save("/home/gn/Desktop/test_data/sample/animation.mp4", fps=60)
plt.show()