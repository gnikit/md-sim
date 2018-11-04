import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import importlib.util as util
sys.path.append("/home/gn/Code/MD-simulation/tools/MD-Simulation-Data-Analysis")
from md import FileNaming

"""
    MD run (dir, steps);
    run.VISUALISE = true;
    // Now run MD simulation
    run.Simulation(...)

    Then change dir to the dir of the files
"""

# TODO: Use a more optimised file format/higher read/writes?

# Read file
dir = "/home/gn/Code/MD-simulation/examples/example_data"
PARTICLES = 1000
STEPS = 15000
RHO = 0.5
TEMPERATURE = 0.5
A = 0.5
POWER = 6
BOX_LENGTH = (PARTICLES / RHO) ** (1./3.0)

file_id = FileNaming(STEPS, PARTICLES)
name_id = file_id.file_searcher(RHO, TEMPERATURE, POWER, A)
# Read files into vectors of vectors


def load_data(dir):
    # each line corresponds to a single step
    # hence each line has N particles (1000) by default

    x_all = np.loadtxt(f"{dir}/x_data{name_id}.txt")
    y_all = np.loadtxt(f"{dir}/y_data{name_id}.txt")
    z_all = np.loadtxt(f"{dir}/z_data{name_id}.txt")

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
ani = animation.FuncAnimation(
    fig, update_figure, frames=STEPS, interval=0, blit=False)
# ani.save("/home/gn/Desktop/test_data/sample/animation.mp4", fps=60)
plt.show()
