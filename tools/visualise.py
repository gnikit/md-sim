import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import importlib.util as util
from mdtools import FileNaming

# TODO: Use a more optimised file format/higher read/writes?
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

def command_line_input(argv):
    try:
        opts, args = getopt.getopt(argv, "hd:s:p:R:T:n:A:P:", [
                                   "dir=", "steps=", "particles=", "rho=", 
                                   "temp=", "n=", "A="])
    except getopt.GetoptError:
        sys.exit(2)

    out_dir = None
    steps = None
    particles = None
    rho = None
    temperature = None
    n = None
    par_a = None

    for opt, arg in opts:
        if opt == '-h':
            print("Display this message")
            sys.exit()
        elif opt in ("-d", "--dir"):
            out_dir = arg
        elif opt in ("-s", "--steps"):
            steps = arg
        elif opt in ("-p", "--particles"):
            particles = arg
        elif opt in ("-R", "--rho"):
            rho = arg
        elif opt in ("-T", "--temp"):
            temperature = arg
        elif opt in ("-n", "--n"):
            n = arg
        elif opt in ("-A", "--A"):
            par_a = arg
        else:
            print("Unrecognised option \"" + opt[0] + "\" entered")

    vars = f"\nInput variables\n" \
           f" Output file destination: {out_dir}\n" \
           f" Steps: {steps}\n" \
           f" Particles: {particles}\n" \
           f" -------------------------------------\n" \
           f" RHO: {rho}\n" \
           f" T: {temperature}\n" \
           f" n: {n}\n" \
           f" A: {par_a}\n"

    print(vars)
    return (out_dir, steps, particles, rho, temperature, n, par_a)

if __name__ == "__main__":

    args = command_line_input(sys.argv[1:])
    print(args)

    if (not all(i is None for i in args)):
        # Read file
        dir = str(args[0])
        STEPS = int(args[1])
        PARTICLES = int(args[2])
        RHO = float(args[3])
        TEMPERATURE = float(args[4])
        POWER = int(args[5])
        A = float(args[6])

        BOX_LENGTH = (PARTICLES / RHO) ** (1./3.0)

        file_id = FileNaming(STEPS, PARTICLES)
        name_id = file_id.file_searcher(RHO, TEMPERATURE, POWER, A)
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

    else:
        print("Error: Uninitialised variables in the arguments list.")  
