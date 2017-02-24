import os
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

# A_list = [0.00, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 10, 50, 100]
# A_list = [0.00, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 4.00]
A_list = [1.5, 2.75, 3, 3.5, 4]



for j in range(len(A_list)):

    index = '{:.2f}'.format(A_list[j])
    sep = "~"
    txt = ".txt"

    # directory = "C:/Users/user/Desktop/MD simulation/Archives of Data/Constant Force/Isothermal~step 5000/Particle Data/"
    directory = "D:/DATA/PLOTS/"
    os.chdir(directory)

    """
    Theoretically, the size of the ndarray is known since it depends on the
    number of steps and the number of particles N the simulation has
    """
    particles = 14 ** 3
    steps = 5000
    rho = 0.5
    L = pow((particles / rho), 1.0 / 3.0)

    lrx, lry, lrz = "Loadrx6" +sep+index+txt,"Loadry6" +sep+index+txt,"Loadrz6" +sep+index+txt
    lvx, lvy, lvz = "Loadvx6" +sep+index+txt,"Loadvy6" +sep+index+txt,"Loadvz6" +sep+index+txt
    rx = np.loadtxt(lrx, delimiter='\n')
    ry = np.loadtxt(lry, delimiter='\n')
    rz = np.loadtxt(lrz, delimiter='\n')
    vx = np.loadtxt(lvx, delimiter='\n')
    vy = np.loadtxt(lvy, delimiter='\n')

    np.ndarray.resize(rx, (steps, particles)),np.ndarray.resize(vx, (steps, particles))
    np.ndarray.resize(ry, (steps, particles)),np.ndarray.resize(vy, (steps, particles))
    np.ndarray.resize(rz, (steps,particles))#,np.ndarray.resize(vz, (steps,particles))

    p = 0
    q = 1

    while p < steps:
        x = rx[p][:]
        y = ry[p][:]
        z = rz[p][:]
        u = vx[p][:]
        v = vy[p][:]

        f1 = plt.figure(1)
        Q = plt.quiver(x, y, u, v, z, pivot='mid', color='red', cmap='Reds_r')
        plt.colorbar(Q)
        plt.gca()
        plt.xlim([0,18])
        plt.ylim([0,18])
        name = "image"
        num = str(p)
        num = num.zfill(5)
        name += num + '.png'

        a = "C:/DATA/6~"
        a += index +'/'
        os.chdir(a)
        p += 2
        f1.savefig(name, dpi=900)
        f1.clf()
########################################################################################################################
        xx = rx[q][:]
        yy = ry[q][:]
        zz = rz[q][:]
        uu = vx[q][:]
        vv = vy[q][:]
        f2 = plt.figure(2)
        V = plt.quiver(xx, yy, uu, vv, zz, pivot='mid', color='red', cmap='Reds_r')
        plt.colorbar(V)
        plt.gca()
        plt.xlim([0, 18])
        plt.ylim([0, 18])
        name = "image"
        num = str(q)
        num = num.zfill(5)
        name += num + '.png'

        a = "C:/DATA/6~"
        a += index + '/'
        os.chdir(a)
        q += 2
        f2.savefig(name, dpi=900)
        f2.clf()



