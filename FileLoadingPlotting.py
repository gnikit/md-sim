import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os  # directory change
import gc
from mpl_toolkits.mplot3d import Axes3D

# Remember to change directory to where the data file is
# in this case:
# os.getcwd() #gets directory
#
# def directory_change (potential, parameter_a):
#
#     """
#        This Function changes the directory in order to access generically
#      all the data files produces by the MD Simulation
#
#      int potential:
#         The Potential energy exponent (power) needed to be accessed
#
#      int parameter_a:
#         The Parameter A of the Potential energy
#     """
#
#     a = "Potential " + potential + "/"
#     b = "A = " + parameter_a
#     directory = "C:/Users/user/Dropbox/University/3.1/PH 3110 Project/sim v2.1/sim v2.1/"+a+b
#     return os.chdir(directory)

#power = raw_input("Enter the Potential you want to investigate: ")
#A = raw_input("Enter parameter A: ")

#directory_change(power, A)

# PrC = "PressureC"+power+".txt"
# PrK = "PressureK"+power+".txt"
# PrCK = "PCK"+power+".txt"

sep = "~"
dif_coef = np.array([])
reduced_dif_coef = np.array([])
dif_err = np.array([])
reduced_dif_err = np.array([])
dif_y_int = np.array([])
reduced_dif_y_int = np.array([])


# Colors used for plots
color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                  '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                  '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
color_sequence2 = ['#1f77b4', '#ff7f0e', '#2ca02c',
                  '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                  '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                  '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
# This is an iterator for the color array
p, c = 0, 0
j = 0
v = 0
s = 0

def energy_plots(power, par_a):
    global k
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    KinEn = "KinEn" + power_str + sep+ A+".txt"
    PotEn = "PotEn" + power_str + sep+ A+".txt"
    TotEn = "TotEn" + power_str + sep+ A+ ".txt"

    #  Loads the files from the dir
    num_lines = sum(1 for line in open(KinEn))  # Calculates the num of lines in a file
    KE = np.loadtxt(KinEn, delimiter="\n")
    U = np.loadtxt(PotEn, delimiter="\n")
    Tot = np.loadtxt(TotEn, delimiter="\n")

    #  Plots the Energies
    step = 0.005
    time = num_lines*step
    x = np.linspace(0, time, num_lines)
    fig = plt.figure()

    Kin = plt.subplot2grid((3, 2), (0, 0), colspan=1)
    Pot = plt.subplot2grid((3, 2), (1, 0), colspan=1)
    TOT = plt.subplot2grid((3, 2), (2, 0), colspan=1)
    All = plt.subplot2grid((3, 2), (0, 1), rowspan=3)

    k,u,t = KE[500], U[500], Tot[500]
    Kin.plot(x, KE, 'r')
    Kin.locator_params(axis='y', nbins=4), Kin.set_ylim(ymax=4)
    Pot.plot(x, U, 'g')
    Pot.locator_params(axis='y', nbins=3)
    Pot.set_ylabel("Energy units", size=16)
    TOT.plot(x, Tot, 'b')
    TOT.locator_params(axis='y', nbins=4)
    TOT.set_ylim(ymax=6)

    x_r = time/2 -time/4
    # Kin.set_title('Individual Plots n = %d' %power, fontsize=17)
    Kin.set_title('Individual Plots', fontsize=17)
    # Kin.text(x_r, k*1.5, 'Kinetic Energy', fontsize=10)
    # Pot.text(x_r, u/1.2, 'Potential Energy', fontsize=10)
    # TOT.text(x_r, t*1.1, 'Total Energy', fontsize=10)
    # All.set_title('Total Energy, Kinetic and Potential',fontsize=17)
    All.set_title('Energy Contributions',fontsize=17)
    All.set_xlabel(r"Time $t$",fontsize=16)

    TOT.set_xlabel(r"Time $t$",fontsize=16)
    fig.subplots_adjust(hspace=0)

    for ax in [Kin, Pot]:
        plt.setp(ax.get_xticklabels(), visible=False)
        # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        ax.set_yticks(ax.get_yticks()[1:])

    All.plot(x, KE, 'r', x, U, 'g', x, Tot, 'b')
    All.set_ylim(ymax=5)
# 3D Plot
def particle_plot(power, par_a):

    # Changes the directory
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    RX = "Loadrx"+power_str+ sep+ A+ ".txt"
    RY = "Loadry"+power_str+ sep+ A+ ".txt"
    RZ = "Loadrz"+power_str+ sep+ A+ ".txt"

    rx = np.loadtxt(RX, delimiter='\n')
    ry = np.loadtxt(RY, delimiter='\n')
    rz = np.loadtxt(RZ, delimiter='\n')
    col = np.sqrt(rx**2+ry**2)
    fig = plt.figure(4)
    ax = fig.add_subplot(111, projection='3d')
    S = ax.scatter(rx, ry, rz, c =rz, cmap='gnuplot_r')
    plt.colorbar(S)

def vector_field(power, par_a):
    # Changes the directory
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    RX = "Loadrx" + power_str + sep + A + ".txt"
    RY = "Loadry" + power_str + sep + A + ".txt"
    RZ = "Loadrz" + power_str + sep + A + ".txt"
    VX = "Loadvx" + power_str + sep + A + ".txt"
    VY = "Loadvy" + power_str + sep + A + ".txt"
    VZ = "Loadvz" + power_str + sep + A + ".txt"

    rx = np.loadtxt(RX, delimiter='\n')
    ry = np.loadtxt(RY, delimiter='\n')
    rz = np.loadtxt(RZ, delimiter='\n')
    vx = np.loadtxt(VX, delimiter='\n')
    vy = np.loadtxt(VY, delimiter='\n')
    vz = np.loadtxt(VZ, delimiter='\n')

    # num_lines = sum(1 for line in open(RX))  # Calculates the num of lines in a file
    # plt.figure(1)
    # d = np.linspace(0,num_lines,num_lines)
    # plt.plot(d,vz)
    # erx = np.sqrt(vx**2+vz**2)
    # ery = np.sqrt(vy**3+vz**2)
    name = "n: " + power_str + "  A: " + A
    Q = plt.quiver(rx,ry,vx,vy, rz,label=name, cmap='gnuplot_r')
    #plt.scatter(rx,ry, s=4, c='black')
    plt.colorbar(Q)
    # plt.xlim(xmax=10.05)
    # plt.ylim(ymax=10.05)
    # plt.title('Velocity vector field for n = %d' %power)
    plt.legend(loc="upper right")

# RDF Histogram
def radial_dist_func(power, par_a):
    global c
    power_str = str(int(power))
    A = "{:.2f}".format(par_a)
    HIST = "Hist" + power_str + sep+A+".txt"
    num_lines = sum(1 for line in open(HIST))


    Hist = np.loadtxt(HIST, skiprows=1, delimiter="\n")
    Hist = np.delete(Hist, 99)

    x = np.linspace(0, num_lines, num_lines-2)
    plt.figure(2)
    n, bins, patches = plt.hist(Hist, 100, color=color_sequence[c])
    plt.xlabel(r"$g(r)$", fontsize=18)
    plt.ylabel(r"Number of particles", fontsize=18)
    # hist = "hist"+power_str+"a"+A+".png"
    path = "C:/Users/user/Desktop/MD simulation/Report/Figures/"
    #hist = path+hist
   # plt.savefig(hist)
    ###########################################################################
    N, rho, Nhist = 8**3, 0.5, 100
    rg = ((N/rho)**(1./3.))/2.
    dr = rg/Nhist
    pwr = float(power)
    force_max = (par_a/(1.+pwr))**(0.5)
    x = np.multiply(x,dr)
    name = "n: " + power_str + "  A: " + A

    plt.figure(3)
    y = np.full((99),1,dtype=int)   # straight line at y=1
    plt.plot(x, y, '--',color='black')
    #xx = np.full((num_lines-2), force_max, dtype=float)
    #plt.plot(xx,Hist, '--')
    plt.plot(x, Hist, label=name,color=color_sequence[c])  # consider using xx
    plt.xlabel(r"Radius $r$", fontsize=18)
    plt.ylabel(r"$g(r)$", fontsize=18)
    #plt.title("Radial Distribution Function $RDF$", fontsize=20)
    plt.xlim(xmax=x[98])
    plt.legend(loc="lower right", borderpad=0.1, labelspacing=0.01, columnspacing=0.01,fancybox=True, fontsize=16)
    rdf = "rdf"+power_str+"a"+A+".png"
    rdf = path+ rdf
    #plt.savefig(rdf)
    c +=1
# VAF
def vel_autocor_func(power, par_a ):
    global p
    # Changes the directory
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    vaf = "VAF" + power_str + sep + A + ".txt"

    cr = np.loadtxt(vaf, delimiter='\n')
    num_lines = sum(1 for line in open(vaf))
    time = 0.005 * num_lines
    xxx = np.linspace(0, time, num=num_lines)
    name = ""
    if num_lines < 7000:
        name = "n: " + power_str + "  A: " + A

    plt.figure(5)
    y = np.full((num_lines), 0)
    xx = np.full((num_lines),time)
    yy = np.linspace(5,-0.5,num_lines)
    plt.plot(xx, yy,'--' , color = 'black')
    plt.plot(xxx, y, '--', color='black')
    plt.plot(xxx, cr, label=name)#,color=color_sequence2[p])
    # plt.title("Velocity Autocorrelation Functions (VAF) ")
    plt.xlabel(r"Time $t$", fontsize=18)
    plt.ylabel(r"$C_v$",fontsize=18)
    plt.ylim(ymax=5,ymin=-0.5)
    plt.legend(loc="upper right", ncol=1, borderpad=0.1, labelspacing=0.01, columnspacing=0.01,fancybox=True, fontsize=16)
    p += 1
# MSD
def mean_sqr_disp(power, par_a):
    global p, j
    global dif_coef, dif_err, dif_y_int
    global reduced_dif_coef, reduced_dif_err, reduced_dif_y_int
    # Changes the directory
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    msd = "MSD" + power_str + sep + A + ".txt"

    MSD = np.loadtxt(msd, delimiter='\n')

    num_lines = sum(1 for line in open(msd))
    limit = 0.005*num_lines
    step = int(0.6 * num_lines)
    x = np.linspace(0, limit, num=num_lines)

    if par_a >= 0:
        x_sliced = x[step:]
        MSD_sliced = MSD[step:]
        gradient, intercept, r_value, p_value, std_err = stats.linregress(x_sliced, MSD_sliced)

        reduced_dif_coef = np.append(reduced_dif_coef, gradient)
        reduced_dif_err = np.append(reduced_dif_err,std_err)
        reduced_dif_y_int = np.append(reduced_dif_y_int,intercept)

    # Regular coefs are calculated independent of the if loop
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x, MSD)
    dif_coef = np.append(dif_coef, gradient)
    dif_err = np.append(dif_err,std_err)
    dif_y_int = np.append(dif_y_int,intercept)

    print 'Diffusion coef: ',gradient,'\n',\
        'y-intercept: ',intercept,'\n',\
        'R value: ',r_value,'\n',\
        'Fit Error: ',std_err

    name = "n: " + power_str + "  A: " + A
    plt.figure(29)
    plt.plot(x, MSD, label=name, color=color_sequence[p])
    # plt.title("Mean Square Displacement (MSD)")
    plt.xlabel(r"Time $t$",fontsize=18)
    plt.ylabel(r"Mean Square Displacement $MSD$",fontsize=16)
    plt.legend(loc="upper left", fontsize=16, fancybox=True)
    p += 1
# Pressure C
def pc(power, par_a):
    global p
    # Changes the directory
    power_str = str(power)
    A = "{:.2f}".format(par_a)

    vaf = "PressureC" + power_str + sep + A + ".txt"

    cr = np.loadtxt(vaf, delimiter='\n')
    num_lines = sum(1 for line in open(vaf))
    time  = num_lines*0.005
    xxx = np.linspace(0, time, num=num_lines)

    name = "n: " + power_str + "  A: " + A
    plt.figure(5)

    plt.plot(xxx, cr, label=name)#,color=color_sequence[p])
    plt.xlabel(r"Time $t$", size=18)
    plt.ylabel(r"Configurational Pressure $P_C$", size=18)
    # plt.ylim(0,0.1)
    # plt.title("Pressure of Liquid against time", size=17)
    plt.legend(loc="upper right",prop={'size':12}, borderpad=0.2, labelspacing=0.2, handlelength=1)
    p += 1
# Averages
def avg_pressure(power):
    global p
    power_str = str(power)

    PC = "data" + power_str + ".txt"
    name = "n: " + power_str

    num_lines = sum(1 for line in open(PC))
    a, Pc = np.loadtxt(PC, delimiter='\t', comments='#', usecols=(0, 4), unpack=True)

    plt.figure(6)
    plt.plot(a, Pc, '-o',label=name)#, color=color_sequence[p])
    # plt.title("Configurational Pressure for multiple Potentials")
    plt.xlabel(r"Parameter $A$", size=18)
    plt.ylabel(r"Pressure $Pc$", size=18)
    plt.legend(loc="upper right")
    p +=1
def avg_Kin(power):
    # Changes the directory
    power_str = str(power)

    K = "data" + power_str + ".txt"
    name = "n: " + power_str

    a, k = np.loadtxt(K, delimiter='\t', comments='#', usecols=(0, 1), unpack=True)

    plt.figure(21)
    plt.plot(a, k, label=name)
    plt.title("Kinetic Energy vs multiple values of A")
    plt.xlabel(r"Parameter A")
    plt.ylabel(r"Kinetic Energy $K$")
    plt.legend(loc="lower right")
def avg_Pot(power):
    # Changes the directory
    power_str = str(power)

    K = "data" + power_str + ".txt"
    name = "n: " + power_str

    a, k = np.loadtxt(K, delimiter='\t', comments='#', usecols=(0, 2), unpack=True)

    plt.figure(22)
    plt.plot(a, k, label=name)
    plt.title("Potential Energy against parameter A")
    plt.xlabel(r"Parameter A")
    plt.ylabel(r"Potential Energy $U$")
    plt.legend(loc="upper right")
def avg_en(power):
    # Changes the directory
    power_str = str(power)

    E = "data" + power_str + ".txt"
    name = "n: " + power_str

    num_lines = sum(1 for line in open(E))

    a,K,U,T = np.loadtxt(E, delimiter='\t', comments='#', usecols=(0, 1,2,3), unpack=True)

    fig = plt.figure(111)


    Kin = plt.subplot2grid((3, 2), (0, 0), colspan=1)
    Pot = plt.subplot2grid((3, 2), (1, 0), colspan=1)
    TOT = plt.subplot2grid((3, 2), (2, 0), colspan=1)
    All = plt.subplot2grid((3, 2), (0, 1), rowspan=3)

    Kin.plot(a, K, color='r')
    Kin.set_ylim(ymin=2.0), Kin.locator_params(axis='y', nbins=5, prune='lower')
    Pot.plot(a, U, color='g')
    Pot.locator_params(axis='y', nbins=4), Pot.set_ylabel("Energy units")
    TOT.plot(a, T, color='c')
    TOT.locator_params(axis='y', nbins=3), TOT.locator_params(axis='x', nbins=4)
    TOT.set_xlabel(r"Parameter $A$")

    # k, u, t = K[4], U[4], T[4]
    # Kin.text(a[4], k*1.01, 'Average Kinetic Energy', fontsize=12)
    # Pot.text(a[4], u*0.95, 'Average Potential Energy', fontsize=12)
    # TOT.text(a[4], t, 'Average Total Energy', fontsize=12)
    Kin.set_title('Individual Plots for n = %d' %power, fontsize=14)
    All.set_title('Total Energy, Kinetic and Potential')
    All.set_xlabel(r"Parameter $A$")

    fig.subplots_adjust(hspace=0)

    for ax in [Kin, Pot]:
        plt.setp(ax.get_xticklabels(), visible=False)
        # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
        ax.set_yticks(ax.get_yticks()[1:])

    All.plot(a, K, 'r', a, U, 'g', a, T, 'c')
    All.set_ylim(ymax=6)
    All.locator_params(axis='x', nbins=4)

def diffusion_plot(power):
    global j, p, v, s
    global dif_coef, dif_err, dif_y_int
    global reduced_dif_coef, reduced_dif_err, reduced_dif_y_int

    my_list = [0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4.]
    my_list_reduced = [1.0, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4.]


    for i in range(0, len(my_list)):
        mean_sqr_disp(power, my_list[i])
        print "-----------------------------"

    # File saving Dif coe
    # np.savetxt("../../PlotAnalysis/Diffusion coefficients n=%d.txt"%power,dif_coef, fmt='%.5f')
    # np.savetxt("../../PlotAnalysis/Reduced Diffusion coefficients n=%d.txt"%power,reduced_dif_coef, fmt='%.5f')
    # np.savetxt("../../PlotAnalysis/dif error n=%d.txt"%power, dif_err, fmt='%.5f')
    # np.savetxt("../../PlotAnalysis/reduced dif error n=%d.txt"%power, reduced_dif_err, fmt='%.5f')
    # np.savetxt("../../PlotAnalysis/dif y int n=%d.txt"%power, dif_y_int, fmt='%.5f')
    # np.savetxt("../../PlotAnalysis/reduced dif y hint n=%d.txt"%power, reduced_dif_y_int, fmt='%.5f')
    name = "n: "+str(power)
    steps = ['5k','10k','12.5k','15k','20k']
    #name = steps[s]
    #s += 1

    plt.figure(66)
    plt.plot(my_list, dif_coef, '--o', label=name)#, color=color_sequence[v])
    #plt.plot(my_list, reduced_dif_coef[j:j+15], 'bo')
    #plt.title(r"Diffusion coefficient against $A$") #%power) # for n = %d
    plt.xlabel(r"Parameter $A$", fontsize=18)
    plt.ylabel(r"Diffusion coefficient $D$", fontsize=18)
    plt.legend(loc="lower right",title="Steps", fancybox=True, ncol=2, labelspacing=0.05, handletextpad=0.5,fontsize=16)

    dif_coef, dif_err, dif_y_int = np.array([]), np.array([]), np.array([])
    reduced_dif_coef, reduced_dif_err, reduced_dif_y_int = np.array([]), np.array([]), np.array([])

    p = 0
    j += 15
    v += 1

# No data plots
def potential(power, par_a):
    global p
    power_str = str(power)
    A = str(float(par_a))

    x = np.linspace(0.00,3,100)  # division with 0
    V = 1/pow((x**2 +par_a), power/2)
    plt.figure(6)
    name = "n: " + power_str + " A: " + A
    plt.plot(x,V, label=name, linestyle='-', color=color_sequence[p])
    # plt.title("Potential Energy curves", size=16)
    plt.xlabel(r"Separation distance $r$", size=16)
    plt.ylabel(r"Potential Energy $\Phi$",size=16)
    lgnd = plt.legend(loc="upper right", fancybox=True,ncol=1)
    lgnd.legendHandles[p]._legmarker.set_markersize(10) # not working...

    p += 1
def force(power, par_a):

    power_str = str(power)
    A = str(float(par_a))

    x = np.linspace(0.0,3,100)  # division with 0
    V = power*x*(pow((x**2 + par_a), (-power/2 -1)))
    plt.figure(7)
    name = "n: " + power_str + " A: " + A
    plt.plot(x,V, label=name)
    # plt.title("Pair Forces vs A", size=16)
    plt.xlabel(r"Separation Distance $r$", size=16)
    plt.ylabel(r"Force $f$", size=16)
    plt.legend(loc="upper right", fancybox=True,fontsize=16)
