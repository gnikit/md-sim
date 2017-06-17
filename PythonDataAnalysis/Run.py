from FileLoadingPlotting import *
from PathHandling import *


path = OSPaths()
path.dir('Density 0.5', '5000')

"""
For MSD call the method diffusion_plot()
"""
# path = "C:/Users/user/Desktop/MD simulation/Archives of Data"
# density = "/Density 0.5"
# # density = "/Constant Force"
# density = "/Density 0.8"
# path += density

# Changing the Directory
# Varying T 5k
# directory = path +"/Non-Isothermal~step 5000/"
#---------------------------------------------------------------------
# iso 2.5k
# directory = path + "/Isothermal~step 2500/"
# iso 5k
# directory = path + "/Isothermal~step 5000/"
# # iso 10k
# directory = path + "/Isothermal~step 10000/"
# # iso 12.5k
# directory = path + "/Isothermal~step 12500/"
# # iso 15k
# directory = path + "/Isothermal~step 15000/"
# # iso 20k
# directory = path + "/Isothermal~step 20000/"
# iso 5k steps 8k particles
# directory = path + "/Particles 8k/"
# iso 5-3k 2.744 particles
# directory = path + "/Particles 2,744k/Isothermal~step 5000/"
# Low temperatures
# directory  =path +"/Low Temperature/Isothermal~step 5000/"
# os.chdir(directory)
cd = ["/Isothermal~step 5000/", "/Isothermal~step 10000/", "/Isothermal~step 12500/",
      "/Isothermal~step 15000/", "/Isothermal~step 20000/"]
steps_num = ['5000', '10000', '12500', '15000', '20000']
#############################################################################
#
# potential(6,1), potential(12,1),potential(6,1.5), potential(12, 1.5)#, potential(6,3), potential(6,3.5)
# potential(12,1), potential(12,1.5), potential(12,1.75), potential(12, 2), potential(12,3), potential(12,3.5)
# force(6,1.5), force(6,1.75), force(6,2), force(6,3), force(6,3.5)
# force(8,1.5), force(8,1.75), force(8,2), force(8,3), force(8,3.5)
# force(12,1.5), force(12,1.75), force(12,2), force(12,3), force(12,3.5)

# force(6,1.5), force(150,1.5)
# radial_dist_func(6, 0), radial_dist_func(6, 1), radial_dist_func(6,1.5), radial_dist_func(6,2), radial_dist_func(6,3)
# particle_plot(6,0), particle_plot(6,0.75), particle_plot(6,1.5), particle_plot(6,4)
#
# avg_Kin(6), avg_Kin(7), avg_Kin(8), avg_Kin(9), avg_Kin(10), avg_Kin(12)
# avg_pressure(6), avg_pressure(8), avg_pressure(10), avg_pressure(12)
# avg_Pot(6), avg_Pot(7), avg_Pot(8), avg_Pot(9), avg_Pot(10), avg_Pot(11), avg_Pot(12)
# avg_en(11)
# energy_plots(6, 2)
# energy_plots(5,0.5)

A_list = [0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4.]
b_list = [5,6,7,8,9,10,11,12]

A_list = [0.75, 1.0, 1.25, 1.5]
b_list = [6, 8, 10, 12]

# for i in b_list:
#     # vel_autocor_func(6,A_list[i])
#     # radial_dist_func(i, 0.75)
#     # pc(6,i)   # A_list
#     # mean_sqr_disp(6, i) # A_list
#     diffusion_plot(i)
#
# for w in steps_num:  #  LINUX
#     temp = OSPaths()
#     temp.dir('Density 0.5', w)
#     for i in range(len(A_list)):
# #         vel_autocor_func(6, A_list[i])

# for w in range(len(cd)):
#     temp = path
#     temp += cd[w]
#     os.chdir(temp)
#     for i in range(len(A_list)):
#         vel_autocor_func(6, A_list[i])

# mean_sqr_disp(6,0),mean_sqr_disp(6,0.75),mean_sqr_disp(6,1.50),mean_sqr_disp(6,4)
# diffusion_plot(6), diffusion_plot(7),diffusion_plot(8),diffusion_plot(9), diffusion_plot(10), diffusion_plot(11),diffusion_plot(12)
# for i in range(6, 13, 2):
#     diffusion_plot(i)
# vel_autocor_func(6,0), vel_autocor_func(6,0.75), vel_autocor_func(6,1), vel_autocor_func(6,4)


# vector_field(6, 0.75),
# vector_field(6,4.00)
# vector_field(12,0.75)

# force(6,.5), force(8,.5),force(10,.5), force(12,.5)
## For A <= 1.0, force increases with n (increasing)
### For A > 1, force decreases with n (increasing)
# potential(6, 0.0), potential(8, 0.), potential(10, 0.), potential(12, 0.) # varying n
# RDF2(6, 0.5), RDF2(8, 0.5), RDF2(10, 0.5), RDF2(12, 0.5)
list_a = [0.75, 1, 1.25, 1.5]
list_n = [6, 8, 10, 12]
mean_sqr_disp(6, 0)
# for index in list_n:
#     radial_dist_func(index, 0)

plt.tight_layout()
plt.show()

