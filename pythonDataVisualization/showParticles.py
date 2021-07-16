import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import pandas as pd
import json
import os
import time

from utils import *


fileName = '_particles.csv'
configName = 'config.json'
dirname = os.path.join(os.path.dirname(__file__), 'data')

data_position = []
data_velocity = []
data_mass = []

# import the config file:
if os.path.isfile(os.path.join(dirname, configName)):
    config = pd.read_json(os.path.join(dirname, configName), typ='series')
else:
    config = {'domainSizeX' : 10, 'domainSizeY' : 10, 'sigma' : 0.3, 'DT': 1}

# import all the data:
count = 0
while os.path.isfile(os.path.join(dirname, str(count) + fileName)):
    data = pd.read_csv(os.path.join(dirname, str(count) + fileName))
    X = np.array(data['px'])
    Y = np.array(data['py'])
    Vx = np.array(data['vx'])
    Vy = np.array(data['vy'])
    Mass_list = np.array(data['m'])
    data_position.append(np.array([X, Y]))
    data_velocity.append(np.array([Vx, Vy]))
    data_mass.append(Mass_list)
    count += 1

data_position = np.array(data_position)
data_mass = np.array(data_mass)
data_velocity = np.array(data_velocity)

print(count, " files imported.")

fig1 = plt.figure()
fig1.set_size_inches((15, 15))
plt.plot(data_position[0, 0], data_position[0, 1], 'ro')
plt.xlim(0, config['domainSizeX'])
plt.ylim(0, config['domainSizeY'])
plt.title('Initial Position')

fig = plt.figure()
fig.set_size_inches((15, 15))
line, = plt.plot([], [], 'ro')
plt.xlim(0, config['domainSizeX'])
plt.ylim(0, config['domainSizeY'])
plt.title('Particles')

def init():
    global line
    line.set_data([],[])
    return line,

def animate(i):
    global data_position
    line.set_data(data_position[i, 0], data_position[i, 1])
    if i == 0:
        fig1 = fig
    # print(compute_mean_kinetic_energy(data_velocity[i, 0], data_velocity[i, 1], data_mass[i]))
    return line,

ani = FuncAnimation(fig, animate, init_func=init, frames=count, blit=True, interval=15, repeat=False)
# save the animation:
# writer = PillowWriter(fps=25)  
# ani.save(dirname + "/animation.gif", writer=writer)  
# ani.save()

# show
# plt.show()

# --------------- tests ----------------
steps = config['steps']
Dt = config['Dt']
v = config['domainSizeX']/Dt
E = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])
V = v * E
W = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
msd = compute_MSD(data_position[0][0], data_position[0][1], data_position[-1][0], data_position[-1][1])
D = compute_D(msd, Dt)
rho = compute_rho(data_mass[0], config['domainSizeX'], config['domainSizeY'])
# visc = compute_kinematic_visc(data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0], rho, D, config['sigma'])
# visc = compute_kinematic_visc_auto_corr(config['domainSizeX'], config['domainSizeY'], data_position[0, 0], data_position[0, 1], data_velocity[0, 0], data_velocity[0, 1], data_position[-1, 0], data_position[-1, 1], data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0], rho, Dt)
visc = config['visc']
tau = compute_tau(visc, Dt, v)
print("visc", visc)
print("tau")
print("kinetic mean", compute_mean_kinetic_energy(data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0]))
print("rho", compute_rho(data_mass[0], config['domainSizeX'], config['domainSizeY']))
print("number of particles: ", data_velocity[0, 1].shape[0])


u_in = np.mean(data_velocity[0], axis=1)
u_out = np.mean(data_velocity[-1], axis=1)
jin = rho*u_in
jout = rho*u_out

Pi_list = []
Pi_neq_list = [] # non equilibrium part of Pi
for i in range(data_velocity.shape[0]):
    Pi = compute_Pi(config['domainSizeX'], config['domainSizeY'], data_velocity[i][0], data_velocity[i][1], data_mass[i])
    Pi_list.append(Pi)
    Pi_neq_list.append(Pi - compute_Pi_eq(rho, v, np.mean(data_velocity[i], axis=1)))
Pi_list = np.array(Pi_list)
Pi_neq_list = np.array(Pi_neq_list)

Pi_xx = Pi_list[:, 0, 0]
Pi_xy = Pi_list[:, 0, 1]
Pi_yy = Pi_list[:, 1, 1]
X_range = np.array(range(Pi_xx.size)) * (steps//Pi_xx.size)

# plot PI ..........................................
fig2 = plt.figure()
plt.plot(X_range, Pi_xx, label=r"$\Pi_{xx}$")
plt.plot(X_range, Pi_xy, label=r"$\Pi_{xy}$")
plt.plot(X_range, Pi_yy, label=r"$\Pi_{yy}$")
plt.xlabel(r"$iteration$")
plt.ylabel(r"$\Pi$")
plt.title('Pi')

plt.legend()

# plot Pi_neq
fig3 = plt.figure()
fig3.set_size_inches((8, 5))
# plt.plot(X_range, Pi_neq_list[:, 0, 0], label=r"$\Pi_{xx}^{neq}$")
plt.plot(X_range, Pi_neq_list[:, 0, 1], label=r"$\Pi_{xy}^{neq}$")
plt.xlabel(r"$iteration$")
plt.ylabel(r"$\Pi_{xy}^{neq}$")
plt.title(r'$\Pi_{xy}^{neq}$')
# plt.plot(X_range, Pi_neq_list[:, 1, 1], label=r"$\Pi_{yy}^{neq}$")

plt.legend()

Pi_in = Pi_list[0]
Pi_out = Pi_list[-1]

# tau = 1
feq = compute_feq(rho, jin, Pi_in, V, W, v, tau)
feq_list = np.array([feq for _ in range(Pi_list.shape[0])])
print("test\n", rho, rho*np.mean(data_velocity[0], axis=1), Pi_list[0], V, W, v)
f_list = np.array([compute_fi_moments(rho, rho*np.mean(data_velocity[i], axis=1), Pi_list[i], E, W, v, Dt, tau) for i in range(Pi_list.shape[0])])
fin = f_list[0]
fout_LJ = f_list[-1]
fout_LB = compute_fout_LB(fin, tau, feq)
# tau with pineq:
tau_pineq = 1/(1 - Pi_neq_list[-1, 0, 1]/Pi_neq_list[0, 0, 1])
fout_LB_tau_pineq = compute_fout_LB(fin, tau_pineq, feq)


# print fi
fig4 = plt.figure()
plt.plot(X_range, f_list[:,0], label=r"$f_0$")
# plt.plot(X_range, f_list[:,1], label=r"$f_1")
# plt.plot(X_range, f_list[:,2], label=r"$f_2")
# plt.plot(X_range, f_list[:,3], label=r"$f_3")
# plt.plot(X_range, f_list[:,4], label=r"$f_4")
# plt.plot(X_range, f_list[:,5], label=r"$f_5")
# plt.plot(X_range, f_list[:,6], label=r"$f_6")
# plt.plot(X_range, f_list[:,7], label=r"$f_7")
# plt.plot(X_range, f_list[:,8], label=r"$f_8")
plt.plot(X_range, feq_list[:,0], '--', label=r"$f^{eq}_0$")
# plt.plot(X_range, feq_list[:,1], '--', label=r"$f^{eq}_1")
# plt.plot(X_range, feq_list[:,2], '--', label=r"$f^{eq}_2")
# plt.plot(X_range, feq_list[:,3], '--', label=r"$f^{eq}_3")
# plt.plot(X_range, feq_list[:,4], '--', label=r"$f^{eq}_4")
# plt.plot(X_range, feq_list[:,5], '--', label=r"$f^{eq}_5")
# plt.plot(X_range, feq_list[:,6], '--', label=r"$f^{eq}_6")
# plt.plot(X_range, feq_list[:,7], '--', label=r"$f^{eq}_7")
# plt.plot(X_range, feq_list[:,8], '--', label=r"$f^{eq}_8")
plt.xlabel(r"$iteration$")
plt.ylabel(r"$f_i$")
plt.title('f_i')
plt.legend()

print("jin", jin)
print("jout", jout)
print("Pi_in", Pi_in)
print("Pi_out", Pi_out)

print("\ntau_pineq", tau_pineq)
print("\nfeq:\n",feq)
print("\nfin:\n",fin)
print("\nfout_LJ:\n", fout_LJ)
print("\nfout_LB:\n", fout_LB)
print("\nfout_LB_tau_pineq:\n", fout_LB_tau_pineq)

print("Pi_neq_in", Pi_neq_list[0, 0, 1])
print("Pi_neq_out", Pi_neq_list[-1, 0, 1])


cov_mat = compute_covariance(rho, jin, Pi_in)

nbParticles = 3600;
sigmaX = 0.5773866671701394
sigmaY = 0.5772862749939868
corr = -1.5574669947325287e-05
muX = 0.04846295
muY = -0.00560023

cov2 = np.array([[sigmaX**2, sigmaX*sigmaY*corr],[sigmaX*sigmaY*corr, sigmaY**2]])

print("covariance:\n", cov_mat)
print("cov2\n", cov2)
print("cov numpy:\n", np.cov(np.stack((data_velocity[0, 0], data_velocity[0, 1]), axis=0)))

print(cov2[0,0] + cov2[1,1] + muX**2 + muY**2)

# show ..........................................
if True:
    plt.show()

# save images:
# show
plt.show()

#  save images:
if True:
    img_save_path = "../../.backup_data/images"
    images_names = ["initial_position", "Pi", "Pi_neq", "fi"]
    figures = [fig1, fig2, fig3]
    for i in range(len(figures)):
        figures[i].savefig(img_save_path+'/'+images_names[i]+".png", bbox_inches='tight', pad_inches=0.05, dpi=200)