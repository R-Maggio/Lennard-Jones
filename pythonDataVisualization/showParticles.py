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
    # print(compute_mean_kinetic_energy(data_velocity[i, 0], data_velocity[i, 1], data_mass[i]))
    return line,

ani = FuncAnimation(fig, animate, init_func=init, frames=count, blit=True, interval=15, repeat=False)
# save the animation:
writer = PillowWriter(fps=25)  
ani.save(dirname + "/animation.gif", writer=writer)  
ani.save()

# show
# plt.show()

# --------------- tests ----------------
Dt = config['Dt']
v = config['domainSizeX']/Dt
V = v * np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])
W = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
msd = compute_MSD(data_position[0][0], data_position[0][1], data_position[-1][0], data_position[-1][1])
D = compute_D(msd, Dt)
rho = compute_rho(data_mass[0], config['domainSizeX'], config['domainSizeY'])
# visc = compute_kinematic_visc(data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0], rho, D, config['sigma'])
visc = compute_kinematic_visc_auto_corr(config['domainSizeX'], config['domainSizeY'], data_position[0, 0], data_position[0, 1], data_velocity[0, 0], data_velocity[0, 1], data_position[-1, 0], data_position[-1, 1], data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0], rho, Dt)
print("visc", visc)
print("kinetic mean", compute_mean_kinetic_energy(data_velocity[-1, 0], data_velocity[-1, 1], data_mass[0]))
print("rho", compute_rho(data_mass[0], config['domainSizeX'], config['domainSizeY']))
print("number of particles: ", data_velocity[0, 1].shape[0])

tau = compute_tau(visc, Dt, v)
u_in = np.mean(data_velocity[0], axis=1)
u_out = np.mean(data_velocity[-1], axis=1)
jin = rho*u_in
jout = rho*u_out

Pi_list = []
Pi_neq_list = [] # non equilibrium part of Pi
for i in range(data_velocity.shape[0]):
    Pi = compute_Pi(config['domainSizeX'], config['domainSizeY'], data_velocity[i][0], data_velocity[i][1])
    Pi_list.append(Pi)
    Pi_neq_list.append(Pi - compute_Pi_eq(rho, v, np.mean(data_velocity[i], axis=1)))
Pi_list = np.array(Pi_list)
Pi_neq_list = np.array(Pi_neq_list)

Pi_xx = Pi_list[:, 0, 0]
Pi_xy = Pi_list[:, 0, 1]
Pi_yy = Pi_list[:, 1, 1]
X_range = range(Pi_xx.size)

# plot PI ..........................................
plt.figure()
plt.plot(X_range, Pi_xx, label=r"$\Pi_{xx}$")
plt.plot(X_range, Pi_xy, label=r"$\Pi_{xy}$")
plt.plot(X_range, Pi_yy, label=r"$\Pi_{yy}$")

plt.legend()

# plot Pi_neq
plt.figure()
# plt.plot(X_range, Pi_neq_list[:, 0, 0], label=r"$\Pi_{xx}^{neq}$")
plt.plot(X_range, Pi_neq_list[:, 0, 1], label=r"$\Pi_{xy}^{neq}$")
# plt.plot(X_range, Pi_neq_list[:, 1, 1], label=r"$\Pi_{yy}^{neq}$")

plt.legend()

Pi_in = Pi_list[0]
Pi_out = Pi_list[-1]

# tau = 1
feq = compute_feq(rho, jin, Pi_in, V, W, v, tau)
fin = compute_fi_moments(rho, jin, 500*Pi_in, V, W, v, Dt, tau)
fout_LJ = compute_fi_moments(rho, jout, Pi_out, V, W, v, Dt, tau)
fout_LB = compute_fout_LB(fin, tau, feq)

print("\nfeq:\n",feq)
print("\nfin:\n",fin)
print("\nfout_LJ:\n", fout_LJ)
print("\nfout_LB:\n", fout_LB)

print("Pi_neq_in", Pi_neq_list[0, 0, 1])
print("Pi_neq_out", Pi_neq_list[-1, 0, 1])


# show ..........................................
if True:
    plt.show()