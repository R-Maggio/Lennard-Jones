import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import json
import os
import time

from utils import *



fileName = '_particles.csv'
configName = 'config.json'
dirname = os.path.join(os.path.dirname(__file__), 'dataPiNeq')

data_position_in = []
data_position_out = []
data_velocity_in = []
data_velocity_out = []
data_mass = []
data_config = []

# import all the data:
count = 0
while os.path.isfile(os.path.join(dirname, str(count) + "_in" + fileName)):
    data_config.append(pd.read_json(os.path.join(dirname, str(count) + "_" + configName), typ='series'))
    data_in = pd.read_csv(os.path.join(dirname, str(count) + "_in" + fileName))
    data_out = pd.read_csv(os.path.join(dirname, str(count) + "_out" + fileName))
    X_in = np.array(data_in['px'])
    Y_in = np.array(data_in['py'])
    Vx_in = np.array(data_in['vx'])
    Vy_in = np.array(data_in['vy'])
    X_out = np.array(data_out['px'])
    Y_out = np.array(data_out['py'])
    Vx_out = np.array(data_out['vx'])
    Vy_out = np.array(data_out['vy'])
    Mass_list = np.array(data_in['m'])
    data_position_in.append(np.array([X_in, Y_in]))
    data_velocity_in.append(np.array([Vx_in, Vy_in]))
    data_position_out.append(np.array([X_out, Y_out]))
    data_velocity_out.append(np.array([Vx_out, Vy_out]))
    data_mass.append(Mass_list)
    count += 1

data_position_in = np.array(data_position_in)
data_velocity_in = np.array(data_velocity_in)
data_position_out = np.array(data_position_out)
data_velocity_out = np.array(data_velocity_out)
data_mass = np.array(data_mass)

print(3*count, " files imported.")

#*----------------------------------

Dt = data_config[0]['Dt']
v = data_config[0]['domainSizeX']/Dt
V = v * np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])
W = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
rho = compute_rho(data_mass[0], data_config[0]['domainSizeX'], data_config[0]['domainSizeY'])

visc_list = []
corr_list = []
for i in range(count):
    visc_list.append(data_config[i]['visc'])
    corr_list.append(data_config[i]['correlation'])

corr_list = np.array(corr_list)
visc_list = np.array(visc_list)

Pi_neq_in_list = [] # non equilibrium part of Pi
Pi_neq_out_list = []
tau_list = []

for i in range(data_velocity_in.shape[0]):
    tau_list.append(compute_tau(visc_list[i], Dt, v))
    Pi_in = compute_Pi(data_config[i]['domainSizeX'], data_config[i]['domainSizeY'], data_velocity_in[i][0], data_velocity_in[i][1])
    Pi_out = compute_Pi(data_config[i]['domainSizeX'], data_config[i]['domainSizeY'], data_velocity_out[i][0], data_velocity_out[i][1])
    Pi_neq_in_list.append(Pi_in - compute_Pi_eq(rho, v, np.mean(data_velocity_in[i], axis=1)))
    Pi_neq_out_list.append(Pi_out - compute_Pi_eq(rho, v, np.mean(data_velocity_out[i], axis=1)))
Pi_neq_in_list = np.array(Pi_neq_in_list)
Pi_neq_out_list = np.array(Pi_neq_out_list)
tau_list = np.array(tau_list)

Pi_neq_in_xy = Pi_neq_in_list[:, 0, 1]
Pi_neq_out_xy = Pi_neq_out_list[:, 0, 1]

# plot PI ..........................................
fig1 = plt.figure()
plt.plot(Pi_neq_in_xy, Pi_neq_out_xy, label=r"$\Pi_{xy}^{neq, out}$")
plt.xlabel(r"$\Pi_{xy}^{neq, in}$")
plt.ylabel(r"$\Pi_{xy}^{neq, out}$")
plt.title(r"$\Pi_{xy}^{neq, out}$")

# plot tau:
fig2 = plt.figure()
plt.plot(Pi_neq_in_xy, Pi_neq_out_xy/Pi_neq_in_xy, label=r"$\Pi_{xy}^{neq, out}/\Pi_{xy}^{neq, in}$")
plt.xlabel(r"$\Pi_{xy}^{neq, in}$")
plt.ylabel(r"$\Pi_{xy}^{neq, out}/\Pi_{xy}^{neq, in}$")
plt.title(r"$(1 - \frac{1}{\tau})$")

# on a smaller interval:
# plt.figure()
# plt.plot(Pi_neq_in_xy[34:50], Pi_neq_out_xy[34:50], label=r"$\Pi_{xy}^{neq, out}$")
# plt.xlabel(r"$\Pi_{xy}^{neq, in}$")
# plt.ylabel(r"$\Pi_{xy}^{neq, out}$")
# plt.title(r"$\Pi_{xy}^{neq, out}$")

# remove "outsider:"
# filtered_pi_out = []
# filtered_pi_in = []

# for i in range(Pi_neq_in_xy.shape[0]):
#     if np.abs(Pi_neq_out_xy[i]) < 0.025:
#         filtered_pi_out.append(Pi_neq_out_xy[i])
#         filtered_pi_in.append(Pi_neq_in_xy[i])
#     else:
#         print(corr_list[i])

# plt.figure()
# plt.plot(filtered_pi_in, filtered_pi_out, label=r"$\Pi_{xy}^{neq, out}$")
# plt.xlabel(r"$\Pi_{xy}^{neq, in}$")
# plt.ylabel(r"$\Pi_{xy}^{neq, out}$")
# plt.title(r"$\Pi_{xy}^{neq, out}$ without outsiders")

plt.show()


if True:
    img_save_path = "../../.backup_data/images"
    images_names = ["Pi_comp", "tau"]
    figures = [fig1, fig2]
    for i in range(len(figures)):
        figures[i].savefig(img_save_path+'/'+images_names[i]+".png", bbox_inches='tight', pad_inches=0.05, dpi=200)