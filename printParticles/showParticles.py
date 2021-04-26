import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import json
import os.path
import time

fileName = '_particles.csv'
configName = 'config.json'
dirname = os.path.join(os.path.dirname(__file__), 'data')

data_list = []

# import the config file:
if os.path.isfile(os.path.join(dirname, configName)):
    config = pd.read_json(os.path.join(dirname, configName), typ='series')
else:
    config = {'domainSize' : 10}

# import all the data:
count = 0
while os.path.isfile(os.path.join(dirname, str(count) + fileName)):
    data = pd.read_csv(os.path.join(dirname, str(count) + fileName))
    X = np.array(data['px'])
    Y = np.array(data['py'])
    data_list.append(np.array([X, Y]))
    count += 1

data_list = np.array(data_list)
print(count, " files imported.")

fig = plt.figure()
fig.set_size_inches((15, 15))
line, = plt.plot([], [], 'ro')
plt.xlim(0, config['domainSize'])
plt.ylim(0, config['domainSize'])
plt.title('Particles')

def init():
    global line
    line.set_data([],[])
    return line,

def animate(i):
    global data_list
    line.set_data(data_list[i, 0], data_list[i, 1])
    return line,

ani = FuncAnimation(fig, animate, init_func=init, frames=count, blit=True, interval=15, repeat=False)

# show
plt.show()