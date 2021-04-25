import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas
import os.path
import time

fileName = '_particles.csv'
dirname = os.path.join(os.path.dirname(__file__), 'data')

data_list = []

# import all the data:
count = 0

while os.path.isfile(os.path.join(dirname, str(count) + fileName)):
    data = pandas.read_csv(os.path.join(dirname, str(count) + fileName))
    X = np.array(data['px'])
    Y = np.array(data['py'])
    data_list.append(np.array([X, Y]))
    count += 1

data_list = np.array(data_list)
print(count, " files imported.")

fig = plt.figure()
fig.set_size_inches((15, 15))
line, = plt.plot([], [], 'ro')
plt.xlim(0, 10)
plt.ylim(0, 10)
plt.title('Particles')

def init():
    global line
    line.set_data([],[])
    return line,

def animate(i):
    global data_list
    line.set_data(data_list[i, 0], data_list[i, 1])
    return line,

ani = FuncAnimation(fig, animate, init_func=init, frames=count, blit=True, interval=500, repeat=False)

# show
plt.show()