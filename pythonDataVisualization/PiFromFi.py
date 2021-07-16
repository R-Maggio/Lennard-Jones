import json
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter

from utils import *

# f_ = np.array([0.03154873, 0.12836987, 0.03262291, 0.10908559, 0.44391194, 0.11282776, 0.02359311, 0.09598237, 0.02439005])
f = np.array([0.44391194, 0.12836987, 0.03154873, 0.10908559, 0.02359311, 0.09598237, 0.02439005, 0.11282776, 0.03262291])
domainSize = 1
Dt = 1
v = domainSize/Dt
# V_ = v* np.array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
#             [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
# W_ = np.array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

V = v * np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])
W = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

# permutation = np.array([4, 1, 0, 3, 6, 7, 8, 5, 2])
# print(V_[permutation] - V)
# print(W_[permutation] - W)
# print(f_[permutation] - f)

rho = compute_rho_fi(f)
j = compute_j_fi(f, V)
Pi = compute_Pi_fi(f, V)
# Pi_neq = Pi - compute_Pi_eq(rho, v, j/rho)

f_eq = compute_feq(rho, j, Pi, V, W, v, 1)
f_neq = compute_fneq2(rho, j, Pi, V, W, v)

f_moments = compute_fi_moments(rho, j, Pi, V, W, v, Dt, 1)

print("rho: ", rho)
print("j: ", j)
print("Pi", Pi)
# print("Pineq", Pi_neq)

print('\n')

print("f_eq", f_eq)
print("f_neq", f_neq)

print('\n')

print("f", f)
print("f_moments", f_moments)

N = 3600 # number of particles
print("N: ", N)
print("m: ", rho*domainSize**2 / N)

covariance = compute_covariance(rho, j, Pi)
print("covariance", covariance)

sigmaX = np.sqrt(covariance[0, 0])
sigmaY = np.sqrt(covariance[1, 1])
corr = covariance[0, 1]/(sigmaX*sigmaY)

print("sigmaX, sigmaY, corr", sigmaX, sigmaY, corr)

print("u: ", j/rho)

# print(Pi/np.array([[ 1.43049070e+22, -2.29096562e+19], [-2.29096562e+19, 1.39582919e+22]]))