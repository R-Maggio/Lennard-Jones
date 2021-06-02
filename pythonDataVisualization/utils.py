import numpy as np
import matplotlib.pyplot as plt

def compute_strain_moments(rho, j, Pi, Dt, v, tau):
    S = np.zeros((2, 2))
    kronecker = lambda a, b: 1 if a == b else 0
    for a in range(S.shape[0]):
        for b in range(S.shape[1]):
            S[a, b] = (Pi[a, b] - rho * 1/3 * v**2 * kronecker(a, b) + 1/rho * j[a] * j[b]) / (-2/3 * Dt * tau * v**2 * rho)
    return S

def compute_Qi(Vi, v):
    Q = np.zeros((2, 2))
    kronecker = lambda a, b: 1 if a == b else 0
    for a in range(Q.shape[0]):
        for b in range(Q.shape[1]):
            Q[a, b] = Vi[a] * Vi[b] - 1/3 * v**2 * kronecker(a, b)
    return Q

def compute_feq(rho, j, Pi, V, W, v, tau):
    f = np.zeros(V.shape[0])
    kronecker = lambda a, b: 1 if a == b else 0
    for i in range(f.size):
        Qi = compute_Qi(V[i], v)
        f[i] = rho * W[i] * (1 + 3 * np.einsum('i,i', V[i], j) / (rho * v**2) + 9/(2 * v**4 * rho**2) * np.einsum('ij, i, j', Qi, j, j))
    return f

def compute_fneq(rho, j, Pi, V, W, v, Dt, tau):
    f = np.zeros(V.shape[0])
    kronecker = lambda a, b: 1 if a == b else 0
    S = compute_strain_moments(rho, j, Pi, Dt, v, tau)
    for i in range(f.size):
        Qi = compute_Qi(V[i], v)
        f[i] = - 3 * Dt * tau * W[i]/(v**2) * rho * np.einsum('ij, ij', Qi, S)
    return f

def compute_fi_moments(rho, j, Pi, V, W, v, Dt, tau):
    return compute_feq(rho, j, Pi, V, W, v, tau) + compute_fneq(rho, j, Pi, V, W, v, Dt, tau)

def compute_fout_LB(fin, tau, feq):
    return fin + 1/tau * (feq - fin)

def compute_D(msd, Dt):
    return msd/(4*Dt)

def compute_tau(visc, Dt, v):
    return visc/(Dt * 1/3 * v**2) + 1/2

def compute_MSD(P1x, P1y, P2x, P2y):
    return np.mean(np.sqrt((P1x - P2x)**2 + (P1y - P2y)**2))

def compute_mean_kinetic_energy(Vx, Vy, Masses):
    return 0.5 * np.mean(Masses * (Vx**2 + Vy**2))

def compute_kinematic_visc_auto_corr(domainSizeX, domainSizeY, Px1, Py1, Vx1, Vy1, Px2, Py2, Vx2, Vy2, Masses, rho, Dt):
    return  3/(4 * rho * (domainSizeX * domainSizeY) * compute_mean_kinetic_energy(Vx2, Vy2, Masses) * Px1.shape[0] * Dt) * (np.sum(Px2 * Masses * Vy2 - Px1 * Masses * Vy1)**2)
    

def compute_rho(Mass_list, domainSizeX, domainSizeY):
    return np.sum(Mass_list) / (domainSizeX * domainSizeY)

# def compute_Pi(j, u):
#     Pi = np.zeros((2, 2))
#     for a in range(Pi.shape[0]):
#         for b in range(Pi.shape[1]):
#             Pi[a, b] = j[a]*u[b]
#     return Pi

def compute_Pi(domainSizeX, domainSizeY, Ux, Uy):
    U = np.array([Ux, Uy])
    Pi = np.zeros((2, 2))
    for a in range(Pi.shape[0]):
        for b in range(Pi.shape[1]):
            for i in range(Ux.shape[0]):
                Pi[a, b] += U[a, i]*U[b, i]
    return Pi / (domainSizeX * domainSizeY)

def compute_Pi_eq(rho, v, u):
    Pi_eq = np.zeros((2, 2))
    kronecker = lambda a, b: 1 if a == b else 0
    for a in range(Pi_eq.shape[0]):
        for b in range(Pi_eq.shape[1]):
            Pi_eq[a, b] = rho * 1/3 * v**2 * kronecker(a, b) + rho * u[a] * u[b]
    return Pi_eq
