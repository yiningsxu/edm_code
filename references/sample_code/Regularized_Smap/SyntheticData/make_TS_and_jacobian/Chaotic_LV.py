import pylab as plt
from scipy import integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import sdeint
import scipy
import os, sys
from mpl_toolkits.mplot3d import Axes3D
############################################################# Chaotic Lotka-Volterra ####################################
r = np.array([1., 0.72, 1.53, 1.27])
A = np.matrix([[1.*r[0], 1.09*r[0], 1.52*r[0], 0.*r[0]],  [0.*r[1], 1.*r[1], 0.44*r[1],1.36*r[1]],  [2.33*r[2], 0.*r[2], 1.*r[2], 0.47*r[2]], [1.21*r[3], 0.51*r[3], 0.35*r[3], 1.*r[3]]])

def StampaSerieTemporale(X, nome, tau):
    idex = 0.; tdex = 2.
    n = open(nome, 'w')
    n.write('%f %f %f %f\n' % (X[0, 0], X[0, 1], X[0, 2], X[0, 3]))
    for j in range(1, np.shape(X)[0]):
        idex += 1
        if((idex*dt) % tau == 0):
            n.write('%f %f %f %f\n' % (X[j, 0], X[j, 1], X[j, 2], X[j, 3]))
            tdex += 1
    n.close()
################################################
###### Initial conditions and time steps #######
x0 = 0.2; y0 = 0.2; z0 = 0.3; k0 = 0.3;
####### 
T = 600;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
X_f1 = np.array([x0, y0, z0, k0])
d_Wiener = np.random.normal(0., np.sqrt(dt), (len(t) - 1, 4))

################################################
SystemSize =  10000

######## Auxiliar Functions ####################
def HeavisideTheta(x):
    if x > 0:
        return(1)
    else:
        return(0)
################################################
def f(X, t = 0):
    dydt = np.array([X[s]*(r[s] - np.sum(np.dot(A,X)[0,s]))for s in range(0,len(X))])
    return(dydt)

    return(dydt)
### Noise with BC
def u(X, t):
    dydt = np.array([1./np.sqrt(SystemSize) * np.sqrt(X[s]*(r[s] + np.sum(np.dot(A,X)[0,s])) * HeavisideTheta(X[s])) for s in range(0,len(X))])
    return np.diag(dydt)
################
printing = True
ini_cond = integrate.odeint(f, X_f1, t)
X_f1 = np.array([ini_cond[len(t)-1,0], ini_cond[len(t)-1,1], ini_cond[len(t)-1,2], ini_cond[len(t)-1,3]])

######### Now that you are on the attractor run the true simulation
stochastic_dynamics = sdeint.itoSRI2(f, u, X_f1, t, dW = d_Wiener)
if printing == True:
    tau = 2.
    StampaSerieTemporale(stochastic_dynamics, 'Chaotic_LV.txt', tau)
    StampaSerieTemporale(d_Wiener, 'rumore_chaotic.txt', tau)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot(stochastic_dynamics[:,0], stochastic_dynamics[:,2], stochastic_dynamics[:,3], color = 'b')

plt.show()





