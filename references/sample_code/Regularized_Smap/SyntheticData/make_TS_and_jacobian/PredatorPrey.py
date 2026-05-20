import numpy as np
import sdeint
import pylab as plt
from scipy import integrate


#### System size
SS = 10000
#### Number of dimensions
d = 2
b = np.array([1,-0.5])
##### 
A = np.asmatrix(np.array([[0, .5], [-.3, 0.]]))
T = 800.;
dt = 0.01;
n_steps = T/dt;
t = np.linspace(0, T, n_steps)
#### Create a Wiener process so that you can use it for the computation of the jacobian
d_Wiener = np.random.normal(0., np.sqrt(dt), (len(t) - 1, 2))


def StampaSerieTemporale(X, nome, tau):
    idex = 0.; tdex = 2.
    n = open(nome, 'w')
    n.write('%f %f\n' % (X[0, 0], X[0, 1]))
    for j in range(1, np.shape(X)[0]):
        idex += 1
        if((idex*dt) % tau == 0):
            n.write('%f %f\n' % (X[j, 0], X[j, 1]))
            tdex += 1
    n.close()
def f(X, t):
    dydt = np.array([X[s]*(b[s] - np.sum(np.dot(A,X)[0,s])) for s in range(0,len(X))])
    return dydt
def u(X, t):
    dydt = np.array([1./np.sqrt(SS) * np.sqrt(X[s]*(abs(b[s]) + np.sum(np.dot(abs(A),X)[0,s]))) for s in range(0,len(X))])
    return np.diag(dydt)


x0 = np.array([3.5,3])
ts = sdeint.itoSRI2(f, u, x0, t,dW = d_Wiener)
tau = 1
StampaSerieTemporale(ts, 'PredatorPrey.txt', tau)
StampaSerieTemporale(d_Wiener, 'rumore_pred_prey.txt', tau)

plt.plot(ts[:,0], ts[:,1], c = 'red', lw = 3)
plt.show()
