import numpy as np
import numdifftools as nd

### open the noise
f = open('rumore_pred_prey.txt', 'r')
rumore = np.asmatrix([map(float,line.split(' ')) for line in f ])
rumore = rumore.astype('float32')
### open the time series
g = open('PredatorPrey.txt', 'r')
time_series = np.asmatrix([map(float,line.split(' ')) for line in g ])
time_series = time_series.astype('float32')
### Recall the Wiener process used to generate the noise
d_Wiener = rumore
### Matrix and growth rates
r = np.array([1,-0.5])
A = np.asmatrix(np.array([[0, .5], [-.3, 0.]]))
### System size
SystemSize = 1000
def HeavisideTheta(x):
    if x > 0:
        return(1)
    else:
        return(0)
################################################
### Stochastic Vector field
def dX_dt(X):
    dydt = np.array([X[s]*(r[s] - np.sum(np.dot(A,X)[0,s])) + 1./np.sqrt(SystemSize) * d_Wiener[s,:] * np.sqrt(X[s]*(abs(r[s]) + np.sum(np.dot(abs(A),X)[0,s]))) * HeavisideTheta(X[s]) for s in range(0,len(X))])
    return(dydt)

### Compute the analytical Jacobian
num_species = time_series.shape[1]
jacobian_matrix = open('jacobian_predator.txt', 'w')
for k in range(0,time_series.shape[0]):
	point = np.squeeze(np.asarray(time_series[k,:]))
	f_jacob = nd.Jacobian(dX_dt)(point)
	for u in range(0,num_species):
		for z in range(0,num_species):
		    	jacobian_matrix.write('%lf ' % (f_jacob[u,z,0]))
	jacobian_matrix.write('\n')
jacobian_matrix.close()
