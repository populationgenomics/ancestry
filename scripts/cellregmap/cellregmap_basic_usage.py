import numpy as np
from numpy import ones, concatenate
from numpy.random import RandomState
from numpy_sugar import ddot
from numpy_sugar.linalg import economic_svd

from cellregmap import CellRegMap

random = RandomState(1)
n = 30                               # number of samples (cells)
p = 5                                # number of individuals
k = 4                                # number of contexts
y = random.randn(n, 1)               # outcome vector (expression phenotype)
C = random.randn(n, k)               # context matrix  
W = ones((n, 1))                     # intercept (covariate matrix)
hK = random.randn(n, p)              # decomposition of kinship matrix (K = hK @ hK.T)
g = 1.0 * (random.rand(n, 1) < 0.2)  # SNP vector

W = concatenate([W, g], axis=1)

# get eigendecomposition of CCt
[U, S, _] = economic_svd(C)
us = U * S

# get decomposition of K \odot CCt
Ls = [ddot(us[:,i], hK) for i in range(us.shape[1])]

# fit null model (interaction test)
crm = CellRegMap(y, W, C, Ls)
# Interaction test
pv = crm.scan_interaction(g)[0]
print(pv)