import numpy as np
from numpy import ones
from numpy.linalg import cholesky
from numpy.random import RandomState
from numpy_sugar import ddot
from numpy_sugar.linalg import economic_svd

from cellregmap import CellRegMap
from cellregmap import run_association, run_interaction

def get_L_values(K, E):
    """
    As the definition of Ls is not particulatly intuitive,
    function to extract list of L values given kinship K and 
    cellular environments E
    """
    # decompose K into hK: hK @ hK.T = K
    hK = cholesky(K)

    # get eigendecomposition of EEt
    [U, S, _] = economic_svd(E)
    us = U * S

    # get decomposition of K \odot EEt
    Ls = [ddot(us[:,i], hK) for i in range(us.shape[1])]
    return Ls

random = RandomState(1)
n = 30                               # number of samples (cells)
p = 5                                # number of individuals
k = 4                                # number of contexts
y = random.randn(n, 1)               # outcome vector (expression phenotype)
C = random.randn(n, k)               # context matrix  
W = ones((n, 1))                     # intercept (covariate matrix)
hK = random.rand(n, p)              # decomposition of kinship matrix (K = hK @ hK.T)
g = 1.0 * (random.rand(n, 1) < 0.2)  # SNP vector

# fit null model (association test)
crm0 = CellRegMap(y, W, C, hK=hK)

# Association test
pv0 = crm0.scan_association(g)[0]
print(f'Association test p-value: {pv0}')

# get eigendecomposition of CCt
[U, S, _] = economic_svd(C)
us = U * S

# get decomposition of K \odot CCt
Ls = [ddot(us[:,i], hK) for i in range(us.shape[1])]

# fit null model (interaction test)
crm = CellRegMap(y, W, C, Ls)

# Interaction test
pv = crm.scan_interaction(g)[0]
print(f'Interaction test p-value: {pv}')

K = np.dot(hK,hK.T)
K = K + np.diag(np.repeat(1e-5, n)) 

### test  get L values function

Ls = get_L_values(K,C)

crm = CellRegMap(y, W, C, Ls)
pv = crm.scan_interaction(g)[0]
print(f'Interaction test p-value after Ls: {pv}')

pv0 = run_association(y, W, C, g, hK=hK)[0]
print(f'Association test p-value try: {pv0}')

pv = run_interaction(y, W, C, g, hK=hK)[0]
print(f'Interaction test p-value try: {pv}')

pv0 = run_association(y, W, C, g, K=K)[0]
print(f'Association test p-value try2: {pv0}')

pv = run_interaction(y, W, C, g, K=K)[0]
print(f'Interaction test p-value try2: {pv}')