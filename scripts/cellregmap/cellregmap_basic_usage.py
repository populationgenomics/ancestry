import numpy as np
from numpy import ones
from numpy.linalg import cholesky
from numpy.random import RandomState
from numpy_sugar import ddot
from numpy_sugar.linalg import economic_svd

from cellregmap import CellRegMap
from cellregmap import run_association, run_interaction, estimate_betas, run_association0

def get_L_values(hK, E):
    """
    As the definition of Ls is not particulatly intuitive,
    function to extract list of L values given kinship K and 
    cellular environments E
    """
    # # decompose K into hK: hK @ hK.T = K
    # hK = cholesky(K)

    # get eigendecomposition of EEt
    [U, S, _] = economic_svd(E)
    us = U * S

    # get decomposition of K \odot EEt
    Ls = [ddot(us[:,i], hK) for i in range(us.shape[1])]
    return Ls

random = RandomState(0)
n = 30                               # number of samples (cells)
p = 5                                # number of individuals
k = 4                                # number of contexts
y = random.randn(n, 1)               # outcome vector (expression phenotype)
C = random.randn(n, k)               # context matrix  
W = ones((n, 1))                     # intercept (covariate matrix)
hK = random.rand(n, p)               # decomposition of kinship matrix (K = hK @ hK.T)
g = 1.0 * (random.rand(n, 2) < 0.2)  # SNP vector

##################
### Old approach 

# fit null model (association test)
crm0 = CellRegMap(y, W, C, hK=hK) # hK -> background is K + EEt

# Association test
pv0 = crm0.scan_association(g)[0]
print(f'Association test p-value: {pv0}')

# get eigendecomposition of CCt
[U, S, _] = economic_svd(C)
us = U * S

# get decomposition of K \odot CCt
Ls = [ddot(us[:,i], hK) for i in range(us.shape[1])]

# fit null model (interaction test)
crm = CellRegMap(y, W, C, Ls) # Ls -> background is K*EEt + EEt

# Interaction test
pv = crm.scan_interaction(g)[0]
print(f'Interaction test p-value: {pv}')

# K = np.dot(hK,hK.T)
# K = K + np.diag(np.repeat(1e-5, n)) 

### test  get L values function

# Ls = get_L_values(K,C)
Ls = get_L_values(hK,C)

# as before (fit null, then test)
crm = CellRegMap(y, W, C, Ls)
pv = crm.scan_interaction(g)[0]
print(f'Interaction test p-value after Ls: {pv}')

###################
### Final approach 

# (each time, check that p-values are consistent)

# "get L values" function is incorporated
# null fitting and testing are merged into a "run" function

# which should take (expanded) hK 
pv0 = run_association(y, W, C, g, hK=hK)[0]
print(f'Association test p-value try: {pv0}')

# which should take (expanded) hK 
pv01 = run_association0(y, W, C, g, hK=hK)[0]
print(f'Association0 test p-value try: {pv01}')

pv = run_interaction(y, W, C, g, hK=hK)[0]
print(f'Interaction test p-value try: {pv}')

# Note that K should not be accepted as an input anymore
# pv0 = run_association(y, W, C, g, K=K)[0]
# print(f'Association test p-value try2: {pv0}')

# pv = run_interaction(y, W, C, g, K=K)[0]
# print(f'Interaction test p-value try2: {pv}')

# also added faster estimate betas function
betas = estimate_betas(y, W, C, g, hK=hK)
betaG = betas[0]
betaGxC = betas[1][0]
print(f'betaG: {betaG}')
print(f'betaGxC: {betaGxC}')
