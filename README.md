Python wrapper for survey propagation/belief propagation C-code [see R. Zecchina].
To read solutions:

import pickle
import numpy as np
N=100
M=350
alpha=M/N
fname = 'sol_N=%i_M=%i_alpha=%.2f_K=3.pkl'%(N, M, alpha)
a=pickle.load(open(fname,'rb'))

######### DATA ###########
X = np.unpackbits(a).astype(int).reshape(-1,N)
s = np.unpackbits(a).astype(int).reshape(-1,N).shape