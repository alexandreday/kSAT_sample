from xor import XOR_SOLVE
import sys, os
import time
import pickle
import numpy as np

param = {}
### default parameters
param['N'] = 1000
param['alpha'] = 0.8
param['n_sample'] = 100
param['K'] = 3

### command line specified parameters
if len(sys.argv) > 1: # > > > read the parameters from the command line // overwrites default parameters
    for a in sys.argv[1:]:
        name, value = a.split('=')
        param[name] = float(value) # more stable across platforms to cast str to float 

param['M'] = param['alpha']*param['N']

print('[xor.py]    Running with the following parameters ...')
print(param)

start_time = time.time()
xor = XOR_SOLVE(int(param['N']), int(param['M']), int(param['K']))

X = xor.sample_solution(int(param['n_sample']), verbose=True)

print('Elapsed time:', time.time() - start_time)

# saving the solution and the formula (save the seed used to generate the formula ?) --> or not !
# want to save the solution a data file, but is dependent on the formula used ...
N=int(param['N'])
alpha=param['alpha']
K=int(param['K'])

root = 'data/'
file_formula = root+'formula/formula_idx=%i_N=%i_a=%.3f_K=%i.pkl'
file_solution = root+'sol/xor_sol_idx=%i_N=%i_a=%.3f_K=%i.pkl'

idx = 0
while os.path.isfile(file_formula%(idx, N, alpha, K)):
    idx+=1

file_formula=file_formula%(idx, N, alpha, K)
file_solution=file_solution%(idx, N, alpha, K)

print('[xor.py]   data saved to :')
print(file_formula)
print(file_solution)

pickle.dump([xor.f, xor.y], open(file_formula,'wb'))
pickle.dump(np.packbits(X), open(file_solution,'wb'))





















