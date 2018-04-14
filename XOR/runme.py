from xor import XOR_SOLVE
import sys, os
import time

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























