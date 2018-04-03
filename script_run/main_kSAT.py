import sys, os

############################ READ COMMAND LINE ##########################

argv = sys.argv
assert len(argv) == 4 # need to specify number of clauses, alpha and number of samples
argv = argv[1:]

param_type={
    'M':int,
    'N':int,
    'n_sample':int
}

param = {}

for a in argv:
    k, v = a.split('=')
    v = int(float(v))
    param[k] = v

print(param)
alpha = param['M']/param['N']

new_dir = 'alpha=%.3f'%alpha
cmd = 'mkdir %s'%new_dir
cmd_rm = 'rm -rf %s'%new_dir
os.system(cmd_rm) # create new directory
os.system(cmd) # create new directory
os.system('cp merge sp verify walksat kSAT.py check_sol.py %s/'%new_dir)
os.chdir(new_dir)
new_cmd = '~/.conda/envs/py35/bin/python kSAT.py n=%i alpha=%.3f n_sample=%i'%(param['N'], alpha, param['n_sample'])
#print(new_cmd)
os.system(new_cmd)
