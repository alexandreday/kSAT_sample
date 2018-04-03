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

alpha = param['M']/param['N']

new_dir = 'alpha=%.3f'%alpha
cmd = 'mkdir %s'%new_dir
cmd_rm = 'rm -rf %s'%new_dir
os.system(cmd_rm) # create new directory
os.system(cmd) # create new directory
os.system('cp merge sp verify walksat kSAT.py checks_sol.py %s/'%new_dir)
os.system('~/.conda/envs/py35/bin/python %s/kSAT.py %s %s %s'%(new_dir, argv[0], argv[1] argv[2]))
