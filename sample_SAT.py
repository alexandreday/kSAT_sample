import os 
import sys
import numpy as np
import pickle
from collections import Counter

def runSAT(n_var = 10000, alpha = 3.9, n_sample = 10000, file='sol_SAT.pkl'):

    seed = np.random.randint(0, 2**32-1)
    #cmd = "./sp -n%i -a%.2f -s%i -l%s > out.idc.txt"
    cmd = "./sp -l formula.tmp.cnf -s%i > out.idc.txt"

    solution = []
    interval = int(round(n_sample / 10))
    count = 0
    for i in range(n_sample):
        if (i+1) % interval == 0 or i == (n_sample - 1):
            print(i)
            pickle.dump(np.vstack(solution), open(file,'wb'))

        seed = np.random.randint(2**32-1)
        file_name = "formula.tmp.cnf"
        torun = cmd%(seed)
        print(torun)
        #torun = cmd
        #blockPrint()
        os.system(torun)
        #enablePrint()

        ## read results 
        f = open('wsat.tmp.out','r')
        ii = 0
        for l in f:
            if l == 'ASSIGNMENT FOUND\n':
                break
            ii+=1

        res = np.loadtxt('wsat.tmp.out', skiprows=ii+1, dtype=int, usecols=1)
        solution.append(np.sign(res))

    return np.vstack(solution)

def main():

    parameters = sys.argv[1:]
    n_var = int(parameters[0])
    alpha = float(parameters[1])
    n_sample = int(parameters[2])

    file = 'sol_alpha=%.2f_nvar=%i.txt'%(alpha, n_var)

    solutions = runSAT(n_var=n_var, alpha=alpha, n_sample=n_sample, file=file)

    check_data(file)

def check_data(file):
    data = pickle.load(open(file,'rb'))
    print(data)
    print(np.mean(data,axis=0))
    print(np.std(data,axis=0))
    print(np.unique(np.std(data,axis=0)))
    return data


if __name__ == "__main__":
    main()