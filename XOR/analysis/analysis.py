import numpy as np
import pickle
from matplotlib import pyplot as plt
import sys
sys.path.append('..')
from xor import XOR_SOLVE



def main():

    root_sol = '../data/sol/'
    root_formula = '../data/formula/'

    N=1000
    formula, y = pickle.load(open(root_formula + 'formula_idx=2_N=1000_a=0.800_K=3.pkl','rb'))
    sol = np.unpackbits(pickle.load(open(root_sol + 'xor_sol_idx=2_N=1000_a=0.800_K=3.pkl','rb'))).astype(int).reshape(-1, N)
    xor = XOR_SOLVE(N=1000, M=800, K=3, f=formula, y=y, save_formula=False)

    #print(xor.A.shape)
    print(np.dot(sol[0], xor.A[0])%2 == xor.y[0])
    print(len(sol[0]))
    exit()
    #exit()
    #print(np.where(xor.A[0] == 1)[0])
    for i, s in enumerate(formula):
        tmp_1 = list(s)
        tmp_2 = np.where(xor.A[i] == 1)[0]
        for j in range(3):
            assert tmp_1[j] == tmp_2[j]
        #print(list(s),'\t', np.where(xor.A[i] == 1)[0])
        #if i > 20:
        #    break
    print("DONE")
    exit()
    print((np.sum(xor.A[0]*sol[0])+y[0])%2)
    print(sol.shape)
    exit()
    for s in sol:
        print(xor.check_solution(s))
    exit()
    #for s in sol:

    
    #xor_sol_idx=3_N=100_a=0.500_K=4.pkl'

    alpha = np.arange(0.5,1.001,0.01)
    N=100
    K=3

    s = pickle.load(open('../entropy/entropy_N=100_M=100_K=3.pkl','rb'))
    t = pickle.load(open('../time/time_N=100_M=100_K=3.pkl','rb'))

    plt.scatter(alpha,t)
    plt.show()
    exit()
    s = []
    t = []
    for a in alpha:
        M=int(a*N)
        sign = 'N=%i_M=%i_K=%i.pkl'%(N,M,K)






if __name__ == "__main__":
    main()