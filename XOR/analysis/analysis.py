import numpy as np
import pickle
from matplotlib import pyplot as plt
import sys
sys.path.append('..')
from xor import XOR_SOLVE

def main():

    root_sol = '../data/sol/'
    root_formula = '../data/formula/'
    formula, y = pickle.load(open(root_formula + 'formula_idx=3_N=100_a=0.500_K=4.pkl','rb'))
    sol = np.unpackbits(pickle.load(open(root_sol + 'xor_sol_idx=3_N=100_a=0.500_K=4.pkl','rb'))).astype(int).reshape(-1, 100)
    print(formula)
    xor = XOR_SOLVE(N=100, M=50, K=4, f=formula, y=y, save_formula=False)
    print(sol.shape)
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