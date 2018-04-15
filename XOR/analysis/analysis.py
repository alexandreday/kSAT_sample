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

    formula, y = pickle.load(open(root_formula + 'formula_idx=0_N=1000_a=0.800_K=3.pkl','rb'))
    X = np.unpackbits(pickle.load(open(root_sol + 'xor_sol_idx=0_N=1000_a=0.800_K=3.pkl','rb'))).astype(int).reshape(-1, N)
    xor = XOR_SOLVE(N=1000, M=800, K=3, f=formula, y=y, save_formula=False)



if __name__ == "__main__":
    main()