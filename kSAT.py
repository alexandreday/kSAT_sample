import numpy as np
import pickle, os, sys


class KSAT:

    def __init__(self, N_ = 1000, alpha_ = 3.8, K_=3, random_state = 0):
        np.random.seed(random_state) # by default u always get the same thing for the same parameters !
        self.N = N_
        self.alpha = alpha_
        self.K = K_
        self.M = round(self.N * self.alpha)

    def generate_formula(self, savefile = None):
        M = self.M
        alpha = self.alpha
        N = self.N
        K = self.K
        all_idx = np.arange(1, N+1,dtype=int)
        signs = np.array([-1,1],dtype=int)

        clause_set = set([])
        while len(clause_set) < M:
            literals = np.random.choice(all_idx, size=K, replace=False)
            clause = literals*np.random.choice(signs, size=K)
            clause_set.add(tuple(clause))

        zeros = np.zeros(M,dtype=int).reshape(-1,1)
        self.formula = np.array(list(clause_set))

        if savefile is not None:
            np.savetxt(savefile, np.hstack((self.formula, zeros)), fmt='%i', delimiter=" ", header = 'p cnf %i %i'%(N,M), comments='')

        return self

def main():

    model = KSAT()
    model.generate_formula(savefile="formula.tmp.cnf")

    model.solve_formula(read_file="formula.tmp.cnf", n_sample=100)




if __name__ == "__main__":
    main()