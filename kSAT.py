import numpy as np
import pickle, os, sys


class KSAT:

    def __init__(self, N_ = 1000, alpha_ = 3.8, K_=3, random_state = 0):
        np.random.seed(random_state) # by default u always get the same thing for the same parameters !
        self.N = N_
        self.alpha = alpha_
        self.K = K_
        self.M = round(self.N * self.alpha)
    
    def get_param(self):
        return self.N, self.M, self.alpha, self.K

    def generate_formula(self, savefile = None):
        N, M, alpha, K = self.get_param()
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

    def solve_formula(self, read_file = None, n_sample = 1):
        
        N, M, alpha, K = self.get_param()

        ## reading formula file or generating new random formula 
        if read_file is None:
            formula_file = "formula.tmp_N=%i_M=%i_alpha=%.2f_K=%i.cnf"%(N, M, alpha, K) # meaningful file name !
            self.generate_formula(savefile=formula_file)
        else:
            formula_file=read_file
        
        



        


def main():

    model = KSAT() # kSAT class
    model.generate_formula(savefile="formula.tmp.cnf") # generate random formula (optional)
    model.solve_formula(read_file="formula.tmp.cnf", n_sample=100) # read formula written in read_file and runs sp code





if __name__ == "__main__":
    main()