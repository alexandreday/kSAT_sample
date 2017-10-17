import numpy as np
import pickle, os, sys
import time


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
        seed = np.random.randint(0, sys.maxsize)
        solutions = []
        zeros = np.zeros(M,dtype=int).reshape(-1,1)

        ## reading formula file or generating new random formula 
        if read_file is None:
            formula_file = "formula.tmp_N=%i_M=%i_alpha=%.2f_K=%i.cnf"%(N, M, alpha, K) # meaningful file name !
            self.generate_formula(savefile=formula_file)
            formula = self.formula
        else:
            formula_file = read_file
            formula = np.loadtxt(formula_file, dtype=int, skiprows=1, delimiter=' ')[:,:3]
        
        if n_sample == 1:
            os.system("./sp -l %s -s%i"%(formula_file,seed))

        else:

            for _ in range(n_sample):
                print("----------------> sample # \t", _)
                #For generating SAT solutions sampled uniformly at random !

                idx_ori = np.arange(1, N+1, dtype=int)
                idx_new = np.arange(1, N+1, dtype=int)
                np.random.shuffle(idx_new)
                
                rand_permutation_map = dict(zip(idx_ori,idx_new))
                inv_rand_permutation_map = dict(zip(idx_new,idx_ori))

                isometry_formula = np.array([np.sign(x)*rand_permutation_map[abs(x)] for x in formula.flatten()], dtype=int)
                isometry_formula=isometry_formula.reshape(-1, K)
            
                file_tmp = '.tmp.cnf.formula.permutation'
                np.savetxt('.tmp.cnf.formula.permutation', np.hstack((isometry_formula, zeros)), fmt='%i', delimiter=" ", header = 'p cnf %i %i'%(N,M), comments='')

                os.system("./sp -l %s -s%i > out.txt"%(file_tmp, seed))

                if self.check_solution(solution_file='solution.tmp.lst', formula_file='.tmp.cnf.formula.permutation'):
                    sol_tmp = np.loadtxt('solution.tmp.lst', dtype=int)
                    sol_tmp_2 = np.array([np.sign(v)*inv_rand_permutation_map[abs(v)] for v in sol_tmp], dtype=int)
                    sol_tmp_2 = sol_tmp_2[np.argsort(np.abs(sol_tmp_2))]
                    
                    is_solution = self.check_solution(solution_array=sol_tmp_2)
                    if is_solution:
                        solutions.append(sol_tmp_2)

            solutions = np.sign(np.vstack(solutions))
            np.savetxt('sol_N%i=_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K), solutions,fmt="%i")

            #print(np.vstack(solutions)[:,:10])

    def check_solution(self, solution_file=None, solution_array=None, formula_file = None):

        K = self.K
        if formula_file is None:
            formula = self.formula
        else:
            formula = np.loadtxt(formula_file, dtype=int, skiprows=1, delimiter=' ')[:,:3]
    
        if solution_array is not None:
            solution = solution_array
        else:
            solution = np.loadtxt(solution_file, dtype=int)
            
        variable_label = np.abs(solution)
        var_max = np.max(variable_label)
        solution_map = np.zeros(var_max+1,dtype=int)
        solution_map[variable_label] = np.sign(solution)
        
        abs_formula = np.abs(formula)
        sign_formula = np.sign(formula)
        res = solution_map[abs_formula]*sign_formula

        return np.count_nonzero(np.sum(res, axis=1) == -K) == 0 # check that all clauses are SAT

def main():

    argv = sys.argv[1:]
    type_arg = {'n':int,'alpha':float}
    tmp = {a.split('=')[0]:a.split('=')[1] for a in argv}
    for k,v in tmp.items():
        tmp[k] = type_arg[k](float(v))

    model = KSAT(N_ = tmp['n'], alpha_ = tmp['alpha'], K_ = 3) # kSAT class
    model.generate_formula(savefile="formula.tmp.cnf") # generate random formula (optional)
    model.solve_formula(read_file="formula.tmp.cnf", n_sample= 20000) # read formula written in read_file and runs sp code

    #print(model.check_solution(solution_file='solution.tmp.lst',formula_file='.tmp.cnf.formula.permutation'))

    #t = time.time()
    #res = model.check_solution(solution_file='solution.tmp.lst',formula_file='formula.tmp.cnf')
    #print(time.time()-t)
    #print(res)






if __name__ == "__main__":
    main()