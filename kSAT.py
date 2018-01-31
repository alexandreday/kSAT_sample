import numpy as np
import pickle, os, sys
import time


def main():
    """
    Use in the following way (example):
        python kSAT.py n=1000 alpha=3.5 n_sample=10000

    What the code below does:

    -> Generates a random K-SAT (pure K-SAT) formula, which is save in "formula.tmp_N=%_M=%i_alpha=%.2f_K=%i.cnf" file
    -> Tries to find solutions to that formula. the number of solutions is specified by the n_sample parameter
    -> Solutions are supposed to be sampled uniformly at random in the space of solutions

    """
    argv = sys.argv[1:]
    type_arg = {'n':int,'alpha':float,'n_sample':int}
    tmp = {a.split('=')[0]:a.split('=')[1] for a in argv}
    for k,v in tmp.items():
        tmp[k] = type_arg[k](float(v))

    assert len(tmp) == 3, "need to specify the 3 parameters, example : n=1000 alpha=3.5 n_sample=10000"
    
    model = KSAT(N_ = tmp['n'], alpha_ = tmp['alpha'], K_ = 3, random_state=0) # kSAT class
    N, M, alpha, K = model.get_param()
    formula_file = "formula.tmp_N=%i_M=%i_alpha=%.2f_K=%i.cnf"%(N, M, alpha, K)
    
    model.generate_formula(savefile=formula_file) # generate random formula (optional)
    model.solve_formula(read_file=formula_file, n_sample=tmp['n_sample']) # read formula written in read_file and runs sp code

################################
################################
################################

def save(obj, file):
    f = open(file,'wb')
    pickle.dump(obj,f)
    f.close()

def load(file):
    f = open(file,'rb')
    return pickle.load(f)

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
        """ Generates a random formula based on N,M,alpha,K parameters specified in the constructor
        
        Parameters
        -------------
        savefile: str, optional
            file in which the CNF formula is to be saved. This is done via np.savetxt()
        
        Returns
        -------
        self
        """

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
        restart_count = 0
        n_restart = 50

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
            nn =0 
            while nn < n_sample:
                os.system("rm noconvergence.tmp.cnf")

                if restart_count > n_restart :
                    "X permutations, still not working !"
                    print("Stopping, no solutions found !")
                    break

                print("----------------> sample # \t", nn)
                #print(restart_count)
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

                seed = np.random.randint(0,2**32-1)
                os.system("./sp -l %s -s%i > out.txt"%(file_tmp, seed)) # solves the permuted formula (equivalent !)

                if os.path.exists("noconvergence.tmp.cnf"):
                    "Solution not find, try a different permutation"
                    restart_count +=1
                else:
                    nn+=1
                    restart_count = 0

                    if self.check_solution(solution_file='solution.tmp.lst', formula_file='.tmp.cnf.formula.permutation'):
                        sol_tmp = np.loadtxt('solution.tmp.lst', dtype=int)
                        sol_tmp_2 = np.array([np.sign(v)*inv_rand_permutation_map[abs(v)] for v in sol_tmp], dtype=int)
                        sol_tmp_2 = sol_tmp_2[np.argsort(np.abs(sol_tmp_2))]
                        
                        is_solution = self.check_solution(solution_array=sol_tmp_2)
                        if is_solution:
                            solutions.append(sol_tmp_2)
                    if nn % (n_sample // 10) == 0 and n_sample > 10 and len(solutions) > 0:
                        print(nn, " saving")
                        solution_stack = np.sign(np.vstack(solutions))
                        solution_stack[solution_stack < 0] = 0
                        save(np.packbits(solution_stack, axis=1),'sol_N=%i_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K))
                        #np.savetxt('sol_N=%i_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K), np.packbits(solution_stack, axis=1), fmt="%i")
            
            if len(solutions) > 0:
                solution_stack = np.sign(np.vstack(solutions))
                solution_stack[solution_stack < 0] = 0
                save(np.packbits(solution_stack, axis=1), 'sol_N=%i_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K))
                #np.savetxt('sol_N=%i_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K), np.packbits(solution_stack, axis=1), fmt="%i")
                #np.savetxt('sol_N=%i_M=%i_alpha=%.2f_K=%i.txt'%(N,M,alpha,K), solution_stack,fmt="%i")
            #print(np.vstack(solutions)[:,:10])

    def check_all_solution(self, N, solution_file, formula_file, hist=True):
        formula = np.loadtxt(formula_file, dtype=int, skiprows=1, delimiter=' ')[:,:3]
        self.formula = formula
        tmp = np.unpackbits(load(solution_file),axis=1)[:,N]
        all_solution = tmp.astype(int)
        N_sol = all_solution.shape[0]


        sol_result=[]
        idx_var = np.arange(1,N+1,dtype=int)

        count_true = 0
        print("Checking %i solutions"%N_sol)
        for i in range(N_sol):
            sol = all_solution[i]*idx_var
            res = self.check_solution(solution_array=sol)
            if res is True:
                count_true +=1
            else:
                print("Found wrong solution !?")
        print("%i out of the %i solutions are correct"%(count_true,N_sol))
        n_unique = len(np.unique(all_solution, axis=0))
        print("number of unique solutions :\t %i"%(n_unique))
    
        if hist is True:
            from matplotlib import pyplot as plt
            import seaborn as sns
            mag = np.mean(all_solution,axis=0)
            nbin = max([N_sol/10,20])
            N, M, alpha, K =  self.infer_parameters(solution_file)
            sns.distplot(mag,bins=nbin,kde=False)
            plt.xlabel('magnetization')
            plt.title('nsample = %i, N=%i, M=%i, alpha=%.2f, K=%i'%(N_sol, N, M, alpha, K))
            plt.show()

    def infer_parameters(self, solution_file):
        sol_s = solution_file.strip(".txt").split('_')
        for p in sol_s:
            if '=' in p:
                ps = p.split('=')
                if ps[0] == 'N':
                    N = int(float(ps[1]))
                elif ps[0] == 'M':
                    M = int(float(ps[1]))
                elif ps[0] == 'alpha':
                    alpha = float(ps[1])
                elif ps[0] == 'K':
                    K = int(float(ps[1]))
        
        return N, M, alpha, K
    
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

if __name__ == "__main__":
    main()