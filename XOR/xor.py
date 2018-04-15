import numpy as np
import time
import pickle
from collections import Counter

def main():

    alpha = np.arange(0.5, 1.001, 0.01)
    N=100
    K=3
    time_vs_alpha = []
    entropy_vs_alpha= []


    for a in alpha:
        print('alpha = %.3f'%a)
        M = int(a*N)
        xor = XOR_SOLVE(N, M, K, save_formula = True)
        t_start= time.time()
        X = xor.sample_solution(10000, verbose=True)
        time_vs_alpha.append(time.time() - t_start)
        entropy_vs_alpha.append(xor.entropy())

        for x in X:
            assert xor.check_solution(x)
        save_solution_file = 'sol/solution_N=%i_M=%i_K=%i.pkl'%(N,M,K)
        if len(X) > 0:
            pickle.dump(np.packbits(X), open(save_solution_file,'wb'))

    pickle.dump(entropy_vs_alpha, open('entropy/entropy_N=%i_M=%i_K=%i.pkl'%(N,M,K),'wb'))
    pickle.dump(time_vs_alpha, open('time/time_N=%i_M=%i_K=%i.pkl'%(N,M,K),'wb'))

def sample_tuple(N, K): # sample tuple uniformly at random
    trace = []
    tup = []
    for i in range(K):
        pos = np.random.randint(0, N-i)
        value = pos
        for t in reversed(trace): # iterate from last element
            if value > t-1:
                value +=1
        tup.append(value)
        trace.append(pos)#print(trace)
    tup.sort()
    return tup

def generate_XOR_formula(N=10, M=10, K=3):
    """ generates a XORSAT formula
    """
    formula = set([])

    while len(formula) < M:
        tup = tuple(sample_tuple(N, K))
        formula.add(tup)
    return list(formula) # better to work with lists => order is always preserved ! 

def generate_sparse(N=10, M=10, K=3, formula=None):
    # Can specify formula, but still have to specify N,M,K

    if formula is not None:
        formula = formula
    else:
        formula = generate_XOR_formula(N,M,K)

    A = np.zeros((M,N),dtype=int)
    for i, clause in enumerate(formula):
        for literal in clause:
            A[i, literal]=1
    return A, formula

def verify_solution(A, y, sol):
    nclause = A.shape[0]
    for i in range(nclause):
        if np.dot(A[i, :], sol) % 2 != y[i]:
            return False
    return True

def swap_rows(A, ri, rj):
    tmp = np.copy(A[ri, :]) # careful when slicing, will return a view
    A[ri, :] = A[rj, :]
    A[rj, :] = tmp
   
def swap(y, i, j):
    tmp = y[i]
    y[i] = y[j]
    y[j] = tmp

def add_to_row(A, ri, rj): # A[ri, :] <- A[ri, :] + A[rj, :]
    A[ri, :] = A[ri, :] + A[rj, :]

def make_diagonal(A_, y_, copy=False):
    """ This reduction is unique """ 
    if copy:
        A=np.copy(A_)
        y=np.copy(y_)
    else:
        A = A_
        y = y_
    
    #1 clean up zero columns ? (no !)
    M, N = A.shape

    pos_pivot = 0
    pivot_list = []
    j = 0
    for i in range(M): # go over lines, for each line find pivot !
        for j in range(i, N):
            pos_one = np.where(A[:,j] == 1)[0]
            if len(pos_one) > 0:
                pos_one = pos_one[pos_one > (i - 1)]
            if len(pos_one) > 0:
                if A[i, j] == 0:
                    swap_rows(A, i, pos_one[0])
                    swap(y, i, pos_one[0])

                pos_one = np.where(A[:,j] == 1)[0]
                for k in pos_one:
                    if k > i :
                        A[k] += A[i] # mod 2
                        A[k] = np.remainder(A[k], 2)
                        y[k] = (y[k] + y[i])%2 # mod 2
                pivot_list.append([i, j])
                break

    for pivot in reversed(pivot_list):
        i, j = pivot
        pos_one = np.where(A[:,j] == 1)[0]
        pos_one = pos_one[pos_one < i]
        for k in pos_one:
            A[k] += A[i] # mod 2
            A[k] = np.remainder(A[k], 2)
            y[k] = (y[k] + y[i])%2 # mod 2

    if copy is True:
        return A, y, np.array(pivot_list,dtype=int)
    else:
        return np.array(pivot_list,dtype=int)

def find_pivot(A_UT):
    return np.where(np.diagonal(A_UT) == 1)[0]

def solve_ES(A, y):
    """ Solver using exhaustive search of all configurations """
    nvar = A.shape[1]
    nsol = 2**nvar
    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=nvar)), dtype=np.int)
    sol_list = []
    for i in range(nsol):
        sol = b2_array(i)
        if check_solution(A, y, sol):
            sol_list.append(sol)
    return sol_list

def marginals(solution_set): # unique variable marginals
    if len(solution_set) > 0:
        return np.mean(solution_set, axis=0)
    else:
        []

def enumerate_solution_GE(A, y):
    """ A has to be in a reduced form """
    M, N = A.shape
    
    pivots = np.where(np.diagonal(A) == 1)[0] # pivots are identified
    none_pivots = np.setdiff1d(np.arange(N), pivots)
    none_pivots_2 = np.setdiff1d(np.arange(M), pivots)

    rank = len(pivots) # matrix rank
    xsol = -1*np.ones(N, dtype=int) # unconstrained variables are marked by a -2

    N_free = N - rank # upper bound on number of log of number of solutions (but may be fewer !)
 
    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=N_free)), dtype=np.int)
    all_sol = []

    pivot_set = set(list(pivots))

    for i in range(2**N_free): # HERE REPLACE BY SAMPLING, THIS THE SPACE WE WISH TO SAMPLE !
        is_sol = False
        xsol = -1*np.ones(N, dtype=int) # unconstrained variables are marked by a -2
        xsol[none_pivots] = b2_array(i)
        y_res = np.remainder(np.dot(A[:,none_pivots], xsol[none_pivots].T) + y, 2)

        if np.count_nonzero(y_res[none_pivots_2] == 1) == 0:
            xsol[pivots] = y_res[pivots]
            is_sol = True
        if is_sol:
            all_sol.append(xsol)
        
    return all_sol

def is_SAT(A, y):
    tmp = np.sum(A, axis=1)
    return np.count_nonzero(y[tmp == 0] == 1) == 0 

def sample_solution_GE(A, y, pivot_ls):
    """ A is in row-echelon form with pivot position provided in pivot_ls
    """
    M, N = A.shape
    t_init = time.time()
    n_pivot = len(pivot_ls)
    n_free = N - n_pivot
    pos_pivot = pivot_ls[:,1]
    none_pivot_pos = np.setdiff1d(np.arange(N), pos_pivot)
    xsol = np.ones(N,dtype=int)

    xsol[none_pivot_pos] = np.random.randint(0, 2, n_free)
    for p in reversed(pivot_ls):
        i, j = p
        xsol[j]^= (y[i] + np.dot(A[i, :], xsol)) % 2
    
    assert verify_solution(A, y, xsol)
    return xsol

class XOR_SOLVE:

    def __init__(self, N=100, M=80, K=3, f=None, y=None, save_formula=True):
        self.N = N
        self.M = M
        self.K = K
        self.A, self.f = generate_sparse(N, M, K, formula = f) # A is not reduced at this point
        
        if save_formula:
            file_name = 'formula/formula_N=%i_M=%i_K=%i.pkl'%(N,M,K)
            pickle.dump(self.f, open(file_name,'wb'))

        if y is not None:
            self.y_original = np.copy(y)
        else:
            self.y_original= np.random.randint(0, 2, M) # random constraints

        self.y = np.copy(self.y_original)
        self.is_reduced = False
    
    def reduce_system(self):
        self.pivots = make_diagonal(self.A, self.y)
        self.is_reduced = True
    
    def entropy(self):
        if self.SAT():
            return (self.N - len(self.pivots))/self.N
        else:
            return 0

    def SAT(self):
        if not self.is_reduced:
            self.reduce_system()
        return is_SAT(self.A, self.y)
    
    def sample_solution(self, n_sample = 10, verbose = 0):
        
        if not self.is_reduced:
            self.reduce_system()
        if self.SAT() is False:
            print("No SOLUTION")
            return []

        x_sample_solution = []
        for i in range(n_sample):
            if verbose != 0:
                if i % 500 == 0:
                    print(i)
            x_sample_solution.append(sample_solution_GE(self.A, self.y, self.pivots))

        return x_sample_solution
    
    def check_solution(self, x):
        return verify_solution(self.A, self.y, x)

if __name__ == "__main__":
    main()