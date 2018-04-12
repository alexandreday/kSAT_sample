import numpy as np
import time

def verify_solution(A, y, sol):
    nclause = A.shape[0]
    for i in range(nclause):
        if np.sum(A[i, :] * sol) % 2 != y[i]:
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

def sample_solution_GE(A, y, pivot_ls, n_sample=10):
    """ WALKSAT on the set of clauses that do not contain pivots followed by fixing clauses
    with pivots. Note that A has to begin in a reduced form.
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
    
    return xsol

def main():
    import time
    from xor import generate_sparse

    N = 1000 # number of literals / variables
    M = 910 # number of constraints
    K = 3
    n_sample = 100
    A, f = generate_sparse(N=N, M=M, K=K)
    y = np.random.randint(0, 2, M) # random constraints
    pivots = make_diagonal(A, y)

    

if __name__ == "__main__":
    main()