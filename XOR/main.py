import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

def check_solution(A, y, sol):
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
    for j in range(M): # enumerate over columns
        #print(j)
        pos_one = np.where(A[:, j] == 1)[0]
        if len(pos_one) > 0:
            pos_one_swap = pos_one[(pos_one > j)] # can only swap with rows below j.
            if len(pos_one_swap) != 0:
                
                if A[j,j] == 0:
                    swap_rows(A, j, pos_one_swap[0])
                    swap(y, j, pos_one_swap[0])

                pos_one = np.where(A[:, j] == 1)[0]
                for i in pos_one:
                    if i > j:
                        A[i, :] += A[j, :] # mod 2
                        A[i, :] = np.remainder(A[i, :], 2)
                        y[i] = (y[j] + y[i])%2 # mod 2

    remove_one_UT(A, y)

    return A, y

def remove_one_UT(A, y):
    """Removes the non-zero elements in the upper-triangular part using the pivots"""
    pivot = find_pivot(A)
    for j in pivot:
        pos_one = np.where(A[:, j] == 1)[0]
        for i in pos_one:
            if i < j:
                A[i, :] += A[j, :]  # mod 2
                A[i, :] = np.remainder(A[i, :], 2)
                y[i] = (y[i] + y[j]) % 2  # mod 2

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

        
def enumerate_solution(Areduce, y_):
    """ This is exponential in the number of free variables, which is of order O(2^(M-N)) 
    """
    y = y_
    test = np.sum(Areduce, axis = 1)
    pos_no_constraint = (test == 0)
    if np.count_nonzero(y[pos_no_constraint] == 1) > 0:
        return []
    else:
        #print('removing :\t'
        A = Areduce[pos_no_constraint == False] # remove no constraint lines !
        y = y_[pos_no_constraint == False]

    M, N = A.shape
    N_free = N-M
    N_check = 2**N_free
    
    all_sol = []
    N_check = 2**(N_free)
    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=N_free)), dtype=np.int)
    
    for i in range(2**N_free): # does not always work, if a column is full of zeros
        xsol = np.zeros(N, dtype=int)
    
        xsol[-N_free:] = b2_array(i) # error here.
        is_sol = True
        for i in reversed(range(M)):
            if np.sum(A[i,:]*xsol) % 2 == y[i]:
                bit = 0
            else:
                bit = 1

            if A[i,i] == 1:
                xsol[i] = bit
            else:
                if bit == 1: # no way to satisfy with any value here
                    is_sol = False
                else: # can take any value !
                    xsol[i] = -1 # star variable !

        if is_sol:
            pos_star = np.where(xsol == -1)[0]
            n_star = len(pos_star)# counts the number of unfrozen variables
            #print('number of star variable\t:%i'%n_star)
            b2_array_spec = lambda n10 : np.array(list(np.binary_repr(n10, width=n_star)), dtype=np.int)
            if n_star > 0:
                for iii in range(2**n_star):
                    tmp = b2_array_spec(iii)
                    for jjj,v in enumerate(tmp):
                        xsol[pos_star[jjj]] = tmp[jjj]
                    all_sol.append(np.copy(xsol))
            else:
                all_sol.append(np.copy(xsol))
    
    return all_sol

def swap_back(sol, swap_history):
    nswap = len(swap_history)
    for s in sol:
        for i in range(nswap):
            p1, p2 = swap_history[-(i+1)]
            swap(s, p1, p2)

def solve_GE(A, y):
    """ Solver using Gaussian elimination and enumerating free variables """
    Anew, ynew, swap_history = make_UT(A, y)
    Afinal, yfinal = make_diag(Anew, ynew)
    print(Afinal)
    return enumerate_solution(Afinal, yfinal), swap_history

def marginals(solution_set): # unique variable marginals
    if len(solution_set) > 0:
        return np.mean(solution_set, axis=0)
    else:
        []


def enumerate_solution_v2(A, y):
    """ A has to be in a reduced form """
    M, N = A.shape
    
    pivots = np.where(np.diagonal(A) == 1)[0] # pivots are identified
    rank = len(pivots) # matrix rank
    max_pivot = pivots[-1] 
    xsol = -1*np.ones(N, dtype=int) # unconstrained variables are marked by a -2

    N_free = N - max_pivot -1

    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=N_free)), dtype=np.int)
    all_sol = []

    for i in range(2**N_free): # here this could simply be replaced by some uniform sampling

        xsol = -1*np.ones(N, dtype=int) # unconstrained variables are marked by a -2
        xsol[-N_free:] = b2_array(i)
        y_res = np.remainder(np.dot(A[:,-N_free:],xsol[-N_free:].T) + y, 2)
    
        is_sol = True
        for j in reversed(range(M)):
            if j > max_pivot:
                if y_res[j] == 1:
                    is_sol = False
                    break # not a solution
            else:
                if A[j,j] == 1:
                    xsol[j] = y_res[j]
                else:
                    if y_res[j] == 1:
                        is_sol = False
                        break # not a solution
        if is_sol:
            all_sol.append(xsol)
        
    return all_sol


def main():
    from xor import generate_sparse
    N = 8
    M = 6
    K = 3
    n_col = N
    n_row = M
    np.random.seed(30)
    for i in range(1000):
        print(i)
        A, f = generate_sparse(N=N, M=M, K=K)
        y = np.random.randint(0, 2, n_row)
        sol_init = solve_ES(A, y)

        make_diagonal(A, y)
        sol_final = enumerate_solution_v2(A, y)

        if len(sol_final)>0:
            n_star = np.count_nonzero(sol_final[0] == -1)

        assert len(sol_init) == len(sol_final)*2**n_star

        if i == 9:
            print(A)
            print(y)
            print(sol_init)
            print(sol_final)
            print(marginals(sol_final))
            print(marginals(sol_init))
            exit()

        if len(sol_init) > 0:
            marg_final = marginals(sol_final)
            marg_final[marg_final < -0.0001] = 0.5 # star variables
            if abs(np.linalg.norm(marginals(sol_init) - marg_final)) > 1e-5:
                print(i)
                assert False
        #print(marginals(sol_init))
        #print(marginals(sol_final))

    exit()

    for i in range(200):
    
        print(i)
        A, f = generate_sparse(N=N, M=M, K=K)
        y = np.random.randint(0, 2, n_row)

        if i == 16:
            print(A)
            print(y)
            #exit()

            sol_ES = solve_ES(A, y)
            sol_GE, swap_hist = solve_GE(A, y)
        #if i==16:
            print(sol_ES)
            print(sol_GE)
            exit()

            swap_back(sol_GE, swap_hist)
            diff = np.linalg.norm(marginals(sol_GE)-marginals(sol_ES))
            if abs(diff) > 1e-10:
                print(i)

    exit()

    print(A)
    print('col sum:\t',np.sum(A,axis=0))
    print(y)
    #exit()
    
    #print(A)
    #print(y)   
    print("Number of solutions original:", len(sol_list))
    #exit()
    Anew, ynew, swap_history = make_UT(A, y)
    print(Anew)
    print(ynew)
    Afinal, yfinal = make_diag(Anew, ynew)

    print(Afinal)
    print(yfinal)
    sol_list_new = enumerate_solution(Afinal, yfinal)
    #print(len(sol_list_new))
    #exit()
    print(sol_list_new)
    #exit()
    #exit()
    #exit()
    #print(sol_list_new[2])
    print(Afinal)
    print(yfinal)
    #exit()
    for n,s in enumerate(sol_list_new):

        #print(n)
        for i,c in enumerate(Afinal):
            #print(i)
            if np.sum(s*c)%2 != yfinal[i]:
                assert False, 'Wrong solution'      
    #exit()
    print(Afinal)
    print(yfinal)
    print('swap\t', swap_history)
    sol_list_new = enumerate_solution(Afinal, yfinal)

    print('exact # of solutions =',len(sol_list))
    print('GE # of solutions =',len(sol_list_new))
    if len(sol_list) > 0:
        print('marginals init :\t\t', np.mean(np.vstack(sol_list),axis=0))
    if len(sol_list_new) > 0:
        swap_back(sol_list_new, swap_history)
        print('marginals final:\t\t', np.mean(np.vstack(sol_list_new),axis=0))
    exit()

    #exit()

    #print(Anew)
    #print(ynew)
    #exit()
    #sol_list_new = check_all_solution(Anew, ynew)
    print("Number of solutions:", len(sol_list_new))
    #print("here")
    #print(swap_history)
    for i, s in enumerate(sol_list_new):
        swap_back(s, swap_history)
        equiv = False
        for ss in sol_list:
            if np.array_equal(s,ss):
                equiv=True  
        assert equiv, "NOT EQUIVALENT"

    #print(sol_list_new)
    #print(sol)
    exit()
    for i in range(A.shape[0]):
        tmp = np.sum(A[i, :] * sol) %2
        print(tmp == y[i])
    #print(A)
    exit()
    np.random.seed(1)
    y, A = construct_formula(10, 0.5)

    print(A.toarray())
    solve(A, y)

    # ----------- > < 
    # ----------- > <
    exit()
    #print(A.toarray())
    print((True + True) % 2)
    print(np.array([-1,0,1])%2)
    exit()
    tmp = A[0,:].tocoo()
    print(tmp.data)
    print(tmp.row)
    print(tmp.col)
    """ for i in tmp:
        print(i,'\t',tmp) """
    exit()
    for a in A[0,:]:
        last_element = a
    print(last_element.to)
    exit()
    v = A[3,:]
    print(v)
    print('\n\n\n')
    #print(v[
    #print(v)
    #print(A[3,9:])
    #print(A)
    exit()

    print(A.toarray())
    exit()
    print(y)
    print(A)


def solve(A, y):

    # A is the adjacency matrix
    # Gaussian elimination mod 2
    r, c = A.shape

    for j in range(c):
        if A[j,j] == 0:
            swap_diagonal(A, j)
        
        maxi = find_pivot(A, j) # what if 
        print('first pivot: ',maxi)

        while maxi > j: # meaning there are non-zero elements below A[j, j]
            #print((A[maxi,:] - A[j,:]))
            A[maxi,:] = (A[maxi,:] - A[j,:])
            print('data: ',A[maxi,:].data[0])
            A[maxi,:].data = mod2_data(A[maxi,:].data)
            print('data: ',A[maxi,:].data[0])
            exit()
            print(A.toarray())
            #y[maxi].data = (y[maxi] - y[j]).data %2
            maxi = find_pivot(A, j)
            print('second pivot', maxi)
            exit()

    # at this point A should be an upper diagonal matrix
    print(A)
    return
# define sparse matrix class for bool operations ... we want row add, swaps etc.
# do no store zeros ....

if __name__ == "__main__":
    main()