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

def sample_solution_GE(A, y, n_sample=10, time_max = 86400):
    """ A has to be in a reduced form """
    M, N = A.shape
    t_init = time.time()
    
    pivots = np.where(np.diagonal(A) == 1)[0] # pivots are identified
    none_pivots = np.setdiff1d(np.arange(N), pivots)
    none_pivots_2 = np.setdiff1d(np.arange(M), pivots)

    rank = len(pivots) # matrix rank
    xsol = -1*np.ones(N, dtype=int) # unconstrained variables are marked by a -2

    N_free = N - rank # upper bound on number of log of number of solutions (but may be fewer !)
 
    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=N_free)), dtype=np.int)

    pivot_set = set(list(pivots))

    if n_sample > 2**N_free:
        assert False, "Make sure you don't ask for too many samples !"
    else:
        all_sol = []        
        while len(all_sol) < n_sample:
            print(len(all_sol))
            is_sol = False
            xsol = -1*np.ones(N, dtype=int)
            xsol[none_pivots] = np.random.randint(0, 2, N_free)
            y_res = np.remainder(np.dot(A[:,none_pivots], xsol[none_pivots].T) + y, 2)

            if np.count_nonzero(y_res[none_pivots_2] == 1) == 0:
                xsol[pivots] = y_res[pivots]
                is_sol = True
            if is_sol:
                all_sol.append(xsol)
            if time.time() - t_init > time_max:
                break

                
    return all_sol

def main():
    import time

    t_all = []
    t_me = []
    from xor import generate_sparse
    N = 100 # number of literals / variables
    M = 25 # number of constraints
    K = 3
    A, f = generate_sparse(N=N, M=M, K=K)
    y = np.random.randint(0, 2, M) # random constraints

    make_diagonal(A, y)
    #print(A)
    sol = sample_solution_GE(A, y, n_sample=50)

    print(sol)
    exit()


    n_col = N
    n_row = M
    np.random.seed(93)
    for i in range(100):
        #print(i)
        A, f = generate_sparse(N=N, M=M, K=K)
        y = np.random.randint(0, 2, n_row)

        t = time.time()
        #sol_init = solve_ES(A, y)
        t_all.append(time.time()-t)

        t = time.time()
        make_diagonal(A, y)
        sol_final = enumerate_solution_v2(A, y)
        print(len(sol_final))
        t_me.append(time.time()-t)

    exit()

    if len(sol_final)>0:
        n_star = np.count_nonzero(sol_final[0] == -1)

    assert len(sol_init) == len(sol_final)*2**n_star

    if len(sol_init) > 0:
        marg_final = marginals(sol_final)
        marg_final[marg_final < -0.0001] = 0.5 # star variables
        if abs(np.linalg.norm(marginals(sol_init) - marg_final)) > 1e-5:
            print(i)
            assert False
        #print(marginals(sol_init))
        #print(marginals(sol_final))

    print('all:', np.mean(t_all)) # this scales like 2^N_variable
    print('me:', np.mean(t_me)) 
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