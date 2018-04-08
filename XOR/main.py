import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

def check_solution(A, y, sol):
    nclause = A.shape[0]
    is_solution = True
    for i in range(nclause):
        if np.sum(A[i, :] * sol) % 2 != y[i]:
            is_solution = False
            break
    return is_solution

def swap_rows(A, ri, rj):
    tmp = np.copy(A[ri, :]) # careful when slicing, will return a view
    A[ri, :] = A[rj, :]
    A[rj, :] = tmp
   
def swap(y, i, j):
    tmp = y[i]
    y[i] = y[j]
    y[j] = tmp

def swap_columns(A, ci, cj):
    tmp = np.copy(A[:, ci]) # careful when slicing, will return a view
    A[:, ci] = A[:, cj]
    A[:, cj] = tmp

def add_to_row(A, ri, rj): # A[ri, :] <- A[ri, :] + A[rj, :]
    A[ri, :] = A[ri, :] + A[rj, :]

def make_UT(A_, y_):
    """ Transforms A_ into a upper triangular matrix (where A.x == y) is a linear system of
    equations in base 2. 

    Returns
    -------
    (A, y, swap_history) : tuple
        A in the UT form
        the equivalent y associated with the new A
        The permutations to do on the original index (x1,x2,x3, etc.) to retrive the original basis.
    """


    col_swap = []
    A = np.copy(A_)
    y = np.copy(y_)

    nr, nc = A.shape
    
    for j in range(nc):
        if j < nr:
            #A.get_first_index_non_zero(col = j)
            pos = np.where(A[:,j] == 1)[0]

            if len(pos) > 0: # there are non-zero elements in this column
                if A[j,j] == 0: # diagonal element is zero -> need pivot to be 1
                    if pos[-1] > j:
                        swap_rows(A, j, pos[-1])
                        swap(y, j, pos[-1])
                        pos = pos[:-1] 
                    else: # already UT !
                        pos = []
                        pos_row = np.where(A[j,:] == 1)[0]
                        if len(pos_row) > 0:
                            swap_columns(A, j, pos_row[0])
                            col_swap.append([j, pos_row[0]])
                            pos = np.where(A[:,j] == 1)[0]
                            #swap(y, j, pos_row[0])                        
                for p in pos:
                    if p > j: # making elements zero => and only working with lower triangle
                        add_to_row(A, p, j)
                        y[p] += y[j]
            A = np.remainder(A, 2)
            y = np.remainder(y, 2)

    return A, y, col_swap # do we care about this => maybe, maybe not.

def make_diag(A_, y_):
    A = np.copy(A_)
    y = np.copy(y_)
    M, N = A.shape
    # starting from last row, reduce:
    for i in reversed(range(M)):
        #print(i)
        if A[i,i] == 1:
            pos = np.where(A[:,i] == 1)[0]
            for r in pos:
                if r < i:
                    A[r, :] -= A[i, :]
                    y[r] -= y[i]

        A = np.remainder(A, 2)
        y = np.remainder(y, 2)

    return A, y
        
def check_all_solution(A, y):
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
    A = np.copy(Areduce)
    y = np.copy(y_)
    M, N = Areduce.shape
    N_free = N-M
    N_check = 2**N_free
    b2_array = lambda n10 : np.array(list(np.binary_repr(n10, width=N_free)), dtype=np.int)

    test = np.sum(A, axis = 1)
    free_add = np.count_nonzero(y[test == 0] == 0)
    if np.count_nonzero(y[test == 0] == 1) > 0:
        return []
    
    all_sol = []
    N_check = 2**(N_free+free_add)

    for i in range(N_check): # does not always work, if a column is full of zeros
        xsol = np.zeros(N,dtype=int)
        xsol[-N_free:] = b2_array(i)
        is_sol = True
        for i in reversed(range(M)):
            if np.sum(A[i,:]*xsol) % 2 == y[i]:
                bit = 0
            else:
                bit = 1

            if A[i,i] == 1:
                xsol[i] = bit
            else:
                if bit == 1:
                    is_sol = False
                # cannot be a solution ... 
                xsol[i] = (bit+1)%2 # two possible solutions

        if is_sol:
            all_sol.append(xsol)
        
    return all_sol

def swap_back(sol, swap_history):
    nswap = len(swap_history)
    for i in range(nswap):
        p1, p2 = swap_history[-(i+1)]
        swap(sol, p1, p2)

def main():
    from xor import generate_sparse

    N = 8
    M = 6
    K = 3
    n_col = N
    n_row = M
    np.random.seed(5)
    A, f = generate_sparse(N=N, M=M, K=K)
    y = np.random.randint(0, 2, n_row)
    sol_list = check_all_solution(A, y)
    print('marginals:\t', np.mean(np.vstack(sol_list),axis=0))

    #print(np.sum(A,axis=1))
    #exit()
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
    sol_list_new = enumerate_solution(Afinal, yfinal)
    print(sol_list_new)
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
    print('swap\t',swap_history)
    sol_list_new = enumerate_solution(Afinal, yfinal)

    print('marginals init:\t\t', np.mean(np.vstack(sol_list),axis=0))
    print('marginals final:\t\t', np.mean(np.vstack(sol_list_new),axis=0))

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

class SparseBool:

    def __init__(self):
        print('etc. etc.')
    
def mod2_data(data):
    return [list(list(np.array(data[0],dtype=np.int8) %2))]

def swap_diagonal(A, j):
    i = find_pivot(A,j)
    swap_rows(A,i,j)

def find_pivot(A, j):
    # looks for the first non-zero in column j starting from the bottom
    # returns the row index of that element
    return A[j:,j].tocoo().row[-1]

""" def swap_rows(A, i, j):
    tmp = np.copy(A[i,:])
    A[i,:] = A[j,:]
    A[j,:] = tmp
 """
def construct_formula(N, alpha, p = 3): # 3-XORSAT
    M=int(alpha*N)
    y = np.random.randint(0,2,size=M)

    A = []
    while len(A) < M:
        c = np.random.randint(0,N,p)
        if len(np.unique(c)) == p:
            A.append(c)

    idx_1 = []
    idx_2 = []
    for i, c in enumerate(A):
        for j, v in enumerate(c):
            idx_1.append(i)
            idx_2.append(v)
    data = np.ones(len(idx_2))
    C = coo_matrix((data,(idx_1,idx_2)), dtype=np.int8, shape=(len(A),N))

    return y, C.tolil()

if __name__ == "__main__":
    main()