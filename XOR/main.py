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

    return A, y, col_swap

def walkSAT(A, y): 
    # simple walkSAT to find solutions => fix variables at random 
    # and try to find a solution by satisfying one constraint at the time
    return

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

def swap_back(sol, swap_history):
    nswap = len(swap_history)
    for i in range(nswap):
        p1, p2 = swap_history[-(i+1)]
        swap(sol, p1, p2)

def main():

    #np.random.seed(0)
    nr = 7
    nc = 8

    A = np.random.randint(0, 2, nc*nr).reshape(nr, nc)

    y = np.random.randint(0, 2, nr)
    sol_list = check_all_solution(A, y)
    
    print(A)
    print(y)   
    print(sol_list)
    #exit()
    Anew, ynew, swap_history = make_UT(A, y)


    print(Anew)
    print(ynew)
    sol_list_new = check_all_solution(Anew, ynew)
    for s in sol_list_new:
        swap_back(s, swap_history)

    print(sol_list_new)
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