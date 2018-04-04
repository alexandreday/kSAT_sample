import numpy as np
from scipy.sparse import coo_matrix, csc_matrix


def swap_rows(A, ri, rj):
    tmp = A[ri, :]
    A[ri, :] = A[rj, :]
    A[rj, :] = tmp

def add_to_row(A, ri, rj): # A[ri, :] <- A[ri, :] + A[rj, :]
    A[ri, :] = A[ri, :] + A[rj, :]

def make_UT(A_, y_):
    A = np.copy(A_)
    y = np.copy(y_)

    nr, nc = A.shape

    for j in range(nc):
        #A.get_first_index_non_zero(col = j)
        pos = np.where(A[:,j] == 1)[0]

        if len(pos) > 0: # there are non-zero elements in this column
            if A[j,j] == 0:
                if len(pos) > 0:
                    if pos[-1] > j:
                        swap_rows(A, j, pos[-1])
                        pos = pos[:-1] 
                    else: # already UT !
                        pos = []
            for p in pos:
                if p != j: # making elements zero !
                    add_to_row(A, p, j)
                    y[p] += y[j]
        A = np.remainder(A, 2)
        print(A)
        if j == 4:
            exit()
        #exit()
        y = np.remainder(y, 2)
    
    return A, y

def main():

    np.random.seed(0)
    A = np.random.randint(0, 2, 40).reshape(8, 5)
    y = np.random.randint(0,2, 8)
    print(A)
    print(y)
    A, y = make_UT(A, y)
    print(A)
    print(y)

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

def swap_rows(A, i, j):
    tmp = np.copy(A[i,:])
    A[i,:] = A[j,:]
    A[j,:] = tmp

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