import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

def main():
    np.random.seed(1)
    y, A = construct_formula(10,0.5)
    print(A.toarray())

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
    r,c = A.shape

    for j in range(c):
       element_j = A[j:,j]
       if element_j.shape[1] == 1:
           print("empty")
    maxj = element_j[-1][0]

        


    return



def construct_formula(N, alpha, p = 3): # 3-XORSAT
    M=int(alpha*N)
    y = np.random.randint(0,2,size=M)

    A = []
    while nc < M:
        c = np.random.randint(0,N,p)
        if len(np.unique(c)) == p:
            

    A = np.random.randint(0,N,size=3*M).reshape(M,3)
    idx_1 = []
    idx_2 = []
    for i, c in enumerate(A):
        for j, v in enumerate(c):
            idx_1.append(i)
            idx_2.append(v)
    data = np.ones(len(idx_2))
    C = coo_matrix((data,(idx_1,idx_2)), dtype=np.int8, shape=(len(A),N))

    return y, C.tocsc()

if __name__ == "__main__":
    main()