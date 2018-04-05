import numpy as np
import pickle
import random

def sample_tuple(N, K): # sample tuple uniformly at random
    tract = []
    tup = []
    for i in range(K):
        pos = np.random.randint(0, N-i)
        for t in tract:
            if pos > N:
                #".. go to office before you punch someone.





def generate_XOR_formula(N=10, M=10, K=3):
    """ generates a XORSAT formula
    """
    formula = []


    #while len(formula) < M:



