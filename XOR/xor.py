import numpy as np
import pickle
import random
from matplotlib import pyplot as plt
from collections import Counter

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
    return formula

def generate_sparse(N=10,M=10,K=3):
    formula = generate_XOR_formula(N,M,K)
    A = np.zeros((M,N),dtype=int)
    for i, clause in enumerate(formula):
        for literal in clause:
            A[i, literal]=1
    return A, formula



