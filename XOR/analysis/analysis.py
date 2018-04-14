import numpy as np
import pickle
from matplotlib import pyplot as plt

def main():

    alpha = np.arange(0.5,1.001,0.01)
    N=100
    K=3

    s = pickle.load(open('../entropy/entropy_N=100_M=100_K=3.pkl','rb'))
    t = pickle.load(open('../time/time_N=100_M=100_K=3.pkl','rb'))

    plt.scatter(alpha,t)
    plt.show()
    exit()
    s = []
    t = []
    for a in alpha:
        M=int(a*N)
        sign = 'N=%i_M=%i_K=%i.pkl'%(N,M,K)






if __name__ == "__main__":
    main()