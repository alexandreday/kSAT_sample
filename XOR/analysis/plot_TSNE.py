from matplotlib import pyplot as plt
import numpy as np
import pickle
from fdc import FDC, plotting
from sklearn.preprocessing import StandardScaler as SS

N=4000
K=3
alpha_range = np.arange(0.5,1.001,0.05)

root= '/Users/alexandreday/GitProject/kSAT_sample/XOR/data/TSNE2/'
root_out= '/Users/alexandreday/GitProject/kSAT_sample/XOR/analysis/plots/'
for a in alpha_range:
    print(a)
    f = 'tSNE_N=%i_a=%.3f_K=%i.pkl'%(N,a,K)
    X, pca_r = pickle.load(open(root+f,'rb'))
    X = SS().fit_transform(X)
    modelfdc=FDC(eta = 0.0)
    modelfdc.fit(X)

    plotting.cluster_w_label(X, modelfdc.cluster_label)

    """ plt.scatter(X[:,0], X[:,1], s=3, alpha=0.5)
    plt.title('$N=%i,\\alpha=%.3f, K=%i$'%(N,a,K))
    plt.xticks([])
    plt.yticks([])
    fout = f.strip('.pkl')+'.pdf'
    plt.tight_layout()
    plt.savefig(root_out+fout)
    plt.clf() """
    #plt.show()
