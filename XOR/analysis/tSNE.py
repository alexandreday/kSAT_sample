from tsne_visual import TSNE
import numpy as np
import sys, os
import pickle
from sklearn.decomposition import PCA

i_param = int(float(sys.argv[1].split('=')[1])) # specified through the command line !

root_in= '/projectnb/fheating/SAT_GLASS/XORSAT/data/sol/' # root absolute path insert here ...
root_out = '/projectnb/fheating/SAT_GLASS/XORSAT/analysis/TSNE2/'
file_list = '/projectnb/fheating/SAT_GLASS/XORSAT/data/sol/file_name.txt' # absolute path insert here .. 

i_count = 0
fname_out = None
for f in open(file_list,'r'):
    if i_count == i_param:
        fname_in = root_in+f.strip('\n')
        sp = f.strip('\n').split('_')
        param = sp[3:]
        fname_out = root_out+'tSNE_'+'_'.join(sp[3:])
        break
    i_count+=1

print('Reading from %s'%fname_in)
print('Saving in %s'%fname_out)

# ---------> RUNNING TSNE
if fname_out is not None:
    X = np.unpackbits(pickle.load(open(fname_in,'rb'))).astype(int).reshape(10000,-1) 
    X = np.unique(X, axis=0)
    model = TSNE(n_components=2, n_iter=2000, perplexity=50)
    pca = PCA(n_components=100)
    Xpca = pca.fit_transform(X)
    Xtsne = model.fit_transform(Xpca)
    pickle.dump([Xtsne, np.sum(pca.explained_variance_ratio_),len(X)], open(fname_out,'wb'))

