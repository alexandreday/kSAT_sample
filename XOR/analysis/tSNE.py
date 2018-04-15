from tsne_visual import TSNE
import numpy as np
import sys, os

i_param = int(sys.argv[1].split('=')[1]) # specified through the command line !

root= '' # root absolute path insert here ...
file_list = '' # absolute path insert here .. 

i_count = 0
fname_out = None
for f in open(file_list,'r'):
    if i_count == i_param:
        fname_in = f
        sp = f.split('_')
        param = sp[3:]
        fname_out = 'tSNE_'+'_'.join(sp[3:])
        break

# ---------> RUNNING TSNE
if fname_out is not None:
    X = np.unpackbits(pickle.load(open(fname_in,'rb'))).astype(int).reshape(10000,-1) 
    model = TSNE(n_components=2, n_iter=2000, perplexity=50)
    Xtsne = model.fit_transform(X)
    pickle.dump(Xtsne, open(fname_out,'wb'))

