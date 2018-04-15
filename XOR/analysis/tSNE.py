from tsne_visual import TSNE
import numpy as np
import sys, os
import pickle

i_param = int(float(sys.argv[1].split('=')[1])) # specified through the command line !

root_in= '/projectnb/fheating/SAT_GLASS/XORSAT/data/sol/' # root absolute path insert here ...
root_out = '/projectnb/fheating/SAT_GLASS/XORSAT/analysis/TSNE/'
file_list = '/projectnb/fheating/SAT_GLASS/XORSAT/data/sol/file_name.txt' # absolute path insert here .. 

i_count = 0
fname_out = None
for f in open(file_list,'r'):
    if i_count == i_param:
        fname_in = root_in+f.strip('\n')
        sp = f.split('_')
        param = sp[3:]
        fname_out = root_out+'tSNE_'+'_'.join(sp[3:])
        break
    i_count+=1

print('Reading from %s'%fname_in)
print('Saving in %s'%fname_out)

# ---------> RUNNING TSNE
if fname_out is not None:
    X = np.unpackbits(pickle.load(open(fname_in,'rb'))).astype(int).reshape(10000,-1) 
    model = TSNE(n_components=2, n_iter=2000, perplexity=50)
    Xtsne = model.fit_transform(X)
    pickle.dump(Xtsne, open(fname_out,'wb'))

