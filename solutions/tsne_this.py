from tsne_visual import TSNE
from matplotlib import pyplot as plt
import numpy as np

X=np.loadtxt('sol_N1000=_M=3900_alpha=3.90_K=3.txt',dtype=float)
Xtsne = TSNE(perplexity=50.0, n_iter=1000, angle=0.3).fit_transform(X)
np.savetxt('tsne1.txt', Xtsne)
plt.scatter(Xtsne[:,0], Xtsne[:,1])
plt.show()
