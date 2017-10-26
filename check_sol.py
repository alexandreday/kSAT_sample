from kSAT import KSAT


# -----------> just check solutions --->
model = KSAT()
M=4050
a = 4.05
sfile = 'sol_N=1000_M=%i_alpha=%.2f_K=3.txt'%(M,a)
ffile = 'formula.tmp_N=1000_M=%i_alpha=%a_K=3.cnf'%(M,a)
model.check_all_solution(sfile,ffile)