from kSAT import KSAT


# -----------> just check solutions --->
model = KSAT()
model.check_all_solution('sol_N=1000_M=3500_alpha=3.50_K=3.txt','formula.tmp_N=1000_M=3500_alpha=3.50_K=3.cnf')
