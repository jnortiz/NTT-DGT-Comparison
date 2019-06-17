from dgt import *
from ntt import *

def compare_ntt_dgt():
    print("\nComparing DGT and NTT")
 
    a = gen_polynomial_modq()
    b = gen_polynomial_modq()

    c_NTT = NTT_Mul(a, b)
    c_DGT = dgt_gentlemansande_mul(a, b)
    
   assert(c_NTT == c_transf)

compare_ntt_dgt()