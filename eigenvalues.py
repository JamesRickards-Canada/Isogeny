from scipy.sparse.linalg import eigs
from scipy.sparse import csr_array
import numpy as np

def secondeval(p, l):
	drc = np.loadtxt("./scipy_adj/"+str(p)+"_"+str(l)+".dat", dtype=int)
	nrows = drc[1][-1] + 1
	M = csr_array((drc[0], (drc[1], drc[2])), dtype='d', shape = (nrows, nrows))
	a, _ = eigs(M, k = 2, which="LR", return_eigenvectors = False)
	return np.real(a)