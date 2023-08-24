from scipy.sparse.linalg import eigs
from scipy.sparse import csr_array
from scipy.linalg import eig
import numpy as np

#All of these computations assume that ssl_graph_scipy(p, l) has already been called. l can be an array.

#Returns the adjacency matrix as a sparse array
def sparseadjmtx(p, l):
	filename = "./scipy_adj/" + str(p) + "_"
	if type(l) == type(1):
		filename = filename + str(l)+ ".dat"
	else:
		filename = filename + str(l[0])
		for lval in l[1:]:
			filename = filename + "-" + str(lval)
		filename = filename + ".dat"
	drc = np.loadtxt(filename, dtype=int)
	if type(drc) == type(drc[1]):
		nrows = drc[1][-1] + 1
	else:#Happens when only 1 vertex
		return
	if nrows <= 3:#Too small for sparse array.
		return
	return csr_array((drc[0], (drc[1], drc[2])), dtype='d', shape = (nrows, nrows))
	
#Returns the second largest eigenvalue of the graph
def secondeval(p, l):
	M = sparseadjmtx(p, l)
	a, _ = eigs(M, k = 2, which = "LR", return_eigenvectors = False)
	return np.real(a)

#Returns the second largest eigenvalue in terms of absolute value of the graph
def secondabseval(p, l):
	M = sparseadjmtx(p, l)
	a, _ = eigs(M, k = 2, which = "LR", return_eigenvectors = False)
	b = eigs(M, k = 1, which = "SR", return_eigenvectors = False)[0]
	a = np.abs(np.real(a))
	b = np.abs(np.real(b))
	if a < b:
		return b
	return a

#Returns all the eigenvalues of the graph
def allevals(p, l):
	M = sparseadjmtx(p, l)
	ev = eig(M.toarray(), right=False)
	sorted = np.array([np.real(x) for x in ev])
	sorted.sort()
	return sorted

