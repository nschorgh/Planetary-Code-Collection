#! /usr/bin/env python3

# performs singular value decomposition of square view factor matrix

import numpy
from sklearn.utils.extmath import randomized_svd
import time

# number of facets in domain
Nflat = (81-2)**2

# input file name
fn = 'viewfactors.dat'

print('reading input file',fn)
vf = numpy.loadtxt(fn,usecols=numpy.arange(5,Nflat+5))

assert (vf.size == Nflat**2)

print('performing SVD decomposition')
start = time.time()

# full SVD
#U,sigma,VT = numpy.linalg.svd(vf, full_matrices=True, compute_uv=True, hermitian=False)
# time complexity O(Nflat**3) for full matrices, faster for sparse matrices

TRUNC = 200
#TRUNC = Nflat

# truncated approximate SVD
U, sigma, VT = randomized_svd(vf, n_components=TRUNC, n_iter=5, random_state=None)
# time complexity O(Nflat**2*log(TRUNC))

end = time.time()
print('Time',end-start)

print(max(sigma), min(sigma), sum(sigma))


print('saving outputs U,sigma,Vt')
numpy.savetxt('svd_sigma.dat',sigma, newline=' ')
numpy.savetxt('svd_U.dat',U[:,0:TRUNC], fmt='%11.7f')
numpy.savetxt('svd_V.dat',VT[0:TRUNC,:], fmt='%11.7f')

