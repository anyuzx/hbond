import cython
import numpy as np
cimport numpy as np

cdef extern from "c_hbond.h":
    void c_hbondone(double* r_array, double* dim_array, int order, int n, int natmm, int& hbond)
    void c_hbondtwo(double* r1_array, double* r2_array, double* dim_array, int order1, int order2, int n1, int n2, int natmm1, int natmm2, int& hbond)

@cython.boundscheck(False)
@cython.wraparound(False)

def hbondone(np.ndarray[double,ndim = 1,mode="c"] r_array not None, np.ndarray[double,ndim = 1,mode="c"] dim_array not None, int order, int natmm):
    cdef int n
    cdef int hbond
    n = r_array.shape[0]/(3*natmm)
    hbond = 0

    c_hbondone(&r_array[0], &dim_array[0], order, n, natmm, hbond)
    return hbond

def hbondtwo(np.ndarray[double,ndim = 1,mode="c"] r1_array not None, np.ndarray[double,ndim = 1,mode="c"] r2_array not None, np.ndarray[double,ndim = 1,mode="c"] dim_array not None, int order1, int order2, int natmm1, int natmm2):
    cdef int n1
    cdef int n2
    cdef int hbond
    n1 = r1_array.shape[0]/(3*natmm1)
    n2 = r2_array.shape[0]/(3*natmm2)
    hbond = 0

    c_hbondtwo(&r1_array[0], &r2_array[0], &dim_array[0], order1, order2, n1, n2, natmm1, natmm2, hbond)
    return hbond