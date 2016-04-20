cimport cython
cimport numpy as np
import numpy as np
cdef double f_typed(double x) except? -2:
    return x * (x - 1)
cpdef double integrate_f_typed(double a, double b, int N):
    cdef int i
    cdef double s, dx
    s = 0
    dx = (b - a) / N
    for i in range(N):
        s += f_typed(a + i * dx)
    return s * dx
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[double] apply_integrate_f_wrap(np.ndarray[double] col_a, np.ndarray[double] col_b, np.ndarray[int] col_N):
    cdef int i, n = len(col_N)
    assert len(col_a) == len(col_b) == n
    cdef np.ndarray[double] res = np.empty(n)
    for i in range(n):
        res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
    return res
