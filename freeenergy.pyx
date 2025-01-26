# freeenergy.pyx
# cython: boundscheck=False, wraparound=False

def compute_sum(double[:] array):
    cdef double total = 0.0
    cdef int i
    for i in range(array.shape[0]):
        total += array[i]
    return total

import numpy as np
cimport numpy as np
from libc.math cimport log

# Define types for numpy arrays
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def compute_point_cluster_energy(np.ndarray[DTYPE_t, ndim=4] prob, str lattice="BCC"):
    cdef DTYPE_t para
    if lattice == "BCC":
        para = -0.25
    elif lattice == "FCC":
        para = 1.25
    else:
        raise ValueError("Invalid lattice type. Choose 'BCC' or 'FCC'.")

    cdef DTYPE_t a0, a1, b0, b1, c0, c1, d0, d1

    # Compute sums
    a0 = (prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[0, 0, 1, 0] + prob[0, 0, 1, 1] +
          prob[0, 1, 0, 0] + prob[0, 1, 0, 1] + prob[0, 1, 1, 0] + prob[0, 1, 1, 1])
    a1 = (prob[1, 0, 0, 0] + prob[1, 0, 0, 1] + prob[1, 0, 1, 0] + prob[1, 0, 1, 1] +
          prob[1, 1, 0, 0] + prob[1, 1, 0, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1])

    b0 = (prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[0, 0, 1, 0] + prob[0, 0, 1, 1] +
          prob[1, 0, 0, 0] + prob[1, 0, 0, 1] + prob[1, 0, 1, 0] + prob[1, 0, 1, 1])
    b1 = (prob[0, 1, 0, 0] + prob[0, 1, 0, 1] + prob[0, 1, 1, 0] + prob[0, 1, 1, 1] +
          prob[1, 1, 0, 0] + prob[1, 1, 0, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1])

    c0 = (prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[0, 1, 0, 0] + prob[0, 1, 0, 1] +
          prob[1, 0, 0, 0] + prob[1, 0, 0, 1] + prob[1, 1, 0, 0] + prob[1, 1, 0, 1])
    c1 = (prob[0, 0, 1, 0] + prob[0, 0, 1, 1] + prob[0, 1, 1, 0] + prob[0, 1, 1, 1] +
          prob[1, 0, 1, 0] + prob[1, 0, 1, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1])

    d0 = (prob[0, 0, 0, 0] + prob[0, 0, 1, 0] + prob[0, 1, 0, 0] + prob[0, 1, 1, 0] +
          prob[1, 0, 0, 0] + prob[1, 0, 1, 0] + prob[1, 1, 0, 0] + prob[1, 1, 1, 0])
    d1 = (prob[0, 0, 0, 1] + prob[0, 0, 1, 1] + prob[0, 1, 0, 1] + prob[0, 1, 1, 1] +
          prob[1, 0, 0, 1] + prob[1, 0, 1, 1] + prob[1, 1, 0, 1] + prob[1, 1, 1, 1])

    # Compute result
    cdef DTYPE_t result = (
        a0 * np.log(a0) + a1 * np.log(a1) +
        b0 * np.log(b0) + b1 * np.log(b1) +
        c0 * np.log(c0) + c1 * np.log(c1) +
        d0 * np.log(d0) + d1 * np.log(d1)
    )

    # Create point_prob array
    cdef np.ndarray[DTYPE_t, ndim=2] point_prob = np.array([
        [a0, b0, c0, d0],
        [a1, b1, c1, d1]
    ], dtype=DTYPE)

    return result * para, point_prob

def compute_ternary_cluster_energy(np.ndarray[DTYPE_t, ndim=4] prob, str lattice="BCC"):
    cdef DTYPE_t para
    if lattice == "BCC":
        para = -3.0
    else:
        para = 0.0

    cdef DTYPE_t abc000, abc001, abc010, abc011, abc100, abc101, abc110, abc111
    cdef DTYPE_t abd000, abd001, abd010, abd011, abd100, abd101, abd110, abd111
    cdef DTYPE_t acd000, acd001, acd010, acd011, acd100, acd101, acd110, acd111
    cdef DTYPE_t bcd000, bcd001, bcd010, bcd011, bcd100, bcd101, bcd110, bcd111

    # Use memory view for faster access
    cdef DTYPE_t[:, :, :, :] prob_view = prob

    # Compute sums
    abc000 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 0, 1]
    abc001 = prob_view[0, 0, 1, 0] + prob_view[0, 0, 1, 1]
    abc010 = prob_view[0, 1, 0, 0] + prob_view[0, 1, 0, 1]
    abc011 = prob_view[0, 1, 1, 0] + prob_view[0, 1, 1, 1]
    abc100 = prob_view[1, 0, 0, 0] + prob_view[1, 0, 0, 1]
    abc101 = prob_view[1, 0, 1, 0] + prob_view[1, 0, 1, 1]
    abc110 = prob_view[1, 1, 0, 0] + prob_view[1, 1, 0, 1]
    abc111 = prob_view[1, 1, 1, 0] + prob_view[1, 1, 1, 1]

    abd000 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 1, 0]
    abd001 = prob_view[0, 0, 0, 1] + prob_view[0, 0, 1, 1]
    abd010 = prob_view[0, 1, 0, 0] + prob_view[0, 1, 1, 0]
    abd011 = prob_view[0, 1, 0, 1] + prob_view[0, 1, 1, 1]
    abd100 = prob_view[1, 0, 0, 0] + prob_view[1, 0, 1, 0]
    abd101 = prob_view[1, 0, 0, 1] + prob_view[1, 0, 1, 1]
    abd110 = prob_view[1, 1, 0, 0] + prob_view[1, 1, 1, 0]
    abd111 = prob_view[1, 1, 0, 1] + prob_view[1, 1, 1, 1]

    acd000 = prob_view[0, 0, 0, 0] + prob_view[0, 1, 0, 0]
    acd001 = prob_view[0, 0, 0, 1] + prob_view[0, 1, 0, 1]
    acd010 = prob_view[0, 0, 1, 0] + prob_view[0, 1, 1, 0]
    acd011 = prob_view[0, 0, 1, 1] + prob_view[0, 1, 1, 1]
    acd100 = prob_view[1, 0, 0, 0] + prob_view[1, 1, 0, 0]
    acd101 = prob_view[1, 0, 0, 1] + prob_view[1, 1, 0, 1]
    acd110 = prob_view[1, 0, 1, 0] + prob_view[1, 1, 1, 0]
    acd111 = prob_view[1, 0, 1, 1] + prob_view[1, 1, 1, 1]

    bcd000 = prob_view[0, 0, 0, 0] + prob_view[1, 0, 0, 0]
    bcd001 = prob_view[0, 0, 0, 1] + prob_view[1, 0, 0, 1]
    bcd010 = prob_view[0, 0, 1, 0] + prob_view[1, 0, 1, 0]
    bcd011 = prob_view[0, 0, 1, 1] + prob_view[1, 0, 1, 1]
    bcd100 = prob_view[0, 1, 0, 0] + prob_view[1, 1, 0, 0]
    bcd101 = prob_view[0, 1, 0, 1] + prob_view[1, 1, 0, 1]
    bcd110 = prob_view[0, 1, 1, 0] + prob_view[1, 1, 1, 0]
    bcd111 = prob_view[0, 1, 1, 1] + prob_view[1, 1, 1, 1]

    # Compute the sum of x * log(x) for each term
    cdef DTYPE_t result = (
        abc000 * log(abc000) + abc001 * log(abc001) + abc010 * log(abc010) + abc011 * log(abc011) +
        abc100 * log(abc100) + abc101 * log(abc101) + abc110 * log(abc110) + abc111 * log(abc111) +
        abd000 * log(abd000) + abd001 * log(abd001) + abd010 * log(abd010) + abd011 * log(abd011) +
        abd100 * log(abd100) + abd101 * log(abd101) + abd110 * log(abd110) + abd111 * log(abd111) +
        acd000 * log(acd000) + acd001 * log(acd001) + acd010 * log(acd010) + acd011 * log(acd011) +
        acd100 * log(acd100) + acd101 * log(acd101) + acd110 * log(acd110) + acd111 * log(acd111) +
        bcd000 * log(bcd000) + bcd001 * log(bcd001) + bcd010 * log(bcd010) + bcd011 * log(bcd011) +
        bcd100 * log(bcd100) + bcd101 * log(bcd101) + bcd110 * log(bcd110) + bcd111 * log(bcd111)
    )

    return result * para

def compute_binary_cluster_energy(np.ndarray[DTYPE_t, ndim=4] prob, str lattice="BCC"):
    cdef np.ndarray[DTYPE_t, ndim=1] para
    if lattice == "BCC":
        para = np.array([1.5, 1], dtype=DTYPE)
    elif lattice == "FCC":
        para = np.array([-1, -1], dtype=DTYPE)
    else:
        raise ValueError("Invalid lattice type. Choose 'BCC' or 'FCC'.")

    cdef DTYPE_t ab00, ab01, ab10, ab11
    cdef DTYPE_t ac00, ac01, ac10, ac11
    cdef DTYPE_t ad00, ad01, ad10, ad11
    cdef DTYPE_t bc00, bc01, bc10, bc11
    cdef DTYPE_t bd00, bd01, bd10, bd11
    cdef DTYPE_t cd00, cd01, cd10, cd11

    # Use memory view for faster access
    cdef DTYPE_t[:, :, :, :] prob_view = prob

    # Compute sums
    ab00 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 0, 1] + prob_view[0, 0, 1, 0] + prob_view[0, 0, 1, 1]
    ab01 = prob_view[0, 1, 0, 0] + prob_view[0, 1, 0, 1] + prob_view[0, 1, 1, 0] + prob_view[0, 1, 1, 1]
    ab10 = prob_view[1, 0, 0, 0] + prob_view[1, 0, 0, 1] + prob_view[1, 0, 1, 0] + prob_view[1, 0, 1, 1]
    ab11 = prob_view[1, 1, 0, 0] + prob_view[1, 1, 0, 1] + prob_view[1, 1, 1, 0] + prob_view[1, 1, 1, 1]

    ac00 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 0, 1] + prob_view[0, 1, 0, 0] + prob_view[0, 1, 0, 1]
    ac01 = prob_view[0, 0, 1, 0] + prob_view[0, 0, 1, 1] + prob_view[0, 1, 1, 0] + prob_view[0, 1, 1, 1]
    ac10 = prob_view[1, 0, 0, 0] + prob_view[1, 0, 0, 1] + prob_view[1, 1, 0, 0] + prob_view[1, 1, 0, 1]
    ac11 = prob_view[1, 0, 1, 0] + prob_view[1, 0, 1, 1] + prob_view[1, 1, 1, 0] + prob_view[1, 1, 1, 1]

    ad00 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 1, 0] + prob_view[0, 1, 0, 0] + prob_view[0, 1, 1, 0]
    ad01 = prob_view[0, 0, 0, 1] + prob_view[0, 0, 1, 1] + prob_view[0, 1, 0, 1] + prob_view[0, 1, 1, 1]
    ad10 = prob_view[1, 0, 0, 0] + prob_view[1, 0, 1, 0] + prob_view[1, 1, 0, 0] + prob_view[1, 1, 1, 0]
    ad11 = prob_view[1, 0, 0, 1] + prob_view[1, 0, 1, 1] + prob_view[1, 1, 0, 1] + prob_view[1, 1, 1, 1]

    bc00 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 0, 1] + prob_view[1, 0, 0, 0] + prob_view[1, 0, 0, 1]
    bc01 = prob_view[0, 0, 1, 0] + prob_view[0, 0, 1, 1] + prob_view[1, 0, 1, 0] + prob_view[1, 0, 1, 1]
    bc10 = prob_view[0, 1, 0, 0] + prob_view[0, 1, 0, 1] + prob_view[1, 1, 0, 0] + prob_view[1, 1, 0, 1]
    bc11 = prob_view[0, 1, 1, 0] + prob_view[0, 1, 1, 1] + prob_view[1, 1, 1, 0] + prob_view[1, 1, 1, 1]

    bd00 = prob_view[0, 0, 0, 0] + prob_view[0, 0, 1, 0] + prob_view[1, 0, 0, 0] + prob_view[1, 0, 1, 0]
    bd01 = prob_view[0, 0, 0, 1] + prob_view[0, 0, 1, 1] + prob_view[1, 0, 0, 1] + prob_view[1, 0, 1, 1]
    bd10 = prob_view[0, 1, 0, 0] + prob_view[0, 1, 1, 0] + prob_view[1, 1, 0, 0] + prob_view[1, 1, 1, 0]
    bd11 = prob_view[0, 1, 0, 1] + prob_view[0, 1, 1, 1] + prob_view[1, 1, 0, 1] + prob_view[1, 1, 1, 1]

    cd00 = prob_view[0, 0, 0, 0] + prob_view[0, 1, 0, 0] + prob_view[1, 0, 0, 0] + prob_view[1, 1, 0, 0]
    cd01 = prob_view[0, 0, 0, 1] + prob_view[0, 1, 0, 1] + prob_view[1, 0, 0, 1] + prob_view[1, 1, 0, 1]
    cd10 = prob_view[0, 0, 1, 0] + prob_view[0, 1, 1, 0] + prob_view[1, 0, 1, 0] + prob_view[1, 1, 1, 0]
    cd11 = prob_view[0, 0, 1, 1] + prob_view[0, 1, 1, 1] + prob_view[1, 0, 1, 1] + prob_view[1, 1, 1, 1]

    # Compute results
    cdef DTYPE_t result1 = (
        ab00 * log(ab00) + ab01 * log(ab01) + ab10 * log(ab10) + ab11 * log(ab11) +
        cd00 * log(cd00) + cd01 * log(cd01) + cd10 * log(cd10) + cd11 * log(cd11)
    )

    cdef DTYPE_t result2 = (
        ac00 * log(ac00) + ac01 * log(ac01) + ac10 * log(ac10) + ac11 * log(ac11) +
        ad00 * log(ad00) + ad01 * log(ad01) + ad10 * log(ad10) + ad11 * log(ad11) +
        bc00 * log(bc00) + bc01 * log(bc01) + bc10 * log(bc10) + bc11 * log(bc11) +
        bd00 * log(bd00) + bd01 * log(bd01) + bd10 * log(bd10) + bd11 * log(bd11)
    )

    return para[0] * result1 + para[1] * result2

#python setup.py build_ext --inplace