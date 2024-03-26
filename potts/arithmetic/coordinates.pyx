
import cython
from functools import reduce
from itertools import product, chain


# @jit(forceobj=True)
def coordinates(corners):
    """
    Determines all possible lattice coordinates with the given corners.

    Args:
        corners (list): List of lattice corners.

    Returns:
        List of tuples representing lattice coordinates.
    """
    return list(reduce(
        lambda A, B: [
            (left, right) if isinstance(left, int) else left+(right,)
            for left, right in product(A,B)
        ],
        [list(range(corner)) for corner in corners]
    ))


cpdef list flattenAndSortSetUnion(list A):
    cdef list settable
    settable = list(set(chain.from_iterable(A)))
    return list(sorted(settable))


def elementwiseSubtract(A, B): return [A[i]-B[i] for i in range(len(A))]

def elementwiseAddModuli(A, B, moduli): return [(A[i]+B[i])%moduli[i] for i in range(len(A))]

def elementwiseAdd(A, B): return [A[i]+B[i] for i in range(len(A))]

def elementwiseAddOne(A, k): return [A[i]+k for i in range(len(A))]

def subtractMany(A, b): return [elementwiseSubtract(a, b) for a in A]

def binaryEncode(a, n): return tuple(int(s) for s in bin(a)[2:].zfill(n))

def binaryUnencode(t): return sum(t)

def increment(t, at, by=1):
    t = list(t)
    for a in at: t[a] += by
    return tuple(t)
