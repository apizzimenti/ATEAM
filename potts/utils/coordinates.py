
from functools import reduce
from itertools import product
from numba import jit, njit


@jit(forceobj=True)
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
        [list(range(corner+1)) for corner in corners]
    ))


@njit
def elementwiseSubtract(A, B): return [A[i]-B[i] for i in range(len(A))]

@njit
def elementwiseAdd(A, B): return [A[i]+B[i] for i in range(len(A))]

# @jit
def subtractMany(A, b): return [elementwiseSubtract(a, b) for a in A]

# @jit
def binaryEncode(a): return int(''.join([str(t) for t in a]), base=2)

# @jit
def binaryUnencode(n, l): return tuple([int(d) for d in str(bin(n))[2:].zfill(l)])


# @numba.jit
def increment(t, at, by=1):
    t = list(t)
    for a in at: t[a] += by
    return tuple(t)
