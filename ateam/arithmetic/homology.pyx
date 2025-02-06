
import cython
import numpy as np
from itertools import product


@cython.boundscheck(False)
@cython.wraparound(False)
def isNullHomologous(A, x, I, includes=[]):
    """
    Checks whether the provided sequence of faces is a boundary of some chain of
    plaquettes.

    Args:
        A (galois.Array): Complete boundary matrix.
        x (galois.Array): Chain of faces.
        I (galois.Array): Identity matrix of dimension \((n+1) \\times (n+1)\),
            where \(n\) is the number of columns of \(A\).
        includes (iterable=[]): Column indices to _include_.

    Returns:
        A boolean which is truthy when \(x\) is nullhomologous and falsy otherwise.
    """
    if len(includes): B = A.take(includes, axis=1)
    else: B = A
    
    # Augment the matrix, get the RREF, and check whether we have a solution or
    # not. If not, the RREF of the matrix is the identity of dimension one
    # larger than the boundary submatrix.
    R = (np.c_[B, x].astype(np.int32)).row_reduce()
    _, m = B.shape
    return not (R[:m+1,:] == I[:m+1,:m+1]).all()


@cython.boundscheck(False)
@cython.wraparound(False)
def sampleFromKernel(A, F, includes=[], relativeCells=None, relativeFaces=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.

    Args:
        A (np.ndarray): Complete coboundary matrix.
        F (galois.GF): Finite field representation.
        includes (iterable=[]): Column indices to _include_.
        relativeCells (iterable=None): Column indices to _include_. (Possible
            duplicate of the above?)
        relativeFaces (iterable=None): Row indices to _include_.

    Returns:
        A cochain on the \((k-1)\)-skeleton of the complex that is a cocycle on
        whatever subcomplex is specified by `includes`, `relativeCells`, and
        `relativeFaces`, represented as a `galois.FieldArray`.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if len(includes): B = A.take(includes, axis=0)
    else: B = A

    # Relative homology.
    if relativeCells: C = B.take(relativeCells, axis=0)
    else: C = B

    if relativeFaces: D = C.take(relativeFaces, axis=1)
    else: D = C
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = D.null_space()
    M, _ = K.shape
    Q = F.Random(M)
    return (Q@K)


@cython.boundscheck(False)
@cython.wraparound(False)
def evaluateCochain(boundary, spins):
    """
    Evaluates the coboundary of a cochain on the \(k\)-skeleton.

    Args:
        boundary (np.ndarray): Boundary of the \(k\)-skeleton.
        spins (galois.FieldArray): Cochain (spins) on the \((k-1)\)-skeleton.
    
    Returns:
        A `galois.FieldArray` of coefficients on the \(k\)-skeleton.
    """
    evaluation = spins[boundary]
    evaluation[:, 1::2] = -evaluation[:, 1::2]
    return evaluation.sum(axis=1)


def autocorrelation(data):
    """
    Computes the autocorrelation of a given observable.

    Args:
        data (Iterable): An iterable, indexed by sample times, containing data
            from a given observable.
    
    Returns:
        An `np.ndarray` of autocorrelation data.
    """

    # Expected value (i.e. sample mean, which converges to the expectation by
    # LLN).
    mu = data.mean()
    normalized = data-mu
    N = len(data)

    autocorrs = np.array([
        np.dot(normalized[t:], normalized[:N-t])*(1/N) for t in range(N)
    ])

    return autocorrs/autocorrs[0]


@cython.boundscheck(False)
@cython.wraparound(False)
def essentialCyclesBorn(
        phatBoundary, coboundary, boundary, reindexed, tranches, homology, F, spins,
        times, indices, lower, highest, stop
    ):
    """
    Computes the persistent homology of the given complex, identifying when the
    first nontrivial element of the parent space (generally a torus) is born.

    Args:
        phatBoundary (phat.boundary_matrix): An initialized PHAT boundary matrix
            with un-set columns.
        coboundary (galois.FieldArray): Coboundary matrix for the complex.
        boundary (dict): Indexed boundary matrices associated with a `Lattice`
            object.
        reindexed (dict): Reindexed version of `boundary`.
        tranches (dict): Dictionary specifying the indices at which cubes of
            each dimension appear in the flattened boundary matrix (i.e. the
            columns of `phatBoundary`).
        homology (int): Homology group \(k\).
        F (galois.GF): Finite field representation.
        spins (galois.FieldArray): Cochain on the \((k-1)\)-skeleton.
        times (np.ndarray): Index for the filtration.
        indices (np.ndarray): Indices for plaquettes.
        lower (list): List of (dimension, boundary) pairs for all cubes of
            dimension less than \(k\).
        highest (list): List of (dimension, boundary) pairs for all cubes of
            dimension greater than \(k+1\).
        stop (int): How many essential cycles are identified before we we-sample
            a cochain.

    Returns:
        A triplet: a \((k-1)\)-cochain as a `galois.FieldArray`; a binary
        \(\\beta \\times n\) `numpy.ndarray`, where each row has a \(1\) at the
        indices of occupied plaquettes when homological percolation occurred,
        \(\\beta\) is the rank of the \(k\)th homology group, and \(n\) is the
        number of plaquettes; a length-\(n\) binary `numpy.ndarray` with ones
        at the indices of satisfied plaquettes.
    """
    # See which faces' spins sum to 0.
    cycles = evaluateCochain(boundary[homology], spins)
    satisfiedIndices = (cycles == 0).nonzero()[0]
    unsatisfiedIndices = (cycles > 0).nonzero()[0]

    # Randomize the order of the satisfied indices, but keep the unsatisfied
    # ones fixed. Then, reindex (homology+1)-dimensional objects' boundaries to
    # reflect the shuffle.
    shuffled = np.random.permutation(satisfiedIndices)
    unordered = np.concatenate([shuffled, unsatisfiedIndices])
    _reindexer = np.array([unordered, np.arange(len(cycles))])
    reindexer = _reindexer[:, _reindexer[0].argsort()][1]
    ordered = indices[reindexer]

    # Reorder the actual boundary matrices, then sort both; tack the other ones
    # on, too.
    target = list(product(
        [homology], np.sort(reindexed[homology][unordered]).tolist()
    ))

    higher = list(product(
        [homology+1], np.sort(ordered[boundary[homology+1]]).tolist()
    ))

    # Compute persistence pairs, and find the first time an essential cycle is
    # found.
    phatBoundary.columns = lower + target + higher + highest
    _births, _deaths = zip(*phatBoundary.compute_persistence_pairs())
    births = set(_births)
    deaths = set(_deaths)
    
    essential = list(sorted(
        set(e for e in times-(births|deaths) if tranches[homology][0] <= e < tranches[homology][1])
    ))

    # Now, make sure we know when *all* the essential cycles are born.
    j = 0
    glb = len(lower)

    occupied = np.zeros((len(essential), len(target))).astype(int)

    # For each essential cycle born...
    for t in essential:
        # Determine which cells were included at the time the cycle was born, and
        # construct three assignments: the *occupied* cells, the *satisfied* cells,
        # and a new spin assignment on the *faces* of the cells.
        occupiedIndices = shuffled[:t-glb]
        occupied[j][occupiedIndices] = 1

        # Only sample the next cocycle from the time we homologically percolate,
        # not after.
        if (j+1) == stop:
            spins = sampleFromKernel(coboundary, F, includes=occupiedIndices)

        j += 1

    satisfied = np.zeros(len(target), dtype=int)
    satisfied[shuffled] = 1

    return spins, occupied, satisfied
