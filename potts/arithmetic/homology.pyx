
import numpy as np
from itertools import product


def sampleFromKernel(A, field, includes=[], relativeCells=None, relativeFaces=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if len(includes): B = A.take(includes, axis=0)
    else: B = A

    # # Relative homology.
    # if relativeCells: C = B.take(relativeCells, axis=0)
    # else: C = B

    # if relativeFaces: D = C.take(relativeFaces, axis=1)
    # else: D = C
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = B.null_space()
    M, _ = K.shape
    Q = field.Random(M)
    return (Q@K)


def evaluateCocycle(boundary, spins):
    evaluation = spins[boundary]

    for j in range(len(boundary[0])):
        if j%2: evaluation[j] = -evaluation[j]
    
    return evaluation.sum(axis=1)


def autocorrelation(data):
    """
    Computes the autocorrelation of a given observable over the provided time
    lag.

    Args:
        data (Iterable): An iterable, indexed by sample times, containing data
            from a given observable.
        lag (int): How far back do we look?
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


def essentialCyclesBorn(
        phatBoundary,
        coboundary,
        boundary,
        reindexed,
        tranches,
        homology,
        field,
        spins,
        times,
        indices,
    ):
    """
    Computes the persistent homology of the given complex, identifying when the
    first nontrivial element of the parent space (generally a torus) is born.

    Args:
    """
    # See which faces' spins sum to 0.
    cycles = evaluateCocycle(boundary[homology], spins)
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
    D = max(reindexed.keys())

    # Do all this dumb sorting.
    lower = sum([
        list(product([d], np.sort(reindexed[d]).tolist())) if d > 0 else [(0, [])]*len(reindexed[d])
        for d in range(homology)
    ], [])

    highest = sum([
        list(product([d], np.sort(reindexed[d]).tolist()))
        for d in range(homology+2, D+1)
    ], [])

    target = list(product(
        [homology], np.sort(boundary[homology][unordered]).tolist()
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
        if j < 1:
            spins = sampleFromKernel(coboundary, field, includes=occupiedIndices)

        j += 1

    satisfied = np.zeros(len(target), dtype=int)
    satisfied[shuffled] = 1

    return spins, occupied, satisfied
