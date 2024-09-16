
import numpy as np


def sampleFromKernel(A, field, includes=None, relativeCells=None, relativeFaces=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if len(includes): B = A.take(includes, axis=0)
    else: B = A

    # Relative homology.
    if len(relativeCells): C = B.take(relativeCells, axis=0)
    else: C = B

    if len(relativeFaces): D = C.take(relativeFaces, axis=1)
    else: D = C
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = D.null_space()
    M, _ = K.shape
    Q = field.Random(M)
    return (Q@K)


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
        tranches,
        skeleta,
        homology,
        field,
        spins,
        times,
        indices
    ):
    """
    Computes the persistent homology of the given complex, identifying when the
    first nontrivial element of the parent space (generally a torus) is born.

    Args:
    """
    # Randomize the indices of the target cells, and check that the spins on
    # their faces sum to 0.
    shuffled = np.random.permutation(skeleta[homology])
    shuffledIndices = shuffled-tranches[homology-1]
    adj = (tranches[homology-2] if homology-2>=0 else 0)

    faces = boundary[homology][shuffledIndices]-adj
    cycles = spins[faces].sum(axis=1)
    satisfiedIndices = (cycles==0).nonzero()[0]
    satisfied = faces[satisfiedIndices]+adj
    satisfied.sort()

    # Add in the "unsatisfied" bonds afterward; note the last time a satisfied
    # bond was added.
    unsatisfiedIndices = np.setdiff1d(shuffledIndices,satisfiedIndices)
    unsatisfied = faces[unsatisfiedIndices]+adj
    unsatisfied.sort()

    _reindexer = np.array([
        np.append(shuffled[satisfiedIndices], shuffled[unsatisfiedIndices]),
        indices
    ])

    reindexer = _reindexer[:, _reindexer[0].argsort()][1]

    higher = (boundary[homology+1]-tranches[homology-1]).flatten()
    higher = (reindexer[higher]).reshape((len(skeleta[homology+1]), -1))
    higher.sort()

    # Construct the columns of the boundary matrix.
    lower = sum([list(zip([d]*len(skeleta[d]), boundary[d])) for d in range(homology)], [])
    target = list(zip([homology]*len(satisfied), satisfied)) + list(zip([homology]*len(unsatisfied), unsatisfied))
    reindexed = list(zip([homology+1]*len(higher), higher))
    highest = sum([list(zip([d]*len(skeleta[d]), boundary[d])) for d in range(homology+2, len(skeleta))], [])

    # Compute persistence pairs, and find the first time an essential cycle is
    # found.
    phatBoundary.columns = lower + target + reindexed + highest
    _births, _deaths = zip(*phatBoundary.compute_persistence_pairs())
    births = set(_births)
    deaths = set(_deaths)
    essential = set(e for e in times-(births|deaths) if tranches[homology-1] < e < tranches[homology])
    essential = min(essential)
    
    # Determine which cells were included at the time the cycle was born, and
    # construct three assignments: the *occupied* cells, the *satisfied* cells,
    # and a new spin assignment on the *faces* of the cells.
    occupiedIndices = shuffled[:(essential-tranches[homology-1])+1]-tranches[homology-1]
    occupied = np.zeros(len(target), dtype=int)
    occupied[occupiedIndices] = 1

    satisfiedIndices = shuffled[satisfiedIndices]-tranches[homology-1]
    satisfied = np.zeros(len(target), dtype=int)
    satisfied[satisfiedIndices] = 1

    # Construct a new assignment on the faces.
    spins = sampleFromKernel(coboundary, field, occupiedIndices)

    return spins, occupied, satisfied
