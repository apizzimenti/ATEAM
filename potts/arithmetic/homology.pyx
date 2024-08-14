
import numpy as np


def sampleFromKernel(A, field, includes=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if len(includes): B = A.take(includes, axis=0)
    else: B = A

    # TODO relative homology for this!
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = B.null_space()
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
        boundary,
        targetAbsoluteIndices,
        targetRelativeIndices,
        lowestRelativeIndex,
        highestRelativeIndex,
        lower,
        target,
        higher,
        other,
        shape,
        dimensions,
        times
    ):
    """
    Computes the persistent homology of the given complex, identifying when the
    first nontrivial element of the parent space (generally a torus) is born.

    Args:
        boundary (phat.boundary_matrix): The PHAT boundary matrix class; this
            object is instantiated once and re-used.
        targetAbsoluteIndices (Iterable): Absolute indices of the target cells.
        targetRelativeIndices (Iterable): Relative indices of the target cells.
        lowestRelativeIndex (int): Index of the first target cell.
        highestRelativeIndex (int): Index of the last target cell.
        lower (np.array): Lower-dimensional cells, listed in the PHAT boundary
            matrix format.
        target (np.array): Target-dimensional cells, listed in the PHAT boundary
            matrix format.
        higher (np.array): One-higher-than-target-dimensional cells, listed in
            the PHAT boundary matrix format.
        other (np.array): All higher-dimensional cells, listed in the PHAT
            boundary matrix format.
        shape (tuple): Shape of the re-indexed boundary matrix; should be
            `(<number of higher-dimensional cells>, -1)`, allowing NumPy to infer
            the other dimension.
        dimensions (Iterable): List containing the dimensions of the underlying
            lattice.
        times (set): Indices of the *filtration* (i.e. (0, ... , <number of cells>-1))
            reported as a set.

    Returns:
        A `np.array` of length `<number of target-dimensional cells>` with a `1`
        for each occupied cell and a `0` for each vacant cell.
    """
    # Randomize the indices of the target cells, then shuffle.
    shuffled = np.random.permutation(targetAbsoluteIndices)
    target = target[shuffled-lowestRelativeIndex]

    # Create a "reindexer" object which allows us to quickly determine the old->new
    # index mapping. Recall that the entries of the matrices *must* be sorted!
    _reindexer = np.array([shuffled, targetAbsoluteIndices])
    reindexer = _reindexer[:, _reindexer[0].argsort()][1]
    higher = (reindexer[higher-lowestRelativeIndex]).reshape(shape)
    higher.sort()

    # Now, construct the ordered boundary matrix for the complex.
    skeleta = [target, higher]
    zipped = [
        list(zip([d]*len(skeleton), skeleton))
        for d, skeleton in zip(dimensions, skeleta)
    ]
    columns = lower + sum(zipped, []) + other

    # Compute persistence pairs, and find the first time an essential cycle is
    # found.
    boundary.columns = columns
    _births, _deaths = zip(*boundary.compute_persistence_pairs())
    births = set(_births)
    deaths = set(_deaths)
    essential = (times-(births|deaths))-{0}
    essential = min(essential)

    # Find indices of all target cells included at or before the birth of the
    # first essential cycle.
    included = shuffled[:(essential-lowestRelativeIndex)+1]
    spins = np.zeros(len(target))
    spins[included-lowestRelativeIndex] = 1
    
    return spins

