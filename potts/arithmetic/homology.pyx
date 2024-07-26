
import numpy as np


def essentialCyclesBorn(boundary, targetAbsoluteIndices, targetRelativeIndices, lowestRelativeIndex, highestRelativeIndex, lower, target, higher, other, shape, dimensions, times):
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
    skeleta = [lower, target, higher, other]
    zipped = [
        list(zip([d]*len(skeleton), skeleton))
        for d, skeleton in zip(dimensions, skeleta)
    ]
    columns = sum(zipped, [])

    # Compute persistence pairs, and find the first time an essential cycle is
    # found.
    boundary.columns = columns
    _births, _deaths = zip(*boundary.compute_persistence_pairs())
    births = set(_births)
    deaths = set(_deaths)
    essential = times-(births|deaths)-{0}
    essential = min(essential)

    # Find indices of all target cells included at or before the birth of the
    # first essential cycle.
    included = shuffled[:(essential-lowestRelativeIndex)+1]
    spins = np.zeros(len(target))
    spins[included-lowestRelativeIndex] = 1
    
    return spins

