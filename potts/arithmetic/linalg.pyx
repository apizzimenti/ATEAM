
import numpy as np

def sampleFromKernel(A, field, includes=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if includes: B = A.take(includes, axis=0)
    else: B = A
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = B.null_space()
    M, _ = K.shape
    Q = field.Random(M, dtype=np.int32)
    return (Q@K).astype(np.int32)


def autocorrelation(data, lag=500):
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
        np.dot(normalized[:N-t], normalized[t:])*(1/N) for t in range(N-lag)
    ])

    return autocorrs/autocorrs[0]

