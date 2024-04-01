
import numpy as np


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
