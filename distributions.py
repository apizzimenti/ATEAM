
import numpy as np

def uniform(low, high):
    """
    Discrete uniform distribution on [low, high).

    Args:
        low (int): Lower end of interval.
        high (int): Upper end of interval.
    
    Returns:
        A value selected uniformly randomly from the interval [low, high).
    """
    return np.random.randint(low, high=high)
