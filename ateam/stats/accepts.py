
import math
import numpy as np


def always(model=None, distribution=None, burnIn=None):
    """
    An acceptance function which always accepts the proposed state; we do nothing
    with the arguments to this function.

    Args:
        model (Model): Optional argument; does nothing.
        distribution (Callable): Optional argument; does nothing.
        burnIn (int): Optional argument; does nothing.

    Returns:
        A function which always returns True.
    """
    def _(current, proposed, time): return True
    return _


def MetropolisHastings(model, distribution, burnIn=0):
    """
    Args:
        model (Model): Evolution model we're using.
        distribution (Callable): Distribution function which consumes a cocycle
            (a list of assignments to faces).
        burnIn (int=0): How long we wait to begin the Metropolis-Hastings scheme.

    Returns:
        A closure which decides whether the proposed state is accepted or not.
    """
    def _(current, proposed, time) -> bool:
        """
        Args:
            current (np.array): Array of current spin assignments.
            proposed (np.array): Array of proposed spin assignments.
            time (int): Current step of the chain.

        Returns:
            A decision on whether the proposed state is accepted or not.
        """
        # First, check if we're within the burn-in window.
        if time < burnIn: return True
        
        # Compute the energies of each state.
        oldEnergy = distribution(current)
        newEnergy = distribution(proposed)

        # Compute the differences in the states' energy; if the energy has *increased*
        # (a thing we don't want), then our difference is less than 0, and we
        # navigate there with some small nonzero chance; if the energy is higher,
        # we always go there.
        diff = newEnergy-oldEnergy

        # Otherwise, if the difference is negative --- that is, the move is bad ---
        # we want to travel to this state with probability inversely exponential
        # to the difference (times our temperature parameter).
        p = min(math.exp(model.temperatureFunction(time)*diff), 1)
        q = np.random.uniform()
        
        return q < p
    return _
