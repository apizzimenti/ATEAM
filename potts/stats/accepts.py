
import math
import numpy as np


def always(chain):
    """
    An acceptance function which always accepts the proposed state; we do nothing
    with the arguments to this function.

    Args:
        chain (Chain): Chain object containing the model and current state.

    Returns:
        A function which always returns True.
    """
    def _(proposed): return True
    return _


def metropolis(burnIn=1/100):
    def _metropolis(chain):
        def _(proposed):
            # First, check if we're within the burn-in window.
            wait = np.ceil(burnIn*chain.steps)
            if chain.step < wait: return True
            
            # Compute the energies of each state.
            oldEnergy = chain.model.energy(chain.lattice, chain.state)
            newEnergy = chain.model.energy(chain.lattice, proposed)

            # Compute the differences in the states' energy; if the energy has *increased*
            # (a thing we don't want), then our difference is less than 0, and we
            # navigate there with some small nonzero chance; if the energy is higher,
            # we always go there.
            diff = newEnergy-oldEnergy

            # Otherwise, if the difference is negative --- that is, the move is bad ---
            # we want to travel to this state with probability inversely exponential
            # to the difference (times our temperature parameter).
            p = min(math.exp(-chain.model.temperatureFunction(chain.step)*diff), 1)
            q = np.random.uniform()
            
            return q < p
        
        return _
    return _metropolis
