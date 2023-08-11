
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


def metropolis(chain):
    def _(proposed):
        # Compute the energies of each state.
        oldEnergy = chain.model.energy(chain.lattice, chain.state)
        newEnergy = chain.model.energy(chain.lattice, proposed)

        # Compute the differences in the states' temperatures.
        diff = newEnergy-oldEnergy
        p = min(np.exp(chain.model.temperatureFunction(chain.step) * diff), 1)
        q = np.random.uniform()

        if newEnergy < oldEnergy:
            return True
        else:
            if q < p: return True
        
        return False
    
    return _
