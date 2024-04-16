
import numpy as np
from ..arithmetic import energy


def Hamiltonian(model):
    """
    Computes the typical Hamiltonian on a cocycle for the given lattice.

    Args:
        model (Model): Model defining the Markov chain.

    Returns:
        A closure which takes 
    """
    # Create blank array for fast summations.
    M = len(model.lattice.cubes[0].faces)
    N = len(model.lattice.cubes)

    # Create a big ol' numpy array for fast indexing; create another array of
    # just zeros.
    cubes = np.array([
        np.array([model.lattice.index.faces[f] for f in cube.faces])
        for cube in model.lattice.cubes
    ])

    Z = np.zeros(N, dtype=type(model.lattice.field.elements[1]))

    # Keep track of the indices of occupied edges.
    cubeIndices = set(range(N))

    def _(cocycle):
        return energy(cocycle, list(cubeIndices-set(model.occupied)), Z, cubes)
    return _
