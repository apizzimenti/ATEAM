
import numpy as np
from functools import reduce
from rustworkx import connected_components as components


def WilsonLoop(model, state):
    """
    Sums (over the finite field) the spins of faces of the "generalized loop"
    (which we're taking to mean the largest loop).

    Args:
        model (Model): Computational model we're using.
        state (np.array): Current State (assignment of spins) to faces of the
            Lattice of the Model.

    Returns:
        The evaluation of the cocycle (State) on the faces making up a connected
        component of the Lattice.
    """
    # Uniformly randomly choose a cube which is *not* occupied, and find the sum
    # of spins around its faces; in other words, we find a cycle which is not
    # a boundary and evaluate the cocycle on it.
    try:
        unoccupied = np.random.choice([c for c, i in model.occupied.items() if not bool(i)])
        index = model.lattice.index.cubes[unoccupied]
    except:
        unoccupied = None
    
    return 0 if not unoccupied else (model.lattice.coboundary[index]*state).sum()
