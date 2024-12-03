
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
    # Choose *one* cube and, the entire time, ask what the sum of the coefficients
    # on its faces are.
    cube = model.lattice.cubes[13]
    q = model.lattice.field([state[model.lattice.index.faces[face]] for face in cube.faces]).sum()
    p = -1 if q else 1
    return p


def GraphWilsonLoop(model, state):
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
    # Choose *one* cube and, the entire time, ask what the sum of the coefficients
    # on its faces are.
    edge = model.lattice.cubes[13]
    u, v = edge.at
    return 1 if u.spin == v.spin else -1
