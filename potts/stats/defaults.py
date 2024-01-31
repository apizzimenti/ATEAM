
from .energies import Hamiltonian as H


def Hamiltonian(model, state):
    """
    Wrapper method for computing the Hamiltonian of a particular state.

    Args:
        model (Model): Model on which we're conducting experiments.
        state (np.array): NumPy array.

    Returns:
        Hamiltonian energy for the given state.
    """
    return H(state, model.lattice.cubes, model.lattice.index.faces, model.lattice.field.order)
