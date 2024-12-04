
import numpy as np

def energy(cochain, unoccupied, Z, cubes):
    """
    Computes the Hamiltonian of the spin assignment (cochain).

    Args:
        cochain (galois.FieldArray): FieldArray of spin assignments.
        unoccupied (set): Set of unoccupied indices.
        Z (np.array): Array of Galois arrays to do field arithmetic quickly.
        cubes (np.ndarray): Matrix where the ith row contains the indices of the
            faces of the ith cube.

    Returns:
        A count of plaquettes on which the cochain does not evaluate to zero.
    """
    # Get the sums of coefficients in the field, then count the nonzero entries.
    for u in unoccupied: Z[u] = cochain[cubes[u]].sum()
    return -np.count_nonzero(Z[unoccupied])