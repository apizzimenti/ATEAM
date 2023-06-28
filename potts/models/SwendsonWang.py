
import numpy as np
import pandas as pd
import numba as nb
from pathlib import Path

from ..stats import constant, always, uniform
from .Model import Model


class SwendsonWang(Model):
    def __init__(
            self, temperature=constant(0.6), accept=always(), testing=False
        ):
        """
        Initializes a Swendson-Wang evolution on the Potts model.

        Args:
            temperature (Callable): A function taking an integer argument which
                returns a real-valued temperature. Defaults to the constant
                schedule.
            accept (Callable): A function which takes a Chain as input and returns
                a boolean.
        """
        self.temperature = temperature
        self.accept = accept
        self.testing = testing
        self.log = ""

        self._directorySetup()


    def _directorySetup(self):
        if self.testing:
            # Creates an output directory if none exists.
            outroot = Path("./output/")
            if not outroot.exists():
                outroot.mkdir()
                (outroot/"figures").mkdir()
                (outroot/"matrices").mkdir()


    def _proposal_multiply(self, chain):
        """
        Multiplies the computed boundary matrix B on the left by an edge-inclusion
        matrix E. As the number of edges grows, this matrix multiplication is
        expected to take a really long time. DEPRECATED.
        """
        lattice, state = chain.lattice, chain.state

        # First, create a matrix that includes/excludes the appropriate edges;
        # we scan through them, and "remove" any edges that have no effect.
        M = len(lattice.structure[1])
        included = lattice.field(np.zeros((M, M), dtype=int))

        for edge in lattice.structure[1]:
            u, v = edge.coordinates
            
            # We have three cases:
            #
            #   1. if the spins on u and v are different, we ignore the edge.
            #   2. if the spins on u and v are the same, we
            #       a. include the edge with probability 1-e^beta,
            #       b. ignore the edge with probability e^beta.

            # Case 1; re-set edge spin to 0.
            if state[u.index] != state[v.index]:
                edge.spin = 0
                continue
                
            # Case 2.
            else:
                # Set the probability that an edge is included, and uniformly
                # randomly select a value in [0,1] to determine whether the edge
                # is included.
                p = 1-(np.exp(self.temperature(chain.step)))
                q = np.random.uniform()

                # Case 2(a).
                if q < p:
                    # Testing!
                    if self.testing: self.log += f"included edge ({u},{v}) with probability {p} > {q}\n"
                    
                    # Set the spins of the edges appropriately.
                    included[edge.index, edge.index] = 1
                    edge.spin = 1
                
                # Case 2(b).
                else:
                    edge.spin = 0

        # Now that we know which edges we're ignoring, we can multiply the two
        # matrices together!
        return np.matmul(lattice.boundary, included)
    
    @staticmethod
    def _excludeEdges(state, probability):
        """
        Given a vector of edges, returns a vector of the same shape with values
        in {0,1} indicating whether each edge is included or excluded.
        """
        # @nb.vectorize()
        def _(edge):
            # Obtain the spins on the vertices incident to the edge.
            u = edge.coordinates[0]
            v = edge.coordinates[1]

            if state[u.index] != state[v.index]:
                edge.spin = 0
                return 0
            else:
                q = np.random.uniform()

                if q < probability:
                    edge.spin = 1
                    return 1
                else:
                    edge.spin = 0
                    return 0

        return _


    def _proposal_copy(self, chain):
        """
        Copies the boundary matrix and zeros out rows instead of multiplying
        large matrices.
        """
        lattice, state = chain.lattice, chain.state

        # Copy the boundary matrix for editing; we want to make *one* call that
        # determines which edges are included and which aren't.
        M = lattice.boundary.copy()
        exclusion = self._excludeEdges(state, 1-(np.exp(self.temperature(chain.step))))
        excluded = np.array(list(map(exclusion, lattice.structure[1])))
        exclude = np.where(excluded==0)
        M[:,exclude] = 0

        return M


    def proposal(self, chain):
        """
        Proposal scheme for the Swendson-Wang evolution on the Potts model. 

        Args:
            chain (Chain): Chain object which contains all the information we
                need.

        Returns:
            A proposed state.
        """
        # If we're testing, we want to record probabilities; since this is called
        # at each step, we'll write something to file at each step.
        if self.testing: self.log += f"### STEP {chain.step} ###\n"
        # boundary = self._proposal_multiply(chain)
        boundary = self._proposal_copy(chain)

        # If we're testing, write the matrices to file and indicate which edges
        # are ignored.
        if self.testing:
            pd.DataFrame(chain.state).to_csv(f"./output/matrices/state-{chain.step}.csv", index=False)
            pd.DataFrame(boundary).to_csv(f"./output/matrices/operator-{chain.step}.csv", index=False)
            self.log += "\n\n"
    
        # Now, let's change the states! First, take the transpose of the boundary
        # matrix (i.e. the coboundary matrix). Then, find a basis for the null
        # space of this matrix, which is equivalent to finding a basis for the
        # vector space of cocycles. Then, (uniformly randomly) sample coefficients
        # from the field underlying the vector space, and scale the ith column of
        # the basis matrix (the ith basis vector) by the ith coefficient. Finally,
        # sum the *rows* of the linear combination, equivalent to summing the
        # vectors element-wise to produce a cocycle!
        coboundary = boundary.T
        basis = coboundary.null_space()
        coefficients = np.random.randint(0, chain.lattice.field.order, size=len(basis))
        linearCombination = (basis.T * coefficients).T
        solution = list(np.add.reduce(linearCombination))

        return solution
    

    def initial(self, lattice, distribution=uniform):
        """
        Generates a random initial state.

        Args:
            lattice (Lattice): The integer lattice we're experimenting on.
            distribution (Callable): A distribution from which we'll select;
                this distribution must be _discrete_ and take _two integer inputs_
                representing the lower and higher endpoints of the interval over
                the integers.

        Returns:
            A cocycle with initial states.
        """
        vertices = lattice.structure[0]
        cocycle = np.array([
            distribution(0, lattice.field.order) for _ in range(len(vertices))
        ])

        # If we're testing, give all the vertices the same spin. This way, we
        # should ignore *no* edges.
        if self.testing: return lattice.field(np.array([1]*len(vertices)))
        return lattice.field(cocycle.astype(int))


    def energy(self, lattice):
        """
        Computes the Hamiltonian (energy) of the lattice in its current state.

        Args:
            lattice (Lattice): Lattice over which we're working.

        Returns:
            Integer representing the Hamiltonian.
        """
        s = 0

        for edge in lattice.structure[1]:
            u, v = edge.coordinates
            s += (1 if u.spin == v.spin else 0)

        return -s

