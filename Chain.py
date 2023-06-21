
def always(lattice, state):
    """
    An acceptance function which always accepts the proposed state; we do nothing
    with the arguments to this function.

    Args:
        lattice (Lattice): The Lattice on which we're experimenting.
        state (np.array): A cocycle.
    """
    return True


class Chain:
    """
    Implements the Markov chain underlying the multidimensional Ising model.
    """

    def __init__(self, lattice, proposal, initial, accept=always, steps=10000):
        """
        Initializes the Chain object.

        Args:
            lattice (Lattice): A Lattice object.
            proposal (callable): A function which consumes the current state on
                the lattice and proposes a new state.
            accept (callable): A function which consumes the proposed next state
                and determines whether we travel to it.
            initial (np.array): A NumPy Array (homomorphism from the (k-1)-simplices
                of the lattice to the finite field over which the simplices are
                a vector space, i.e. a functional) assigning spins to simplices.
            steps (int): The number of iterations in the chain.
        """
        self.lattice = lattice
        self.proposal = proposal
        self.accept = accept
        self.initial = initial
        self.steps = steps


    def __iter__(self):
        """
        Initializes the Chain object as a generator.
        """
        self.step = 0
        self.state = self.initial
        return self
    

    def __next__(self):
        """
        Performs the computations specified by the proposal and acceptance schemes.
        """
        # Kick off the iterations!
        if self.step == 0:
            self.step += 1
            return self.state
        
        # While we haven't reached the max number of steps, propose a new plan,
        # check whether it's acceptable/valid, and continue.
        while self.step < self.steps:
            proposed = self.proposal(self.lattice, self.state)
            self.state = (proposed if self.accept(self.lattice, proposed) else self.state)
            self.step += 1
        
        # Done!
        raise StopIteration
    
    def progress(self):
        """
        Progress bar.
        """
        from tqdm.auto import tqdm
        return tqdm(self)
            