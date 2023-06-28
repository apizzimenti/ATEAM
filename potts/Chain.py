
class Chain:
    """
    Implements the Markov chain underlying the multidimensional Ising model.
    """

    def __init__(self, lattice, model, initial, statistics={}, steps=10000):
        """
        Initializes the Chain object.

        Args:
            lattice (Lattice): A Lattice object.
            proposal (callable): A function which consumes this Chain object and
                proposes a new state.
            accept (callable): A function which consumes this Chain object and
                determines whether we travel to it.
            initial (np.array): A NumPy Array (homomorphism from the (k-1)-simplices
                of the lattice to the finite field over which the simplices are
                a vector space, i.e. a functional) assigning spins to simplices.
            statistics (dict): A mapping of names to functions which take the lattice
                as an argument. The Chain keeps track of these at each iteration
                and stores whatever output is given.
            steps (int): The number of iterations in the chain.
        """
        self.lattice = lattice
        self.model = model
        self.initial = initial
        self.steps = steps

        # Store stats and things.
        self.functions = statistics
        self.functions["energy"] = model.energy
        self.statistics = { name: [] for name in self.functions.keys() }


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
        # While we haven't reached the max number of steps, propose a new plan,
        # check whether it's acceptable/valid, and continue.
        while self.step < self.steps:
            # Propose the next state and check whether it's valid.
            proposed = self.model.proposal(self)
            self.state = (proposed if self.model.accept(self) else self.state)

            # Assign things to the lattice and collect statistics.
            self.lattice.assign(self.state)
            for name, function in self.functions.items(): self.statistics[name].append(function(self.lattice))
            
            # Iterate.
            self.step += 1
            
            return self.state
        
        # Done!
        raise StopIteration
    

    def progress(self):
        """
        Progress bar.
        """
        from tqdm.auto import tqdm
        return tqdm(self, total=self.steps)
            