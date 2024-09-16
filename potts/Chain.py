
from .stats import always


class Chain:
    """
    Simulates a Markov chain on the given Model.
    """

    def __init__(
            self, model, initial=None, accept=always(), sampleInterval=0, statistics={},
            steps=10000
        ):
        """
        Initializes the Chain object.

        Args:
            proposal (callable): A function which consumes this Chain object and
                proposes a new state.
            initial (np.array): A NumPy Array (homomorphism from the (k-1)-simplices
                of the lattice to the finite field over which the simplices are
                a vector space, i.e. a functional) assigning spins to simplices.
            accept (callable): A function which consumes the lattice, model, and
                state to determine whether we're going to a good place.
            sampleInterval (int): Number representing the number of spin assignments we
                should save.
            statistics (dict): A mapping of names to functions which take the lattice
                as an argument. The Chain keeps track of these at each iteration
                and stores whatever output is given.
            steps (int): The number of iterations in the chain.
        """
        self.model = model
        self.initial = initial if initial else model.initial()
        self.steps = steps
        self.accept = accept

        # Store stats and things.
        self.functions = statistics
        self.statistics = { name: [] for name in self.functions.keys() }

        # Assignment-related things.
        self.sampleInterval = sampleInterval
        self.assignments = []


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
            # Propose the next state and check whether we want to accept it as
            # the next state or not; assign whichever state is chosen to the
            # Model.
            proposedStates = self.model.proposal(self.step)
            proposed = proposedStates[0]
            self.state = proposed if self.accept(self.state, proposed, self.step) else self.state
            self.model.assign(self.state)

            # Compute statistics.
            for name, function in self.functions.items():
                self.statistics[name].append(function(self.model, self.state))
            
            # Iterate.
            self.step += 1
            
            return proposedStates
        
        raise StopIteration
    
    
    def progress(self):
        """
        Progress bar.
        """
        from tqdm.auto import tqdm
        return tqdm(self, total=self.steps)
