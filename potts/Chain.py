
from .stats import always

class Chain:
    """
    Implements the Markov chain underlying the multidimensional Ising model.
    """

    def __init__(
            self, model, initial=None, accept=always, sampleInterval=0, statistics={},
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
        self.accept = accept(self)

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
            # Propose the next state and check whether it's valid; assign the
            # state to the Model.
            proposed = self.model.proposal(self.step)
            self.state = (proposed if self.accept(proposed) else self.state)
            self.model.assign(self.state)

            # Compute statistics.
            for name, function in self.functions.items():
                self.statistics[name].append(function(self.model, self.state))

            # If we're collecting samples, collect!
            # try:
            #     if self.step % self.sampleInterval == 0:
            #         self.assignments.append(list(self.state))
            # except: pass
            
            # Iterate.
            self.step += 1
            
            return self.state
        
        # If we haven't returned, we're done; convert the assignments to JSON-ifiable
        # types, and stop iteration.
        # self.assignments = [[int(s) for s in assignment] for assignment in self.assignments]
        
        raise StopIteration
    

    def progress(self):
        """
        Progress bar.
        """
        from tqdm.auto import tqdm
        return tqdm(self, total=self.steps)
            