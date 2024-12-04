
from .stats import always


class Chain:
    def __init__(self, model, accept=always(), statistics={}, steps=10000):
        """
        Simulates a Markov chain on the given Model.

        Args:
            model (Model): An instantiated descendant of `Model` (e.g. `HomologicalPercolation`)
                generating the Markov chain.
            accept (Callable): A function which consumes the lattice, model, and
                state to determine whether we're going to a good place.
            statistics (dict): A mapping of names to functions which take the lattice
                as an argument. The Chain keeps track of these at each iteration
                and stores whatever output is given.
            steps (int): The number of iterations in the chain.
        """
        self.model = model
        self.steps = steps
        self.accept = accept

        # Store stats and things.
        self.functions = statistics
        self.statistics = { name: [] for name in self.functions.keys() }


    def __iter__(self):
        """
        Initializes the Chain object as a generator.
        """
        self.step = 0
        self.state = self.model.spins
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
            proposed = self.model.proposal(self.step)
            self.state = proposed if self.accept(self.state, proposed, self.step) else self.state
            self.model.assign(self.state[0])

            # Compute statistics.
            for name, function in self.functions.items():
                self.statistics[name].append(function(self.model, self.state))
            
            # Iterate.
            self.step += 1
            
            return self.state
        
        raise StopIteration
    
    
    def progress(self):
        """
        Progress bar.

        Returns:
            `tqdm` iterable.
        """
        from tqdm.auto import tqdm
        return tqdm(self, total=self.steps)
