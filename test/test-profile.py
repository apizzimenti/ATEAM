
import cProfile
from potts import Lattice, SwendsonWang, Chain

usecProfile = False

if usecProfile:
    # Begin profiling.
    pr = cProfile.Profile()
    pr.enable()

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([20, 20], field=5)
model = SwendsonWang()
initial = model.initial(L)

# Create and run the chain.
chain = Chain(L, model, initial, steps=10)
for state in chain.progress(): pass

if usecProfile:
    # Disable profiling and write to file.
    pr.disable()
    print(f"{'x'.join([str(t) for t in L.corners])} integer lattice over GF({L.field.order})\n\n")
    pr.print_stats(sort="cumulative")
