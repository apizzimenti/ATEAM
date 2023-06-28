
from potts import Lattice, SwendsonWang, Chain

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([20, 20], field=5)
model = SwendsonWang()
initial = model.initial(L)

# Create and run the chain.
chain = Chain(L, model, initial, steps=1000)

for state in chain.progress(): pass
