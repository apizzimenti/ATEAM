
from potts.models import GraphIsing
from potts.structures import GraphLattice
from potts.stats import critical, metropolis
from potts import Chain

GL = GraphLattice([40, 40], field=2)
model = GraphIsing(temperatureFunction=critical(GL.field.order))
initial = model.initial(GL) 

# Create and run the chain.
chain = Chain(GL, model, initial, steps=10, accept=metropolis)

for state in chain.progress(): pass
