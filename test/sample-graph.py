
from potts.models import GraphSwendsonWang
from potts.structures import GraphLattice
from potts.stats import critical
from potts import Chain

GL = GraphLattice([40, 40], field=7)
model = GraphSwendsonWang(temperatureFunction=critical(GL.field.order))
initial = model.initial(GL) 

# Create and run the chain.
chain = Chain(GL, model, initial, steps=10)

for state in chain.progress(): pass
