
from potts.models import GraphPercolation
from potts.structures import GraphLattice
from potts.stats import critical
from potts import Chain
import pandas as pd

GL = GraphLattice([20,20], field=2)
model = GraphPercolation(temperatureFunction=critical(GL.field.order))
initial = model.initial(GL) 

# Create and run the chain.
chain = Chain(GL, model, initial, steps=2)

for state in chain.progress(): pass
