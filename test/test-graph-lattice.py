
from potts.models import GraphIsing
from potts.structures import GraphLattice
from potts.stats import critical
from potts import Chain
import pandas as pd
from pympler import asizeof

GL = GraphLattice([512,512], field=2)
print(asizeof.asizeof(GL))
model = GraphIsing(temperatureFunction=critical(GL.field.order))
initial = model.initial(GL) 

# Create and run the chain.
chain = Chain(GL, model, initial, steps=10000)

for state in chain.progress(): pass
