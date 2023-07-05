
from potts import Lattice, SwendsonWang, Chain
from potts.stats import critical

GL = Lattice([40, 40], field=7)
model = SwendsonWang(temperature=critical(GL.field.order))
initial = model.initial(GL) 

# Create and run the chain.
chain = Chain(GL, model, initial, steps=10)

for state in chain.progress(): pass
