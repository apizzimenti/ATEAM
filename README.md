
## Installation
Clone this repository by `git clone https://github.com/apizzimenti/ateam.git` and run `sh setup.sh` in your favorite terminal. This installs the `ateam` package (and its dependencies) globally in development mode, so any changes you make to the source are reflected system-wide. **Please read the 

## Example Use

Below, we provide sample code for simulating a 10-iteration Markov chain determined by the plaquette invaded-cluster algorithm.

```python
from ateam.structures import Lattice
from ateam.models import InvadedCluster
from ateam import Chain

# Create a 64x64 square lattice with antipodal points identified (i.e. a torus),
# construct a model for the invaded cluster algorithm on the first homology
# group, and initialize a Markov chain simulating the algorithm for 10 steps.
L = Lattice().fromCorners([64,64], field=3)
HP = InvadedCluster(L, homology=1)
M = Chain(HP, steps=10)

for state in M:
    <do whatever>
```

## Citing

### BibTeX
```bibtex
@software{ATEAM,
    title={{ATEAM: Algebraic Topology-Enabled Algorithms for Mechanics}},
    author={Pizzimenti, Anthony E.},
    url={https://github.com/apizzimenti/ateam},
    version={1.0.0}
}
```
