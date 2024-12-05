
## Installation
Clone this repository by `git clone https://github.com/apizzimenti/ateam.git`, navigate into the project's root, and run `sh setup.sh` in your favorite terminal. This installs the `ateam` package (and its dependencies) globally in development mode, so any changes you make to the source are reflected system-wide. **Please read the installation note in `setup.sh` to ensure PHAT is correctly installed on your system.**

## Example Use

This sample code runs the plaquette invaded-cluster algorithm on a 64x64 cubical 2-torus.

```python
from ateam.structures import Lattice
from ateam.models import InvadedCluster
from ateam import Chain

L = Lattice().fromCorners([64,64], field=3)
HP = InvadedCluster(L, homology=1)

for state in Chain(HP, steps=10):
    <do whatever>
```

## Citing

### BibTeX
```bibtex
@software{ATEAM,
    title={{ATEAM: Algebraic Topology-Enabled Algorithms for Mechanics}},
    author={Pizzimenti, Anthony E.},
    url={https://github.com/apizzimenti/ateam},
    version={1.0.0},
    doi={10.5281/zenodo.14284172}
}
```
