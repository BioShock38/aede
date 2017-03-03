# Aede

Local Ancestry Simulation

# Tutorial

You need an array of distance in morgan called a genetic map.
So for each SNP you have the distance in morgan in `gen_map`.

The most simple way to run a simulation, with default parameters
```python
import aede.simulator as simu

loc_anc = simu.LocalAncestrySimulator(gen_map, simu.default_pipeline_hapmix)
```

Then we can build admixed haplotypes with parental haplotypes
```python
h_adm, jumps, clusters = loc_anc.run(pop_1=H_ceu, pop_2=H_yri)
```
We obtain `h_adm` a numpy array of admixed haplotypes, `jumps` is a boolean array with True where the jumps are.
`clusters` is an array indicating from which population is the SNP selected from.

# Full Example

```python
import numpy as np
from functools import partial

H_ceu = np.load("H_ceu.npy")
H_yri = np.load("H_yri.npy")

gen_map = np.load("gen_map.npy")

import aede.simulator as simu

# Custom lambda parameter for exponential law
jumps_10 = partial(jumps_builder, lambd=10)

loc_anc = simu.LocalAncestrySimulator(gen_map, simu.pipeline_hapmix)
h_adm, jumps, clusters = loc_anc.run(pop_1=H_ceu, pop_2=H_yri,
                                     jumps_builder=jumps_10,
                                     cluster_builder=simu.cluster_default) 
```

# Simulation Design

To customize the LocalAncestrySimulator one need to define a pipeline taking at least a 
genetic map in input.
When you run the simulator every parameter of the `run` will be pass to the pipeline you wrote.
