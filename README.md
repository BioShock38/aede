# Aede

Local Ancestry Simulation

# Install

From github
```bash
git clone https://github.com/BioShock38/aede.git

cd aede

python setup.py install --user
```


# Tutorial

You need an array of distances in morgan also called a genetic map.
So for each SNP you have the distance in morgan in `gen_map`.

The most simple way to run a simulation, with default parameters
```python
import aede.simulator as simu

loc_anc = simu.LocalAncestrySimulator(gen_map, simu.default_pipeline_hapmix)
```

Then we can build admixed haplotypes with parental haplotypes
```python
h_adm, clusters = loc_anc.run(pop_1=H_ceu, pop_2=H_yri)
```
We obtain `h_adm` a numpy array of admixed haplotypes and `clusters` is an array indicating from which population is the SNP selected from.

You can wrap the functions, to build one function.
```python
import aede.simulator as simu
from functools import partial

def h_adm_simu(lambd, H_1, H_2, genmap, mu=0.8, sigma2=0.01):
    jumps_fun = partial(simu.jumps_builder, lambd=lambd)
    cluster_fun = partial(simu.cluster_builder, a_beta=mu, b_beta=sigma2)

    loc_anc = simu.LocalAncestrySimulator(genmap, simu.pipeline_hapmix)
    h_adm, clusters = loc_anc.run(pop_1=H_1, pop_2=H_2,
                                  jumps_builder=jumps_fun,
                                  cluster_builder=cluster_fun)
    return h_adm, clusters"
```
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
h_adm, clusters = loc_anc.run(pop_1=H_ceu, pop_2=H_yri,
                              jumps_builder=jumps_10,
                              cluster_builder=simu.cluster_default) 
```
