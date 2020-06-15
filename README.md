# Implementations of stochastic cell rescaling

This repository contains the reference implementations for the stochastic cell rescaling barostat presented in Bernetti and Bussi, in preparation.

The modified SimpleMD code is contained in [this directory](./simplemd).

The modified GROMACS code is contained in [this repository](https://github.com/bussilab/crescale-gromacs).

### Notes on the GROMACS implementation

Several branches and tags can be found:
- Tag `v2019.4-crescale` is the version used in the manuscript and is based on GROMACS 2019.4.
- Tag `v2019.6-crescale ` integrates GROMACS updates up to 2019.6.
- Tag `v2020.2-crescale` is based on GROMACS 2020.2.
- Branch `release-2019-crescale` is based on the `release-2019` branch of the official GROMACS repository.
- Branch `release-2020-crescale` is based on the `release-2020` branch of the official GROMACS repository.
- Branch `master-crescale` is based on the `master` branch of the official GROMACS repository.

The changes in version 2019 are [relatively small](https://github.com/bussilab/crescale-gromacs/compare/328a18d71dda42ca67edf76b2f93781dab6fdf9d..43b914078bf7ee83afe17f53e559e02932af7ae6). The added code is basically a copy of the Berendsen code, with some modification. The changes in version 2020 are [slighly more complicated](https://github.com/bussilab/crescale-gromacs/compare/5e788350ad75c15ba91d2ba02779f1f8200f61ee..0fb489c52ca1732f6dd2a26ce988e2b630796965) since in this version the coordinate update might be done on the GPU only. In this case, it was necessary to add a function to allow scaling of the velocities to be performed on the GPU.
