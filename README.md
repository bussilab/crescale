# Implementations of stochastic cell rescaling

This repository contains the reference implementations for the stochastic cell rescaling barostat presented in Bernetti and Bussi, Pressure control using stochastic cell rescaling, [arXiv:2006.09250](https://arxiv.org/abs/2006.09250) (2020).

The modified SimpleMD code is contained in [this directory](./simplemd).

The modified GROMACS code is contained in [this repository](https://github.com/bussilab/crescale-gromacs).

### Notes on the GROMACS implementation

Several branches and tags can be found:
- Tag `v2019.4-crescale` is the version used in the manuscript and is based on GROMACS 2019.4.
- Tag `v2019.6-crescale ` integrates GROMACS updates up to 2019.6.
- Tag `v2020.2-crescale` is based on GROMACS 2020.2.
- Branches `release-2019`, `release-2020`, and `master` are in sync with the official GROMACS repository.
- Branch `release-2019-crescale` is based on branch `release-2019`.
- Branch `release-2020-crescale` is based on branch `release-2020`.
- Branch `master-crescale` is based on branch `master`.

The changes in version 2019 are [relatively small](https://github.com/bussilab/crescale-gromacs/compare/release-2019..release-2019-crescale). The added code is basically a copy of the Berendsen code, with some modification. The changes in version 2020 are [slighly more complicated](https://github.com/bussilab/crescale-gromacs/compare/release-2020..release-2020-crescale) since in this version the coordinate update might be done on the GPU only. In this case, it was necessary to add a function to allow scaling of the velocities to be performed on the GPU. Changes on master branch are [similar](https://github.com/bussilab/crescale-gromacs/compare/master..master-crescale) and just required some adjustment to be compatible with master branch.
