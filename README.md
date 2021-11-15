# Stochastic cell rescaling

This repository contains additional information related to the stochastic cell rescaling barostat presented in Bernetti and Bussi, Pressure control using stochastic cell rescaling, [arXiv:2006.09250](https://arxiv.org/abs/2006.09250) (2020). It also contains materials associated to the anisotropic version of the barostat (manuscript currently in preparation). More materials associated to the anisotropic version will be added soon (work in progress).

## Reference implementations

The modified SimpleMD code is contained in [this directory](./simplemd).

The modified SimpleMD code with the anisotropic implementation is in [this directory](./simplemd_anisotropic).

The modified GROMACS code, implementing both the isotropic and the anisotropic versions, is contained in [this repository](https://github.com/bussilab/crescale-gromacs).

The modified LAMMPS code, implementing both the isotropic and the anisotropic versions, is contained in [this repository](https://github.com/bussilab/crescale-lammps).

## Input files

Input files for the isotropic version are available in [this directory](./input_file).

## Generated datasets

The generated datasets for the isotropic version are available on [Zenodo](https://doi.org/10.5281/zenodo.3921885).

## Analysis scripts

A notebook to perform the analysis reported in the manuscript related to the isotropic version is available [here](./Supporting_Info_figures.ipynb).
