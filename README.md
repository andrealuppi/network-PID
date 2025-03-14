# network-PID: Quantifying synergy and redundancy between networks
Authors: A.I. Luppi, E. Olbrich, C. Finn, L.E. Suárez, F.E. Rosas, P.A.M. Mediano, J. Jost

This repository provides code to illustrate the central method in Luppi et al. (2024) "Quantifying synergy and redundancy between networks" ([_Cell Reports Physical Science_](10.1016/j.xcrp.2024.101892)).

## Code description
The code consists of the MATLAB function `PartialNetworkDecomposition_GlobEff.m`
Given two binary undirected networks NETX and NETY, this code computes the partial network decomposition of global efficiency (i.e. inverse shortest path length). It returns a structure with the redundant, unique, and synergistic contributions corresponding to both networks.
Partial Network Decomposition is inspired by the PArtial Information Decomposition of Williams and Beer (2010) ([arXiv](
https://doi.org/10.48550/arXiv.1004.2515))

## Contact Information
This code was developed in MATLAB 2019a by Andrea Luppi and Pedro Mediano.
The code relies on MATLAB code from the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet) for MATLAB by Rubinov and Sporns (2010) NeuroImage.
For questions, please email: [al857@cam.ac.uk](al857@cam.ac.uk).

