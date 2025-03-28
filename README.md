# Introduction

This script facilitates the convertion of MD trajectories (produced by GPUMD) into one unitcell structure.
One can see the Atomic distributions from MD as **NSR** [National Science Review 11 nwae216 (2024)].

# Input Files

Prepare two input files (inp.dat and dump.xyz):

1. **`inp.dat`**: Some basic information about the size of supercell, number of atomic type, and so on.
2. **`dump.xyz`**: MD trajectories.

An example, `inp.dat` and `dump.xyz` of Cs3Bi2I6Cl3, is provided with all necessary files.
