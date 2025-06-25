# Electro-magneto-elasticity in the eddy current approximation

This repository contains the code base for the FEniCSx finite-element implementation of the large deformation theory we developed for electro-magneto-elasticity in the eddy current approximation. 

This repository contains the relevant FEniCSx Jupyter notebooks,  FEniCSx Python script files, mesh files, and experimental data files which were used in simulating the benchmark TEAM problems 12 and 28, as reported in Part II of the associated two-part paper: 

- L. Anand and E. M. Stewart. Electro-magneto-elasticity in the eddy current approximation. Part I: Theory. *Journal of the Mechanics and Physics of Solids*, Submitted.

- E. M. Stewart and L. Anand. Electro-magneto-elasticity in the eddy current approximation. Part II: Applications. *Journal of the Mechanics and Physics of Solids*, Submitted.

# Running the codes

These codes can be run using either FEniCSx v0.9.0 or v0.8.0. FEniCSx can be installed somewhat simply in a Conda environment as described on the FEniCS project website [here](https://fenicsproject.org/download/), following which the codes can be run in that environment using VSCode. 

Alternatively, we also provide a detailed guide for installing FEniCSx v0.8.0 in a Docker container and running the notebooks using VSCode, both for [Mac](https://github.com/ericstewart36/finite_viscoelasticity/blob/main/FEniCSx_v08_Docker_install_mac.pdf) and [Windows](https://github.com/ericstewart36/finite_viscoelasticity/blob/main/FEniCSx_v08_Docker_install_windows.pdf). 

These are our preferred methods for editing and running FEniCSx codes, although [many other options exist](https://fenicsproject.org/download/).

We have provided a Python script version of the TEAM12 simulation which is meant to be run using MPI parallelization. To run this script in parallel on e.g. 16 cores use the following command syntax in the terminal:  

```
mpirun -n 16 python3 TEAM12_3D_eddy_bending_bx0p5T.py
```

![](https://github.com/ericstewart36/eddy_current_elasticity/blob/main/example_movies.gif)

# Citations

- L. Anand and E. M. Stewart. Electro-magneto-elasticity in the eddy current approximation. Part I: Theory. *Journal of the Mechanics and Physics of Solids*, Submitted.

- E. M. Stewart and L. Anand. Electro-magneto-elasticity in the eddy current approximation. Part II: Applications. *Journal of the Mechanics and Physics of Solids*, Submitted.
