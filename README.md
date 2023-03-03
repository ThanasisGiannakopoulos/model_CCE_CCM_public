# Toy models for Cauchy-Characteristic Extraction (CCE) and Matching (CCM)
 
This repository contains the code used in the paper [arxiv XXX.XXXXX](...).
The goal of this work is to analyze the numerical convergence of CCE and CCM
for toy models that mimic the hyperbolic structure of the general relativistic PDE systems
used in CCE and CCM. The implemented models are the following:

- for the Initial Boundary Value Problem (IBVP):
$$\partial_t \phi_1 = - v_{\phi_1} \partial_\rho \phi_1 + a_{z11} \partial_z \phi_1 + a_{z12} \partial_z \psi_{v 1} + a_{z13} \partial_z \psi_1 + b_{11} \phi_1 + b_{12} \psi_{v 1} + b_{13}\psi_1 $$ 

$$\partial_t \psi_{v1} = - v_{\psi_{v1}} \partial_\rho \psi_{v1} + a_{z21} \partial_z \phi_1 + a_{z22} \partial_z \psi_{v 1} + a_{z23} \partial_z \psi_1 + b_{21} \phi_1 + b_{22} \psi_{v 1} + b_{23}\psi_1 $$

$$\partial_t \psi_1 = v_{\psi_1} \partial_\rho \psi_1 + a_{z31} \partial_z \phi_1 + a_{z32} \partial_z \psi_{v 1} + a_{z33} \partial_z \psi_1 + b_{31} \phi_1 + b_{32} \psi_{v 1} + b_{33}\psi_1 $$

- for the Characteristic Initial Boundary Value Problem (CIBVP):
$$\partial_x \phi_2 = a_{z11} \partial_z \phi_2 + a_{z12} \partial_z \psi_{v 2} + a_{z13} \partial_z \psi_2 + b_{11} \phi_2 + b_{12} \psi_{v 2} + b_{13}\psi_2 $$ 

$$\partial_x \psi_{v2} = a_{z21} \partial_z \phi_2 + a_{z22} \partial_z \psi_{v 2} + a_{z23} \partial_z \psi_2 + b_{21} \phi_2 + b_{22} \psi_{v 2} + b_{23}\psi_2 $$

$$\partial_u \psi_2 = v_{\psi_2} \partial_\rho \psi_2 + a_{z31} \partial_z \phi_2 + a_{z32} \partial_z \psi_{v 2} + a_{z33} \partial_z \psi_2 + b_{31} \phi_2 + b_{32} \psi_{v 2} + b_{33}\psi_2 $$

The fields $\phi_1, \psi_{v1}, \phi_2, \psi_{v2}$ are left-moving, whereas $\psi_1, \psi_2$ right-moving. The paremeters $v_{\phi_1}, v_{\psi_{v1}}, v_{\psi_1}, v_{\psi_2}$ control the speeds of the fields and should receive only positive values to maintain the direction of propagation of the fields, as well as the correct prescription of boundary data. The speeds of $\phi_2, \psi_{v2}$ are fixed to $1$. The parameters $a_{zij}, b_{ij}$, with $i,j=1,2,3$, control the angular principal part and sources of the systems, respectively, and can be tuned independently for the IBVP and CIBVP.

In the paper, the speeds of propagation are fixed to $\phi_1 = \psi_{v1} = \psi_1 = 1$ and $\psi_2 = 0.5$. Furthermore, the models explored are the symmetric hyperbolic $a_{z12}=a_{z21}=a_{z33}=1$, and the weakly hyperbolic $a_{z21}=a_{z33}=1$, with the rest of the $a_{zij}$ vanishing in each case. Regarding source terms, the following cases are explored: 
+ homogeneous i.e. $b_{ij} = 0$ for all $i,j$
+ inhomogeneous with only non-vanishing $b_{13}=1$
+ inhomogeneous with only non-vanishing $b_{32}=1$

The code is written in the [Julia programming languange](https://julialang.org/) as a module,
and tested in Julia version 1.8.5.

## Installation

Change to your local directory where the repository "model_CCE_CCM_public" is
saved. The module can be installed using Julia's REPL mode:
```
julia> ]
pkg> add .
```
