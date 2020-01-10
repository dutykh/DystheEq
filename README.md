![Inelastic collision of four Dysthe-Lo-Mei solitary waves](/pic/4Coll.png "Inelastic collision of four Dysthe-Lo-Mei solitary waves")

# Dysthe-Lo-Mei equation

The present repository contains a Fourier-type pseudo-spectral solver for Dysthe-Lo-Mei equation as described in the references cited below. More precisely, we solve the dimensionless version of this equation, which corresponds to Equations (2.6) - (2.9) from the Lo & Mei (1985) paper. A very high order Runge-Kutta scheme is used for time integration. We employ also the integrating factor technique to slightly remove the stiffness. The work of this code is illustrated on a simple evolution of the ground state to the low order NLS equation. This solution can be computed analytically.

## Acknowledgements

The Author is grateful [Chia-Cheng (Finite) Tsai](https://finitetsai.github.io) for inciting him to share the present Matlab (TM) code.

## References

* E. Lo, C.C. Mei. [A numerical study of water-wave modulation based on a higher-order nonlinear Schrödinger equation](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/numerical-study-of-waterwave-modulation-based-on-a-higherorder-nonlinear-schrodinger-equation/2326FB7BFD6E0BBE4499E59F3FE80D39), J. Fluid Mech., **150**, 395-416, 1985

* F. Fedele, D. Dutykh. [Hamiltonian form and solitary waves of the spatial Dysthe equations](https://link.springer.com/article/10.1134%2FS0021364011240039), JETP Letters, **94**(12), 840-844, 2011

* F. Fedele, D. Dutykh. [Hamiltonian description and traveling waves of the spatial Dysthe equations](https://arxiv.org/abs/1110.3605), Research report, arXiv:1110.3605, 2012
