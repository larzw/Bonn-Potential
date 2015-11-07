#One-Boson-Exchange potential and Nucleon-Nucleon scattering
This FORTRAN 77 code encompasses the relativistic momentum space One-Boson-Exchange potential (OBEP) and Nucleon-Nucleon scattering. It was written by R. Machleidt and appeared in *[R. Machleidt, Computational Nuclear Physics 2 Nuclear Reactions, (Springer-Verlag, 1993), K. Langanke, J.A. Maruhn, S.E. Koonin eds., Chap. 1, pp. 1-29]*. Additional resources can be found in

- *R. Machleidt, K. Holinde, and C. Elster, Phys. Rept. 149, 1 (1987)*
- *R. Machleidt, Adv. Nucl. Phys. 19, 189 (1989)*

##Quick Start
To compile the code use `f77` or `g77`
```
f77 bonn.f phases.f cph.f
```
Inputs for the phase shift calculation are contained in *dphbonnbpv.d* and the outputs in *phbonnbpv.d*

##Bonn Potential
The subroutine `bonn.f` computes the momentum-space OBEP in terms of partial-wave helicity matrix elements. Although the code uses several subroutines, the user needs to call only `bonn` without worrying about the other computer programs contained in the package.

##Phase Shifts
Phase shifts are calculated from the *R-matrix*. The *R-matrix* is calculated using matrix inversion *[Michael I. Haftel and Frank Tabakin, Nucl. Phys. A 158, 1 (1970)]* and the Bonn Potential.
