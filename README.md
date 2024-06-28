# libNonUniformROUND: An OpenFOAM Library of the ROUND scheme on non-uniform meshes

The library provides the ROUND schemes on Non-uniform meshes (ROUND-N). The library can offer an improved structure-preserving property that gives essentially oscillation-free solutions and preserves the structures of the passive transported scalar.

## Compilation

The libROUNDSchemes library does not require any third-party dependency.
After sourcing the OpenFOAM version, simply execute:

```
./Allwmake
```
## Usage

The library can work with any OpenFOAM solvers containing convection terms.

* The library should be linked to the solver. Add the following to your system/controlDict file:

```
libs
(
    "libROUNDSchemes.so" 
);
```
## Validation and Verification
### 1. Advection of complex profiles on non-uniform meshes.
As shown in the figure below, the ROUND scheme produces low-dissipative results and preserves the structure of the passive transported scalar.
![ROUNDN](https://github.com/advanCFD/libNonUniformROUND/assets/118991833/c512d4e3-d565-4d3f-a18f-aec53f7292b5)

## Citation
Deng, X., Jiang, Z.H., and Yan, C., Efficient ROUND schemes on non-uniform grids applied to discontinuous Galerkin schemes with Godunov-type finite volume sub-cell limiting.
