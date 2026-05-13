# CFD Assignment 2

This repository contains the core MATLAB implementations for CFD Assignment 2a and 2b.

## Files

- `conv_ee_fillin.m`: 1D unsteady pure advection equation with central differences in space and explicit Euler time integration.
- `conv_ie_fillin.m`: 1D unsteady pure advection equation with central differences in space and implicit Euler time integration.
- `conv_diff_cn_fillin.m`: 1D unsteady advection-diffusion equation with central differences in space and Crank-Nicolson time integration.

## Numerical Problems

Assignment 2a solves

```text
d(phi)/dt = -U0 d(phi)/dx
```

with periodic boundary conditions on `x in [0, 2*pi]` and initial condition `phi(x,0) = sin(x)`.

Assignment 2b solves

```text
d(phi)/dt = -U0 d(phi)/dx + Gamma d2(phi)/dx2
```

with the same periodic domain and initial condition.

Each script plots the numerical solution against the analytical solution so that stability, dissipative error, and dispersive error can be inspected visually.
