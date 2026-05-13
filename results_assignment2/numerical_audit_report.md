# Numerical Audit Report for CFD Assignment 2a and 2b

## 1. Purpose of the numerical audit

The goal of this workflow is to verify the analytical solution implementation, the consistency of periodic boundary conditions, and the expected stability, dissipation, and dispersion behaviour of the three numerical schemes.  The tests also document the sensitivity to advection velocity `U0`, diffusion coefficient `Gamma`, grid spacing `dx`, and timestep `dt`.

## 2. Verification metrics

- **Relative L2 error**: `norm(phi - phi_a, 2) / norm(phi_a, 2)`. This measures the global profile error relative to the analytical solution.
- **Linf error**: `max(abs(phi - phi_a))`. This catches the largest pointwise deviation.
- **Periodic boundary mismatch**: `abs(phi(1) - phi(end))`. This checks whether the duplicated periodic endpoint remains consistent.
- **Numerical amplitude**: `(max(phi) - min(phi)) / 2`. This is used to diagnose artificial growth or damping.
- **Analytical amplitude**: `1` for pure advection and `exp(-Gamma*t)` for advection-diffusion.
- **Courant number**: `C = U0*dt/dx`. This controls the advective timestep size.
- **Fourier number**: `F = Gamma*dt/dx^2`. This measures the diffusive timestep size.

## 3. Implementation sanity checks

The periodic mismatch should remain close to machine precision because `x=0` and `x=2*pi` represent the same physical point.  For `U0 > 0`, the sine wave should propagate in the positive x direction.  For `U0 = 0`, no translation should occur.  For `Gamma = 0`, the analytical solution has no physical amplitude decay, while for `Gamma > 0` the analytical amplitude is `exp(-Gamma*t)`.

- Representative EE-2 final periodic mismatch range: 0.000e+00 to 4.441e-16.
- Representative EE U0-study stopping times: 1 to 1.
- Representative IE-2 final periodic mismatch range: 0.000e+00 to 6.661e-16.
- Representative IE U0-study final amplitude range: 0.97967 to 0.998355.
- CN-1 final amplitude: numerical 0.904703, analytical 0.904837.
- CN-3 pure diffusion final amplitude: numerical 0.904849, analytical 0.904837.
- CN U0-study final analytical amplitude is 0.904837 for all U0 values.

## 4. Results for Explicit Euler + Central Difference

The explicit central advection scheme exhibits unstable amplitude growth. A smaller Courant number delays visible blow-up, but it does not remove the fundamental instability of the FTCS-like central scheme for pure advection. The larger Courant-number cases show faster growth in relative L2 error and `max(abs(phi))`.

Relevant outputs:
- `figures/EE_short_time_profile.png`
- `figures/EE_Courant_L2_history.png`
- `figures/EE_Courant_maxabs_history.png`
- `figures/EE_Courant_1_profile.png`
- `figures/EE_U0_L2_history.png`
- `figures/EE_U0_maxabs_history.png`
- `figures/EE_grid_final_L2.png`
- `tables/table_EE_Courant_study.csv`
- `tables/table_EE_U0_study.csv`
- `tables/table_EE_grid_study.csv`

## 5. Results for Implicit Euler + Central Difference

The implicit Euler scheme remains bounded for all tested Courant numbers. However, larger Courant number increases numerical dissipation, which is visible as stronger amplitude decay and larger final L2 error.  The U0 sweep shows the same trend because larger U0 increases the Courant number when dt and dx are fixed.  This confirms that a stable scheme can still be inaccurate when the timestep is too large.

Relevant outputs:
- `figures/IE_baseline_profile.png`
- `figures/IE_Courant_L2_history.png`
- `figures/IE_Courant_amplitude_history.png`
- `figures/IE_Courant_profile_C010_C200.png`
- `figures/IE_U0_L2_history.png`
- `figures/IE_U0_amplitude_history.png`
- `figures/IE_grid_final_L2.png`
- `tables/table_IE_Courant_study.csv`
- `tables/table_IE_U0_study.csv`
- `tables/table_IE_grid_study.csv`

## 6. Results for Crank-Nicolson Advection-Diffusion

### 6.1 Coupled advection-diffusion baseline

The baseline case translates and decays as expected.  The final numerical amplitude 0.904703 is close to the analytical amplitude 0.904837, and the final relative L2 error is 2.586e-04.

### 6.2 Pure advection limit Gamma = 0

With `Gamma = 0`, no physical amplitude decay is present.  The solution stays bounded in all tested Courant-number cases, while larger Courant number makes phase/dispersive error more visible.

### 6.3 Pure diffusion limit U0 = 0

With `U0 = 0`, the wave does not translate.  The sinusoidal shape is preserved while the amplitude decays approximately as `exp(-Gamma*t)`.  The pure diffusion final relative L2 error is 1.285e-05.

### 6.3b U0 sensitivity at fixed Gamma

With fixed `Gamma = 0.1`, all U0 cases have the same analytical amplitude decay `exp(-0.1)`.  Changing U0 changes the translation distance and Courant number, so phase/dispersive error can vary while the final amplitudes remain close to the same analytical value.

### 6.4 Gamma sensitivity

Increasing `Gamma` increases physical damping.  The final numerical amplitudes track the analytical amplitudes in `figures/CN_Gamma_final_amplitude.png`.

### 6.5 dt and grid sensitivity

The timestep study shows that larger `dt` increases the final error, while the grid study shows decreasing final L2 error under mesh refinement.  These trends are consistent with Crank-Nicolson time integration and central spatial differences.

Relevant outputs:
- `figures/CN_baseline_profile.png`
- `figures/CN_baseline_L2_history.png`
- `figures/CN_pure_advection_C100_profile.png`
- `figures/CN_pure_diffusion_amplitude_history.png`
- `figures/CN_pure_diffusion_profile.png`
- `figures/CN_U0_profiles.png`
- `figures/CN_U0_L2_history.png`
- `figures/CN_Gamma_final_amplitude.png`
- `figures/CN_dt_final_L2.png`
- `figures/CN_grid_final_L2.png`
- `tables/table_CN_baseline_and_limits.csv`
- `tables/table_CN_U0_study.csv`
- `tables/table_CN_Gamma_study.csv`
- `tables/table_CN_dt_study.csv`
- `tables/table_CN_grid_study.csv`

## 7. Suggested report-ready discussion text

The explicit Euler scheme with central spatial differencing is not suitable for stable long-time integration of the pure advection equation.  Although small Courant numbers and short integration times can make the numerical profile appear reasonable, the error and amplitude histories reveal the unstable amplification expected from the FTCS-like central scheme.  Increasing the Courant number makes the instability visible more quickly.

The implicit Euler central scheme removes the blow-up observed in the explicit scheme, but it does so at the price of numerical damping.  Larger timesteps remain bounded, yet the wave amplitude decays too strongly and the phase error increases. Therefore stability alone is not sufficient for accuracy.

For the advection-diffusion equation, the Crank-Nicolson central scheme gives a more balanced result.  The coupled baseline agrees well with the analytical solution, including both translation and exponential amplitude decay.  The limiting cases confirm that the implementation reduces correctly to pure advection when `Gamma = 0` and pure diffusion when `U0 = 0`.

## 8. Final checklist

[ ] All periodic boundary mismatches are small
[ ] Numerical propagation direction is correct
[ ] EE instability observed
[ ] IE numerical damping observed
[ ] CN pure advection case behaves reasonably
[ ] CN pure diffusion case matches exp(-Gamma*t)
[ ] Parameter tables exported
[ ] Figures exported
