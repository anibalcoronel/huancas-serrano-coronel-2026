# huancas-serrano-coronel-2026
 Simulation scripts and minimum working examples for Section 4 (Huancas, Serrano &amp; Coronel, 2026).
# Non-Autonomous Eco-Climatic Plant–Pollinator Dynamics (MATLAB Codes)

This repository contains the MATLAB simulation scripts and minimum working examples (MWE) for Section 4 of the manuscript:

> **Huancas, F., Serrano, F., & Coronel, A. (2026).** *Existence, Multiplicity, and Stability of Periodic Solutions in a Non-Autonomous Eco-Climatic Plant–Pollinator System*.

---

## Requirements

* **MATLAB** (tested on R2023b / R2024a; base installation only, no external toolboxes required).
* Standard built-in ODE solvers: `ode45`, `ode15s`, and a self-contained fixed-step `RK4` shooting solver.

---

## Repository Structure

All simulation scripts automatically create the output directory `FIGURASERA/` (if it does not exist) and export all figures in high-resolution `.png` (300 DPI) and vector `.eps` formats.

| Script | Description | Numerical Method & Tolerances | Generated Output (`FIGURASERA/`) |
| :--- | :--- | :--- | :--- |
| **`main_examples_for_ERA.m`** | **Master driver script** that executes all numerical experiments sequentially. | Calls all scripts below | All `.png` and `.eps` figures |
| **`exampleunoERA.m`** | **Scenario A:** Convergence to a unique globally attracting $\omega$-periodic orbit ($\lambda = 1.0$, $t \in [0, 25]$ yr). | `ode45` (`RelTol = 1e-6`, `AbsTol = 1e-8`) | `convergenia_periodic_orbit.{png,eps}` |
| **`exampledosERA.m`** | **Scenario B:** Topological multiplicity and coexistence of low- and high-amplitude regimes under superlinear growth ($\lambda = 0.05$, $t \in [0, 1.5]$ yr). | `ode15s` (`RelTol = 1e-6`, `AbsTol = 1e-8`) | `scenario_B_multiplicity.{png,eps}` |
| **`exampletresERA.m`** | **Scenario C (Part 1):** Homotopic bifurcation comparing sustainable ($\lambda = 0.5 < \lambda_c$) vs. over-forced ($\lambda = 2.5 > \lambda_c$) regimes ($t \in [0, 25]$ yr). | `ode45` (`RelTol = 1e-6`, `AbsTol = 1e-8`) | `scenario_C_bifurcation.{png,eps}` |
| **`collapse.m`** | **Scenario C (Part 2):** Physical tipping point in the original model ($\lambda = 1.0$), comparing baseline ($Q_0 = 2.0, h_0 = 0.4$) vs. thermal collapse ($Q_0 = 8.0, h_0 = 25.0$). | `ode45` (`RelTol = 1e-6`, `AbsTol = 1e-8`) | `scenario_C_physical_collapse.{png,eps}` |
| **`example_Q0.m`** | **Emission Cap Validation:** Compares the analytical threshold formula $Q_0^*(\gamma_3)$ against calibrated full non-autonomous ODE trajectories ($t \in [0, 50]$ yr). | `ode45` (`RelTol = 1e-6`, `AbsTol = 1e-8`) | `scenario_E_validation.{png,eps}` |
| **`example_sensitividad_partA.m`** | **Sensitivity Curves (Part A):** Multi-parameter sweep of periodic orbit statistics ($\bar{y}, y_{\min}, \bar{T}, T_{\max}$) across $Q_0$, $\alpha_2$, and $\gamma_3$. | Custom `RK4` ($h = 2\times 10^{-3}$) + Newton shooting (`tol = 1e-10`) | `sensitivity_curves_scan.{png,eps}` |
| **`example_sensitividad_partB.m`** | **Sensitivity Curves (Part B):** Bisection scan (16 iterations) computing the critical bifurcation threshold $\lambda_c$ across $Q_0$, $\alpha_2$, and $\gamma_3$. | `ode45` (`RelTol = 1e-5`, `AbsTol = 1e-7`) + Bisection | `scenario_D_sensitivity.{png,eps}` |

---

## How to Run

1. Clone or download this repository and open MATLAB in the repository root folder.
2. Run the master script from the MATLAB Command Window:
   ```matlab
   main_examples_for_ERA
   ```
3. All figures will be generated and saved automatically inside the `FIGURASERA/` folder. Individual scenarios can also be reproduced independently by running any of the `.m` scripts listed above.
