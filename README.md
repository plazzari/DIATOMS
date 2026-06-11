# Diatom Size-Structured Population Model

A Python implementation of a size-structured diatom population model describing the evolution of diatom populations along the **Cell Size Reduction and Restitution Cycle (CSRR)**.

The model combines biological allometries, mortality processes, sexual reproduction, and density-dependent regulation within a transport-reaction framework. It provides both a numerical solution of the governing partial differential equation and a characteristic-based analytical approximation for comparison.

---

## Features

- Size-structured population dynamics
- Allometric growth rates based on cell volume
- Size-dependent mortality
- Sexualization and size restitution processes
- Logistic population regulation
- Numerical PDE integration using an upwind scheme
- Characteristic-based analytical solution
- Visualization of population dynamics in size and time
- Comparison between numerical, analytical, and equilibrium solutions

---

## Mathematical Model

The population density $n(t,s)$ evolves according to

$\frac{\partial n}{\partial t}+\frac{\partial}{\partial s}\left(\rho(s)n\right)=\left[\mu(s)-\delta(s)-\beta(s)\right]n-\frac{n^2}{K(s)}$

where:

| Symbol | Description |
|----------|-------------|
| $n(t,s)$ | Population density |
| $s$ | Cell size (µm) |
| $t$ | Time (days) |
| $\mu(s)$ | Growth rate |
| $\delta(s)$ | Mortality rate |
| $\beta(s)$ | Sexualization/restitution rate |
| $K(s)$ | Carrying capacity |
| $\rho(s)$ | Velocity in size space |

The advection term represents the progressive reduction in average cell size resulting from successive vegetative divisions.

---

## Biological Assumptions

### Growth

Growth rates are derived from empirical allometric relationships relating cell volume to population growth:

- Tang (1995)
- Irwin et al. (2006)

### Mortality

Mortality is modeled as a size-dependent process based on volume-dependent relationships inspired by:

- Shimoda et al. (2016)

### Sexualization and Size Restitution

Diatoms progressively decrease in size during vegetative reproduction. When cells approach a critical size threshold, the probability of sexual reproduction increases, allowing restoration of maximal cell size.

The model represents this process through a smooth size-dependent sexualization function $\beta(s)$.

### Density Dependence

Population regulation is introduced through a size-dependent carrying capacity $K(s)$, which may contain preferred size classes represented by Gaussian peaks.

---

## Numerical Method

The governing equation is solved on a uniform size grid using:

- Explicit upwind discretization for transport
- Logistic reaction terms
- Time stepping on a regular temporal grid
- Non-local boundary condition representing the production of rejuvenated cells

The code automatically evaluates the Courant–Friedrichs–Lewy (CFL) condition to verify numerical stability.

---

## Characteristic Solution

To validate the numerical solution, the model also integrates the characteristic equation

$\frac{dn}{ds}=\frac{(\mu(s)-\delta(s)-\beta(s))n-n^2/K(s)}{\rho(s)}$

using SciPy's implicit Radau solver.

This provides a direct comparison between:

- Numerical PDE solution
- Characteristic solution
- Local equilibrium solution

---

## Requirements

Install dependencies:

```bash
pip install numpy scipy matplotlib
```

---

## Running the Model

Execute:

```bash
python diatom_model.py
```

Choose the simulated population by setting:

```python
id = "small_diatoms"
```

or

```python
id = "large_diatoms"
```

The two configurations correspond to different size ranges and sexualization thresholds.

---

## Output

### 1. Rates and Characteristic Comparison

```text
rates_and_characteristic_<id>.png
```

Contains:

- Growth rate $\mu(s)$
- Mortality rate $\delta(s)$
- Sexualization rate $\beta(s)$
- Carrying capacity $K(s)$
- Numerical solution
- Characteristic solution
- Equilibrium solution

### 2. Size-Time Heatmap

```text
heatmap_<id>.png
```

Visualizes the evolution of the size distribution through time.

### 3. Time Series

```text
time_series_<id>.png
```

Shows temporal dynamics for selected size classes.

---

## Scientific Applications

This framework can be used to investigate:

- Diatom cell-size reduction and restitution cycles
- Size-structured phytoplankton populations
- Cohort dynamics
- Effects of size-dependent mortality
- Sexual reproduction thresholds
- Structured population theory
- Validation of analytical characteristic solutions

---

## References

- Tang, E. P. Y. (1995). *The allometry of algal growth rates.*
- Irwin, A. J., Finkel, Z. V., Schofield, O. M., & Falkowski, P. G. (2006).
- Shimoda, Y. et al. (2016).
- Sinko, J. W., & Streifer, W. (1967). *A new model for age-size structure of a population.*

---



---

## License

MIT License.

Feel free to use, modify, and distribute this software with appropriate attribution.
