# vchem-insurance-simulation

# Methods

## Data Sources and Study Population

We analyzed health insurance transitions using longitudinal data from the Medical Expenditure Panel Survey (MEPS) 
for two periods: 2011-2013 (pre-ACA) and 2021-2023 (post-ACA). MEPS provides monthly insurance status through 
overlapping two-year panels of approximately 15,000 households per year, with oversampling of low-income and 
minority populations. We categorized insurance coverage into four mutually exclusive 
states: employer-sponsored insurance, other private insurance, public insurance (Medicaid/Medicare), 
and uninsured. 

Cross-sectional insurance coverage distributions by age were obtained from the American Community Survey (ACS) 
for 2012 and 2022 to serve as calibration targets, with approximately 3 million observations per year.

## Statistical Analysis

We developed a continuous-time Markov chain model with state space 

$$
S = \{\text{Employer}, \text{OthPrivate}, \text{Public}, \text{Uninsured}, \text{Death}\}
$$

where Death serves as an absorbing state. For each age $a$ and period $p$ (pre/post-ACA), we estimated a
$5 \times 5$ transition intensity matrix $\mathbf{Q}(a,p)$ where element $q_{ij}$ represents the 
instantaneous rate of transition from state $i$ to state $j$. The diagonal elements $q_{ii}$ were
defined as $-\sum_{j \neq i} q_{ij}$ to ensure proper row sums. Background mortality rates 
$\mu(a,p)$ were incorporated from the Human Mortality Database using period life tables.
Initial maximum likelihood estimates of transition intensities were computed as:

$$
\hat{q}_{ij}(a,p) = \frac{n_{ij}(a,p)}{n_i(a,p) \times \Delta t}
$$

where $n_{ij}(a,p)$ represents the observed number of transitions from state $i$ to $j$ at age $a$ in period 
$p$, $n_i(a,p)$ is the total person-time observed in state $i$, and $\Delta t$ is the observation interval (1 month).

To address sampling variability while preserving meaningful age-related changes, we employed generalized additive 
models with transition-specific smoothing parameters:

$$
\log(q_{ij}(a)) = s(a; k = k_{ij}, b = b_{ij})
$$

where $s(\cdot)$ represents a cubic regression spline, $k_{ij}$ is the number of basis functions, 
and $b_{ij}$ specifies the basis type. For transitions involving policy-relevant age thresholds (e.g., 
Medicaid eligibility at age 19), we used adaptive smoothing ($b = \text{"ad"}$) with increased flexibility 
($k = 30$). Employer-sponsored insurance transitions were modeled using cubic regression splines 
($b = \text{"cr"}$, $k = 25$) to capture graduation and early career patterns, while other transitions 
used cubic smoothing splines ($b = \text{"cs"}$, $k = 20$).
The smoothed transition rates were calibrated to match ACS cross-sectional distributions using 
Metropolis-Hastings MCMC. For each iteration $t$:

1. Select a random transition type $(i,j)$
Propose new rates: $q^{ij}(a) = q_{ij}(a) \times \exp(\varepsilon)$, where $\varepsilon \sim N(0, 0.01)$
3. Calculate state occupancy probabilities $\boldsymbol{\pi}(a)$ by solving the forward Kolmogorov equations:

$$
\frac{d\boldsymbol{\pi}(a)}{da} = \boldsymbol{\pi}(a)\mathbf{Q}(a)
$$

4. Accept proposal with probability $\min(1, \exp(-[SSE^ - SSE]/\tau))$
where $SSE$ represents the sum of squared errors between predicted and ACS-observed state occupancies, and $\tau = 0.01$ 
is the temperature parameter. The chain was run for 10,000 iterations with convergence monitored via acceptance rates 
(target: 0.5-0.8) and sequential improvements in fit metrics.

All analyses were performed using R software (version 4.2.0, R Foundation for Statistical Computing) with the 
mgcv package for GAM estimation and expm for matrix exponentials. The study protocol was deemed exempt by
[Institution] IRB as it used only publicly available data.