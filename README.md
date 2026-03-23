# Quick-Algorithm-Pricing-BarrierOptions

This repository contains high-performance C++ implementations of advanced Monte Carlo algorithms for pricing **down-and-out barrier options** under the **Merton Jump Diffusion (MJD)** model.  
The project was developed as part of my MSc dissertation in Mathematical Finance.
[Please read the dissertation PDF](docs/Efficient_Barrier_Pricing_Daneel_Patel_final.pdf)

---

## Project Overview

Accurate pricing of path-dependent derivatives remains a central challenge in quantitative finance — especially when markets exhibit sudden jumps that classical models cannot capture. Barrier options are particularly sensitive to the trajectory of the underlying asset, making their valuation both mathematically and computationally demanding.

This project explores and extends two specialised Monte Carlo algorithms designed for pricing barrier options in the MJD framework. These algorithms were originally proposed with only sketch-level derivations; the dissertation develops them rigorously from first principles and implements them in C++ for high-performance simulation.

---

## Summary of the Research

The work focuses on overcoming the difficulties traditional Monte Carlo methods face when pricing barrier options under jump–diffusion dynamics.

### Key Contributions

- **Rigorous derivation** of two efficient barrier-crossing algorithms:
  - **Uniform Sampling Method** – samples a random point in an extended interval and uses Brownian bridge logic to detect crossings.
  - **Taylor Expansion Method** – approximates the barrier-crossing probability using a truncated Taylor series and numerical integration.
- **Use of Brownian bridge interpolation** to capture barrier hits between time steps.
- **Variance reduction via a tailored control variate**, achieving up to **50% reduction in estimator variance**.
- **Full C++ implementation** for speed and reproducibility.

These techniques provide faster and more reliable pricing, especially when barrier hits are rare — a setting where standard Monte Carlo often performs poorly.

---



- **[code/](code/README.md)**  
  Contains an ongoing refactor using **Strategy** and **Factory** patterns for improved modularity and extensibility. This will include the build instructions as well.

---

## Technologies

- **C++**
- Monte Carlo simulation
- Stochastic calculus & jump-diffusion modelling

---

## Future Work

- Extend to *up-and-out* and *double* barrier options  
- Explore quasi-Monte Carlo and additional variance reduction  
- GPU-accelerated versions (CUDA/OpenCL)

