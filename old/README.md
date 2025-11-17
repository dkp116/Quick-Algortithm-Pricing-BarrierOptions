# Old Build —  Documentation 

This directory contains the original C++ implementation from my MSc dissertation.  
The `main.cpp` file acts as the driver for running and comparing multiple Monte Carlo–based pricing algorithms for **down-and-out barrier options** under the **Merton Jump Diffusion (MJD)** model.

The program benchmarks three methods:

1. **Uniform Sampling Algorithm**
2. **Taylor Expansion Algorithm**
3. **Standard Monte Carlo**
4. **Control Variate–enhanced estimator** (using the MJD European call price as the control)

It outputs the estimated option price and standard errors for each method across a range of simulation sizes.

---

## What `main.cpp` contains

### 1. **Implements Helper Analytics**
The file includes:
- A standard normal CDF (`norm_cdf`)
- Black–Scholes European call price
- Closed-form approximations for down-and-out barrier options (with and without rebate)
- A Poisson-summation pricing function for European calls under **Merton Jump Diffusion**

These serve both as validation tools and as the **control variate** foundation.

---

##  Merton Jump Diffusion Pricing

The function `PriceMJD()` computes the European call price under MJD by summing over the Poisson distribution of jump counts:

- For each possible number of jumps `n`, it adjusts:
  - Effective drift
  - Effective volatility
- Computes the Poisson probability of having exactly `n` jumps
- Prices using a Black–Scholes call with the adjusted parameters
- Aggregates the weighted prices

This value is used later as the **control variate target**.

---

## Barrier Option Simulation Methods

For each simulation run, the program calls:

### **Uniform Sampling Method**
`Derivative.UniformVarRedCall(stock)`  
Returns:
- Estimated discounted payoff  
- The control variate value  

This method uses Brownian bridge logic with a uniform-sampled point inside an extended time interval to detect barrier crossings more efficiently.

---

### **Taylor Expansion Method**
`Derivative.PriceByMJD_Taylor(stock)`  
Computes the crossing probability using a truncated Taylor expansion of the integral expression.  
Useful when barriers are rarely hit.

---

### **Standard Monte Carlo**
`Derivative.StandardMonteCarlo(stock)`  
A baseline Euler-discretised Monte Carlo simulation.

---

## Control Variate Estimator

The program computes:
- Covariance between the uniform-sampling estimator and control variate
- Optimal β coefficient  
- Final variance-reduced estimator:  
  \[
  \hat{X}_{CV} = \hat{X} - \beta \left( \hat{Y} - E[Y] \right)
  \]

Where `E[Y]` is the analytic MJD call price from `PriceMJD()`.

This typically cuts variance by **≈50%**, especially when the barrier is far from the spot price.

---

## Simulation Loop

The program runs the pricing algorithms for increasing simulation counts:

