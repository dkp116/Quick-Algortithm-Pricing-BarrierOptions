# New refactored Build — Documentation

This directory contains the refactored implementation of the barrier option pricing system.  
The `main.cpp` file demonstrates how to price a **down-and-out European call option** using a **Merton Jump Diffusion model** and a **Strategy-based pricing engine**.

---

##  High-Level Architecture

The system separates responsibilities into **dynamics**, **products**, and **pricing engines**, allowing for clean, modular, and extensible code.

### Dynamics
- Asset price evolution is encapsulated in classes implementing the `IDynamics` interface.
- Current implementation uses **Merton Jump Diffusion**.
- Supports swapping to other models, e.g., Black–Scholes, without changing the rest of the code.

### Stock
- Represents the underlying asset.
- Contains the initial price and a reference to the chosen dynamics model.
- Separates the **asset identity** from its **price evolution**.

### Options
- Base `Option` class with specialized derivatives.
- Example here: **Down-and-Out European Call**.
- Defines strike price, barrier level, maturity, and exercise type.

### Pricing Engine
- Implemented using the **Strategy pattern**.
- Pricing algorithms are separate from products and dynamics.
- Example engines: `StandardMonteCarlo`, `TaylorApproximation`.
- Enables plugging in new pricing methods (quasi-Monte Carlo, Fourier, PDE solvers, etc.) without changing the product or dynamics code.

---

##  Program Flow

1. **Select a dynamics model** (e.g., Merton Jump Diffusion).  
2. **Construct a Stock** with the chosen dynamics.  
3. **Define the option** (type, strike, barrier, maturity).  
4. **Select a pricing engine** implementing `IPricing`.  
5. **Compute the option price**.  
6. **Output the numerical result**.

---

## Design Patterns Used

### Strategy Pattern
- `IDynamics` encapsulates price path simulation.
- `IPricing` encapsulates the pricing method.
- 'Option' ecapsulates the different options that can be priced.
- Enables interchangeable models and pricing engines.

### Composition
- `Stock` contains a dynamics object.
- Pricing engine contains the stock and option objects.
- Reduces inheritance complexity and improves maintainability.

### Factory Pattern (planned)
- Architecture is ready for factories:
  - `DynamicsFactory`
  - `OptionFactory`
  - `PricingEngineFactory`
- Would allow object creation without exposing `make_shared` calls in user code.

---

##  Benefits of Refactor

- **Modularity**: Easily extend or replace any component (dynamics, options, pricing engines).  
- **Readability**: Clear separation of responsibilities.  
- **Testability**: Each component can be independently unit-tested.  
- **Extensibility**: Future integration of new pricing algorithms or models requires minimal changes.  

---

## Notes

- The system currently demonstrates pricing a down-and-out European call using **Taylor Approximation**.  
- Additional engines (e.g., Monte Carlo, Uniform Sampling) are available and can be swapped with minimal code changes.  
- Designed with future Factory pattern integration in mind for fully automated object creation.
