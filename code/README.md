# New refactored Build — Documentation

This directory contains the refactored implementation of the barrier option pricing system.  
The `main.cpp` file demonstrates how to price a **down-and-out European call option** using a **Merton Jump Diffusion model** and a **Strategy-based pricing engine**.

---

# Build Instructions

```bash
mkdir build
cd build
cmake ..
make
./Main

```
---

# Example Pricing Function

The function below prices a down-and-out European call option using a Taylor approximation pricing engine under Merton Jump Diffusion dynamics.

```cpp
double price_down_and_out_call_with_taylor_series()
{
    // Define Merton Jump Diffusion dynamics
    auto dynamics = std::make_shared<MertonJumpDynamics>(
        0.05,  // risk-free rate
        0.25,  // diffusion volatility
        2.0,   // jump intensity (lambda)
        0.0,   // mean jump size
        0.1    // jump volatility
    );

    // Create underlying asset
    auto stock = std::make_shared<Stock>(
        100.0,     // initial price
        dynamics
    );

    // Define Down-and-Out European call option
    auto option = std::make_shared<DownAndOut>(
        ExerciseType::European,
        OptionType::Call,
        110.0,   // strike price
        85.0,    // barrier level
        1.0      // maturity in years
    );

    // Select pricing engine
    TaylorApproximation pricing(
        stock,
        option,
        100000,  // number of simulations
        StandardErrorCalculation::NotIncluded,
        Time::NotIncluded
    );

    // Compute option price
    return pricing.Price();
}
```

## Architecture Overview

The system separates responsibilities into three independent components:

- **Dynamics** → how the underlying asset evolves
- **Products** → what payoff is being priced
- **Pricing Engines** → how the price is computed

This separation promotes **modularity**, **testability**, and **extensibility**, allowing components to be modified or replaced without affecting the rest of the system.

---

### Dynamics

The dynamics layer defines the stochastic process governing asset price evolution.

All models implement the `IDynamics` interface.

Current implementation:

- **Merton Jump Diffusion**
  - continuous diffusion component (Brownian motion)
  - discontinuous jump component (Poisson arrivals)

Key benefits:

- models are interchangeable
- pricing logic remains unchanged when switching models
- new processes (e.g. Black–Scholes, Heston) can be added easily

---

### Stock

The `Stock` class represents the underlying asset.

Responsibilities:

- stores the initial asset price
- holds a reference to the chosen dynamics model

This design separates:

- asset identity
- price evolution behaviour

allowing the same product to be evaluated under different stochastic models.

---

### Options (Products)

The `Option` hierarchy represents derivative contracts.

Base class:

- `Option`

Example implementation:

- **Down-and-Out European Call**

Each product defines:

- strike price
- maturity
- barrier level (if applicable)
- exercise type (European, American, etc.)
- payoff structure

New products can be introduced by extending the base `Option` class without modifying pricing logic.

---

### Pricing Engines

Pricing engines implement numerical valuation methods via the `IPricing` interface.

Example implementations:

- `StandardMonteCarlo`
- `TaylorApproximation`

Key properties:

- pricing algorithms are independent of product definition
- pricing algorithms are independent of stochastic model
- engines can be swapped without changing product or model code

Future extensions may include:

- quasi-Monte Carlo methods
- Fourier transform methods
- PDE solvers
- adjoint differentiation methods

---

## Program Flow

Typical usage follows the workflow below:

1. Select a stochastic model  
   (e.g. Merton Jump Diffusion)

2. Construct the underlying asset
   (binds initial price and dynamics)

3. Define the derivative product
   (option type, strike, maturity, barrier)

4. Select a pricing engine
   (implementation of `IPricing`)

5. Compute the option price

6. Output the numerical result

---

## Design Patterns

### Strategy Pattern

Encapsulates interchangeable algorithms behind common interfaces.

Interfaces:

- `IDynamics` → defines asset price evolution
- `IPricing` → defines pricing methodology
- `Option` → defines payoff structure

Benefits:

- models are interchangeable
- pricing methods are interchangeable
- new functionality can be introduced with minimal modification to existing code

---

### Composition

Objects are constructed using composition rather than deep inheritance hierarchies.

Relationships:

- `Stock` contains a dynamics object
- pricing engines contain both stock and option objects

Advantages:

- improved flexibility
- simpler class hierarchies
- easier unit testing
- reduced coupling

---

### Factory Pattern (Planned)

The architecture is designed to support factory-based object creation.

Planned factories:

- `DynamicsFactory`
- `OptionFactory`
- `PricingEngineFactory`

Factories will:

- centralise object creation logic
- reduce direct use of `std::make_shared`
- simplify configuration-based workflows
- improve usability for library users

Example future usage:

```cpp
auto option = OptionFactory::Create("DownAndOutCall", parameters);
