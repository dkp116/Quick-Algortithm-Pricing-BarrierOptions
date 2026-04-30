# Quick-Algorithm-Pricing-BarrierOptions

A C++ implementation of fast Monte Carlo pricing methods for **down-and-out barrier options** under the **Merton Jump Diffusion** model.
This repository contains a modular pricing system, example executable, and a test suite used to validate the pricing algorithms.

## Repository contents

- `code/` — main C++ implementation, headers, source files, and build configuration
- `code/src/main.cpp` — example program that constructs a barrier option and computes prices
- `code/tests/` — Catch2-based tests for pricing methods
- `docs/` — dissertation document with mathematical derivation and research results

> See `docs/Efficient_Barrier_Pricing_Daneel_Patel_final.pdf` for the full dissertation and theoretical background.

## What this code does

The code implements pricing for barrier options using several numerical methods.
It is designed for modularity and extensibility, separating:

- asset dynamics (`Dynamics`)
- underlying security (`Stock`)
- option payoffs and barriers (`Option`, `Barrier`)
- pricing algorithms (`IPricing` implementations)

Supported pricing engines include:

- `StandardMonteCarlo`
- `UniformSample`
- `TaylorApproximation`

## Build and run

In the `code/` directory:

```bash
cd code
mkdir -p build
cd build
cmake ..
cmake --build .
```

Then run the example:

```bash
./main
```

Run the test suite:

```bash
ctest --output-on-failure
```

## Example usage

The example program demonstrates a workflow with:

1. `MertonJumpDynamics(r, sigma, lambda, mu_jump, sigma_jump)`
2. `Stock(S0, dynamics)`
3. `DownAndOut(ExerciseType::European, OptionType::Call, strike, barrier, maturity)`
4. a pricing engine such as `UniformSample`, `TaylorApproximation`, or `StandardMonteCarlo`
5. calling `Price()` or `PriceWithVarianceReduction()`

The `main` executable prints a price estimate and the standard error.

## Tests and validation

The project includes integration tests that verify down-and-out call pricing for different methods.
The tests use a reference price and validate the algorithms with Catch2.

## Project structure

- `code/include/` — header files for dynamics, options, stock, and pricing algorithms
- `code/src/` — implementation of the library and the example program
- `code/tests/` — test implementations and helper functions

## Extending the project

To add new models or products:

- implement a new `IDynamics` model for a different stochastic process
- extend the `Option` hierarchy for new derivative payoffs
- reuse existing pricing engines or add a new `IPricing` implementation

## Notes

- The project is built with C++23.
- The example program currently prices a down-and-out European call under a Merton jump diffusion model.
- The modular design allows easy swapping of models, products, and pricing engines.
