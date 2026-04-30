 # Barrier Option Pricing — Code Documentation

 This folder contains a refactored C++ implementation of barrier option pricing under the Merton Jump Diffusion model.
 It uses a modular design to separate stochastic dynamics, option definitions, and pricing algorithms.

 ## Build Instructions

 ```bash
 cd code
 mkdir -p build
 cd build
 cmake ..
 cmake --build .
 ```

 The build produces:

 - `main` — the example pricing executable
 - `tests` — the Catch2 test binary (enabled by default)

 ## Running the example

 ```bash
 ./main
 ```

 ## Running tests

 ```bash
 ctest --output-on-failure
 ```

 ## What the example does

 The executable in `src/main.cpp` demonstrates pricing a down-and-out European call option using:

 - `MertonJumpDynamics` for asset dynamics
 - `Stock` for the underlying asset
 - `DownAndOut` for the barrier option product
 - `UniformSample` as the pricing engine

 The example computes:

 - `PriceWithVarianceReduction()`
 - `Price()`
 - standard error via `GetStandardError()`

 ## Available components

 ### Dynamics

 Implemented dynamics models live under `src/Dynamics/` and are exposed through `include/Dynamics/`.
 
 Current implementations:

 - `MertonJumpDynamics` — jump-diffusion with Poisson jumps
 - `BlackScholesDynamics` — continuous diffusion without jumps

 ### Stock

 `Stock` binds an initial price to a dynamics model so pricing engines can simulate the underlying asset.

 ### Options

 The option hierarchy defines payoffs and barrier parameters.
 Current implementation:

 - `DownAndOut` — down-and-out barrier option

 ### Pricing engines

 Implemented pricing methods include:

 - `StandardMonteCarlo`
 - `UniformSample`
 - `TaylorApproximation`

 Each pricing engine implements `IPricing` and can price any supported product under any supported dynamics model.

 ## How to use the code

 1. Choose a dynamics model, e.g. `std::make_shared<MertonJumpDynamics>(...)`
 2. Create a `Stock` using the dynamics model
 3. Create an `Option` product such as `DownAndOut`
 4. Construct a pricing engine with the stock, option, and simulation settings
 5. Call `Price()` or `PriceWithVarianceReduction()`

 ## Why this design

 The architecture is built for flexibility:

 - pricing methods are independent of option product definitions
 - option products are independent of model dynamics
 - new models or pricing algorithms can be added with minimal changes

 ## Notes

 - The code is written in modern C++ (C++23).
 - The current example focuses on a down-and-out European call under a Merton jump diffusion process.
