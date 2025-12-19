#pragma once
#include <memory>
#include "Stock/Stock.h"
#include "Options/Option.h"
#include "Dynamics/IDynamics.h"
#include "PricingAlgorithm/UniformSampleEstimate.h"
#include "Dynamics/MertonJumpDynamics.h"

// A helper function to price a Down-and-Out call option with Monte Carlo
double price_down_and_out_call_with_taylor_series();
double price_down_and_out_call_with_uniform_distribution();
