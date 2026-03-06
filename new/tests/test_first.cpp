#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "test_utils/TestDownAndOutPricing.h"

TEST_CASE("Down-and-Out call using Taylor Series approximation price check", "[montecarlo]") {
    double price = price_down_and_out_call_with_taylor_series();
    double expected = 9.0;

    CAPTURE(price, expected);  // automatically prints both if REQUIRE fails
    INFO("Testing Taylor Series Monte Carlo pricing for Down-and-Out Call");

    REQUIRE(price == Catch::Approx(expected).epsilon(0.1));
}

TEST_CASE("Down-and-Out call using Uniform distribution price check", "[montecarlo]") {
    double price = price_down_and_out_call_with_uniform_distribution();
    double expected = 9.0;

    CAPTURE(price, expected);
    INFO("Testing Uniform sampling Monte Carlo pricing for Down-and-Out Call");

    REQUIRE(price == Catch::Approx(expected).epsilon(0.1));
}

TEST_CASE("Down-and-Out call using Standard Monte Carlo price check", "[montecarlo]") {
    double price = price_down_and_out_call_with_standard_monte_carlo();
    double expected = 9.0;

    CAPTURE(price, expected);
    INFO("Testing Standard Monte Carlo pricing for Down-and-Out Call");

    REQUIRE(price == Catch::Approx(expected).epsilon(0.5));
}


TEST_CASE("Down-and-Out call using Standard Monte Carlo price check with Varience", "[montecarlo]") {
    double price = price_down_and_out_call_with_standard_monte_carlo_and_varience();
    double expected = 9.0;

    CAPTURE(price, expected);
    INFO("Testing Standard Monte Carlo pricing for Down-and-Out Call");

    REQUIRE(price == Catch::Approx(expected).epsilon(0.5));
}



