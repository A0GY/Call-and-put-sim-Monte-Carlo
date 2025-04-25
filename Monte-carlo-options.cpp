// A program to simulate strike prices of a stock using Monte Carlo simulation
// we output a call and put depending on factores such as the stock price, strike price, risk free rate, volatility, time to maturity, and number of simulations
// to then give a value of the option
// this is a basic example of how to use Monte Carlo simulation to price options desgin in cpp
// Author : A Y

#include <iostream>
#include <cmath>
#include <random>
#include <utility>

// MSVC portability helper
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440   //#
#endif

// Generates a single random draw from a Normal(mean, stddev) distribution
double generateGaussianNoise(double mean, double stddev) {
    static std::mt19937 generator(std::random_device{}());  // Mersenne Twister engine 
    std::normal_distribution<double> normalDist(mean, stddev);
    return normalDist(generator);
}

// Payoff function for a Call option: max(stockPrice - strikePrice, 0)
double getCallPayoff(double stockPriceAtMaturity, double strikePrice) {
    return std::max(stockPriceAtMaturity - strikePrice, 0.0);
}

// Payoff function for a Put option: max(strikePrice - stockPrice, 0)
double getPutPayoff(double stockPriceAtMaturity, double strikePrice) {
    return std::max(strikePrice - stockPriceAtMaturity, 0.0);
}

// helper function  normal distribution
inline double norm_cdf(double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); }

// closed form black scholes analytic call
double blackScholesCall(double S0, double K, double r, double sigma, double T) {
    double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) /
        (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S0 * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}   // function now correctly closed

struct McPrice {              
    double fairValue;      // discounted payoff
    double expectedPayoff; // mean at maturity
    double stdErr;         // one sigma standard error
};

// Monte Carlo simulation for pricing a European option (Call or Put)
// takes an unsigned seed param  0 means randomize with device
McPrice monteCarloOptionPricing(
    double initialStockPrice,
    double strikePrice,
    double riskFreeRate,
    double volatility,
    double timeToMaturity,
    int numberOfSimulations,
    bool isCallOption,
    unsigned seed = 0                 // fallback 
) {
    double totalPayoff = 0.0;
    double totalSqPayoff = 0.0;

    // If seed==0  use random_device else use fixed seed for reproducibility
    std::mt19937 gen(seed ? seed : std::random_device{}());
    std::normal_distribution<double> Z(0.0, 1.0);

    for (int i = 0; i < numberOfSimulations; ++i) {
        double z = Z(gen);
        // 
        double stockPriceAtMaturity = initialStockPrice * std::exp(
            (riskFreeRate - 0.5 * volatility * volatility) * timeToMaturity +
            volatility * std::sqrt(timeToMaturity) * z);
        double payoff = isCallOption
            ? getCallPayoff(stockPriceAtMaturity, strikePrice)
            : getPutPayoff(stockPriceAtMaturity, strikePrice);

        totalPayoff += payoff;
        totalSqPayoff += payoff * payoff;
    }

    double averagePayoff = totalPayoff / static_cast<double>(numberOfSimulations);
    double variance = (totalSqPayoff / numberOfSimulations) - averagePayoff * averagePayoff;
    double stdErr = std::sqrt(variance / numberOfSimulations);
    double discountedPayoff = std::exp(-riskFreeRate * timeToMaturity) * averagePayoff;

    return { discountedPayoff, averagePayoff, stdErr };
}

int main() {
    // Option parameters
    double initialStockPrice = 100.0;
    double strikePrice = 100.0;
    double riskFreeRate = 0.05;
    double volatility = 0.2;
    double timeToMaturity = 1.0;
    int    numberOfSimulations = 100000;

    // Pass seed=42 for reproducibility; set to 0 to randomize
    McPrice call = monteCarloOptionPricing(
        initialStockPrice, strikePrice, riskFreeRate,
        volatility, timeToMaturity,
        numberOfSimulations, true, 42);

    McPrice put = monteCarloOptionPricing(
        initialStockPrice, strikePrice, riskFreeRate,
        volatility, timeToMaturity,
        numberOfSimulations, false, 42);

    std::cout << "The Put fair value is " << put.fairValue
        << "  (95% CI +/-" << 1.96 * put.stdErr << ")  and expected payoff at maturity is "
        << put.expectedPayoff << '\n';

    std::cout << "The Call fair value is " << call.fairValue
        << "  (95% CI +/-" << 1.96 * call.stdErr << ")  and expected payoff at maturity is "
        << call.expectedPayoff << '\n';

    // Compare Monte Carlo call to closed-form Black–Scholes value
    double bsCall = blackScholesCall(
        initialStockPrice, strikePrice,
        riskFreeRate, volatility, timeToMaturity);

    std::cout << "Black–Scholes analytic call = " << bsCall
        << "  |  relative error = "
        << std::fabs(call.fairValue - bsCall) / bsCall * 100 << " %\n";

    return 0;
}

//S0 stock price at time 0 (now)
//ST stock price at time T (future)/ aka = maturity
//Maturity(T) how long until option expires
//Stike Price(K) fixed price which you can call or put the asset at maturity
// European option can only be exercised at maturity
// risk free rate(r) rate of return on risk free investment
// volitivity(sigma/o) high volitivity = high risk low volitivity = low risk
//
// (Black–Scholes model) ST formula = S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)
