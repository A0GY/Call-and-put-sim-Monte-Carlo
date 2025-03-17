// A program to simulate strike prices of a stock using Monte Carlo simulation
// we output a call and put depending on factores such as the stock price, strike price, risk free rate, volatility, time to maturity, and number of simulations
// to then give a value of the option
// this is a basic example of how to use Monte Carlo simulation to price options desgin in cpp
// Author : A Y




#include <iostream>     
#include <cmath>       
#include <random>       
#include <vector>  
#include <utility> 

// Generates a single random draw from a Normal(mean, stddev) distribution
double generateGaussianNoise(double mean, double stddev) {
    static std::mt19937 generator(std::random_device{}());  // Mersenne Twister engine (seeded once)
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

// Monte Carlo simulation for pricing a European option (Call or Put)
std::pair<double,double> monteCarloOptionPricing(double initialStockPrice,double strikePrice,double riskFreeRate,double volatility,double timeToMaturity,int numberOfSimulations,bool isCallOption)
{
    double totalPayoff = 0.0;

    for (int i = 0; i < numberOfSimulations; ++i) {
        // Z is a standard normal random variable (mean=0, stddev=1)
        double Z = generateGaussianNoise(0.0, 1.0);

        // Calculate the stock price at maturity using the Black–Scholes model assumption:
        //     S_T = S_0 * exp( (r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z )
        double stockPriceAtMaturity = initialStockPrice * std::exp((riskFreeRate - 0.5 * volatility * volatility) * timeToMaturity + volatility * std::sqrt(timeToMaturity) * Z);
		// ST the stock price predicted at maturity over 10,000 simulations
       
        
        
        // Determine payoff depending on whether it's a Call or Put
        double payoff = 0.0;
        if (isCallOption) {
            payoff = getCallPayoff(stockPriceAtMaturity, strikePrice);
        }
        else {
            payoff = getPutPayoff(stockPriceAtMaturity, strikePrice);
        }

        // Accumulate sum of all payoffs
        totalPayoff += payoff;
    }

    // Average payoff across all simulations
    double averagePayoff = totalPayoff / static_cast<double>(numberOfSimulations);

    // Discount the average payoff back to present value using riskFreeRate
    double discountedPayoff = std::exp(-riskFreeRate * timeToMaturity) * averagePayoff;
	double notDiscountedPayoff = averagePayoff;
    return { discountedPayoff,notDiscountedPayoff };
}

int main() {
    // Option parameters
    double initialStockPrice = 100.0;
    double strikePrice = 100.0;
    double riskFreeRate = 0.05;
    double volatility = 0.2;
    double timeToMaturity = 1.0;
    int numberOfSimulations = 100000;

    // Compute Call Option Prices
    std::pair<double, double> callOptionPriceboth = monteCarloOptionPricing(
        initialStockPrice, strikePrice, riskFreeRate, volatility, timeToMaturity, numberOfSimulations, true
    );

    double callOptionPrice = callOptionPriceboth.first;   // Discounted Call Price
    double callOptionPricedis = callOptionPriceboth.second; // Non-discounted Call Price

    // Compute Put Option Prices
    std::pair<double, double> putOptionPriceBoth = monteCarloOptionPricing(
        initialStockPrice, strikePrice, riskFreeRate, volatility, timeToMaturity, numberOfSimulations, false
    );

    double putOptionPrice = putOptionPriceBoth.first;   // Discounted Put Price
    double putOptionPricedis = putOptionPriceBoth.second; // Non-discounted Put Price

	std::cout << "The Put fair value is " << putOptionPrice << " and expected payoff at maturtiy is " << putOptionPricedis << std::endl;
	std::cout << "The call fair vale is " << callOptionPrice << " and expected payoff at maturuity is  " << callOptionPricedis << std::endl;
	// using std::pair to give the end user fiar value and expected payoff by taking in a extra result in monte carlo function
	// which removes the discounting of the payoff calculations to give the expected payoff at maturity

    return 0;
}

//S0 stock price at time 0 (now)
//ST stock price at time T (future)/ aka = maturity
//Maturity(T) how long until option expires
//Stike Price(K) fixed price which you can call or put the asset at maturity
// European option can only be exercised at maturity
// risk free rate(r) rate of return on risk free investment
// volitivity(sigma/o) high volitivity = high risk low volitivity = low risk


// (Black–Scholes model) ST formula = S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)