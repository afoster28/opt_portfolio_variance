#include <iostream>
#include <cmath>
#include <random>

class BlackScholes {
private:
    std::string dir;
    std::string optType;
    double S;  // Spot price of the underlying asset
    double K;  // Strike price of the option
    double r;  // Risk-free interest rate
    double sigma;  // Volatility of the underlying asset
    double T;  // Time to expiration of the option
    double Q; // Trade volume

public:
    // Constructor to initialize the parameters
    BlackScholes(std::string direction, std::string optionType, double spotPrice, double strikePrice, double interestRate, double volatility, double timeToExpiration, double quantity)
        : dir(direction), optType(optionType), S(spotPrice), K(strikePrice), r(interestRate), sigma(volatility), T(timeToExpiration), Q(quantity) {}

    // Method to calculate the Black-Scholes option price
    double calculateOptionPrice() {
        T = T / 365;

        double d1 = (log(S / K) + (r + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);

        double dir_num;

        if (dir == "L") {
            dir_num = 1;
        } else if (dir == "S") {
            dir_num = -1;
        }

        if (optType == "C") {
            return Q * dir_num * (S * normalDistribution(d1) - K * exp(-r * T) * normalDistribution(d2));
        } else if (optType == "P") {
            return Q * dir_num * (K * exp(-r * T) * normalDistribution(-d2) - S * normalDistribution(-d1));
        }

        // Default case (should not reach here)
        return 0.0;
    }

private:
    // Helper method to calculate the cumulative distribution function of the standard normal distribution
    double normalDistribution(double x) {
        static const double invSqrt2Pi = 0.3989422804014337;  // 1 / sqrt(2 * pi)
        return 0.5 * (1.0 + erf(x * invSqrt2Pi));
    }
};

int main() {
    // Example usage
    std::string direction = "L";
    std::string optionType = "C";
    double spotPrice = 5000.0;
    double strikePrice = 5000.0;
    double interestRate = 0.05;
    double volatility = 0.2;
    double timeToExpiration = 60.0;
    double quantity = 10.0;

    BlackScholes option(direction, optionType, spotPrice, strikePrice, interestRate, volatility, timeToExpiration, quantity);

    // Calculate and print call option price
    double callOptionPrice = option.calculateOptionPrice();
    std::cout << "Black-Scholes Option Price: " << callOptionPrice << std::endl;

    return 0;
}
