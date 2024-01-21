#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <numeric>

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
    BlackScholes(std::string direction, std::string optionType, double spotPrice, double strikePrice, double interestRate, double volatility, double timeToExpiration, double quantity)
        : dir(direction), optType(optionType), S(spotPrice), K(strikePrice), r(interestRate), sigma(volatility), T(timeToExpiration), Q(quantity) {}

    double calculateOptionPrice() const {
        double adjustedT = T / 365;  // Adjust T for calculation

        double d1 = (log(S / K) + (r + 0.5 * pow(sigma, 2)) * adjustedT) / (sigma * sqrt(adjustedT));
        double d2 = d1 - sigma * sqrt(adjustedT);

        double dir_num = (dir == "L") ? 1 : -1;

        if (optType == "C") {
            return Q * dir_num * (S * normalDistribution(d1) - K * exp(-r * adjustedT) * normalDistribution(d2));
        } else if (optType == "P") {
            return Q * dir_num * (K * exp(-r * adjustedT) * normalDistribution(-d2) - S * normalDistribution(-d1));
        }

        return 0.0;
    }

private:
    double normalDistribution(double x) const {
        static const double invSqrt2Pi = 0.3989422804014337;
        return 0.5 * (1.0 + erf(x * invSqrt2Pi));
    }
};

class Portfolio {
private:
    std::vector<BlackScholes> options;

public:
    void addOption(const BlackScholes& option) {
        options.push_back(option);
    }

    void resetPortfolio() {
        options.clear();
    }

    double getVariance() const {
        if (options.empty()) return 0.0;

        double mean = std::accumulate(options.begin(), options.end(), 0.0, 
                                      [](double sum, const BlackScholes& op) { return sum + op.calculateOptionPrice(); }) 
                      / options.size();

        double variance = std::accumulate(options.begin(), options.end(), 0.0, 
                                          [mean](double sum, const BlackScholes& op) {
                                              double diff = op.calculateOptionPrice() - mean;
                                              return sum + diff * diff;
                                          }) 
                         / options.size();

        return variance;
    }
};

int main() {
    Portfolio myPortfolio;
    std::string userInput;

    // Default values for S, r, and sigma
    double defaultSpotPrice = 5000.0;
    double defaultInterestRate = 0.05;
    double defaultVolatility = 0.2;

    while (true) {
        std::cout << "Enter option parameters (optType, dir, K, T, Q), or 'exit' to finish: ";
        std::getline(std::cin, userInput);

        if (userInput == "exit") {
            break;
        }

        std::istringstream input_stream(userInput);
        std::string dir, optType;
        double K, T, Q;
        char delimiter;

        input_stream >> optType;
        input_stream >> delimiter;
        input_stream >> dir;
        input_stream >> delimiter;
        input_stream >> K;
        input_stream >> delimiter;
        input_stream >> T;
        input_stream >> delimiter;
        input_stream >> Q;

        if (input_stream && delimiter == ',') {
            myPortfolio.addOption(BlackScholes(dir, optType, defaultSpotPrice, K, defaultInterestRate, defaultVolatility, T, Q));
            std::cout << "Option added to portfolio." << std::endl;
        } else {
            std::cerr << "Invalid input. Please enter all parameters correctly, separated by commas." << std::endl;
        }
    }

    double portfolioVariance = myPortfolio.getVariance();
    std::cout << "Portfolio Variance: " << portfolioVariance << std::endl;

    return 0;
}
