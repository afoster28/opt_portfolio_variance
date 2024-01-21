#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <limits>

class BlackScholes {
private:
    std::string dir;
    std::string optType;
    double S;  // Spot price
    double K;  // Strike price
    double r;  // Risk-free rate
    double sigma;  // Volatility
    double T;  // Time to expiration
    double Q;  // Trade volume
    bool debugMode;  // Debug mode flag

public:
    BlackScholes(std::string direction, std::string optionType, double spotPrice, double strikePrice, double interestRate, double volatility, double timeToExpiration, double quantity, bool debug = false)
        : dir(direction), optType(optionType), S(spotPrice), K(strikePrice), r(interestRate), sigma(volatility), T(timeToExpiration), Q(quantity), debugMode(debug) {}

    double normalDistribution(double x) const {
        if (std::abs(x) > 6) {  // Handling extreme values
            return (x > 0) ? 1.0 : 0.0;
        }
        static const double invSqrt2Pi = 0.3989422804014337;
        double nd = 0.5 * (1.0 + erf(x * invSqrt2Pi));
        return nd;
    }

    double calculateOptionPrice() const {
        double adjustedT = T / 365;  // Adjust T for calculation

        double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * adjustedT) / (sigma * sqrt(adjustedT));
        double d2 = d1 - sigma * sqrt(adjustedT);

        double term1 = S * normalDistribution(d1);
        double term2 = K * exp(-r * adjustedT) * normalDistribution(d2);

        double optionPrice = (optType == "C") ? term1 - term2 : term2 - term1;
        optionPrice *= ((dir == "L") ? 1 : -1) * Q;

        if (debugMode) {
            std::cout << "Calculated Option Price: " << optionPrice << std::endl;
        }
        return optionPrice;
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

    std::cout << "Enable debug mode? (y/n): ";
    char debugChoice;
    std::cin >> debugChoice;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


  bool debugMode = (debugChoice == 'y' || debugChoice == 'Y');

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

      input_stream >> optType >> delimiter >> dir >> delimiter >> K >> delimiter >> T >> delimiter >> Q;

      if (input_stream && delimiter == ',') {
          myPortfolio.addOption(BlackScholes(dir, optType, defaultSpotPrice, K, defaultInterestRate, defaultVolatility, T, Q, debugMode));
          std::cout << "Option added to portfolio." << std::endl;
      } else {
          std::cerr << "Invalid input. Please enter all parameters correctly, separated by commas." << std::endl;
      }
  }

  double portfolioVariance = myPortfolio.getVariance();
  std::cout << "Portfolio Variance: " << portfolioVariance << std::endl;

  return 0;
}
