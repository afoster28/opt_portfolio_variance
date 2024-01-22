#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

class FinancialModel {
private:
    double mu;  // Annual mean return
    double sigma;  // Annual standard deviation
    std::default_random_engine generator;
    std::lognormal_distribution<double> distribution;
    std::vector<double> logReturns;

public:
    FinancialModel(double meanReturn, double stddev)
        : mu(meanReturn), sigma(stddev) {
        // Calculate mu and sigma for the lognormal distribution
        double muLog = std::log(1 + mu) - 0.5 * std::log(1 + (sigma * sigma) / ((1 + mu) * (1 + mu)));
        double sigmaLog = std::sqrt(std::log(1 + (sigma * sigma) / ((1 + mu) * (1 + mu))));

        // Set up the lognormal distribution
        distribution = std::lognormal_distribution<double>(muLog, sigmaLog);

        // Seed the random number generator
        std::random_device rd;
        generator.seed(rd());
    }

    // Function to generate annual return
    double generateAnnualReturn() {
        return distribution(generator);
    }

    // Function to generate a list of annual returns
    std::vector<double> generateAnnualReturnList(int numReturns) {
        std::vector<double> returns;
        returns.reserve(numReturns);

        for (int i = 0; i < numReturns; ++i) {
            double annualReturn = generateAnnualReturn();
            returns.push_back(annualReturn);
        }

        return returns;
    }

    // Function to process annual returns and calculate log returns
    void processReturns(const std::vector<double>& annualReturns) {
        logReturns.clear();
        logReturns.reserve(annualReturns.size());

        for (const auto& annualReturn : annualReturns) {
            // double logReturn = std::log(1 + annualReturn);
            double logReturn = std::log(annualReturn);
            logReturns.push_back(logReturn);
        }
    }

    // Function to calculate prices from log returns and initial spot price
    std::vector<double> calculatePrices(double spotPrice) const {
        std::vector<double> prices;
        // prices.reserve(logReturns.size() + 1);
        prices.reserve(logReturns.size());

        // double currentPrice = spotPrice;

        // The first price is the spot price
        // prices.push_back(currentPrice);

        // Calculate subsequent prices based on log returns
        // for (const auto& logReturn : logReturns) {
        //     currentPrice *= std::exp(logReturn);
        //     prices.push_back(currentPrice);
        // }
        
        // Calculate subsequent prices based on log returns
        for (const auto& logReturn : logReturns) {
            double scenarioPrice = spotPrice * (1 + logReturn);
            prices.push_back(scenarioPrice);
        }

        return prices;
    }

    // Function to save data to a CSV file
    void saveToCSV(const std::vector<double>& data, const std::string& filename) {
        std::ofstream outputFile(filename);

        if (outputFile.is_open()) {
            for (const auto& value : data) {
                outputFile << value << "\n";
            }
            outputFile.close();
            std::cout << "Data has been saved to: " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    // Function to run the financial model
    void runModel(int numReturns, double spotPrice) {
        // Generate annual returns
        std::vector<double> returnList = generateAnnualReturnList(numReturns);

        // Save annual returns to a CSV file
        saveToCSV(returnList, "annual_returns.csv");

        // Process annual returns to calculate log returns
        processReturns(returnList);

        // Calculate prices based on log returns
        std::vector<double> priceList = calculatePrices(spotPrice);

        // Save log returns to a CSV file
        saveToCSV(logReturns, "log_returns.csv");

        // Save prices to a CSV file
        saveToCSV(priceList, "asset_prices.csv");
    }
};

int main() {
    // Set asset return parameters
    double meanReturn = 0.05;      // 5% annual mean return
    double stdDev = 0.2;           // 20% annual standard deviation
    double spotPrice = 5000;       // typical S&P500 spot price

    // Create FinancialModel object
    FinancialModel financialModel(meanReturn, stdDev);

    // Run the financial model
    financialModel.runModel(10000, spotPrice);

    return 0;
}
