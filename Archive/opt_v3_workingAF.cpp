#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>
#include <algorithm> // Include for std::remove_if
#include <cctype>  // Include for std::isspace

class Option {
public:
    std::string type;      // Call or Put
    std::string direction; // Long or Short
    double strike;
    double timeToMaturity;
    int quantity;

    Option(const std::string& t, const std::string& d, double s, double ttm, int q)
        : type(t), direction(d), strike(s), timeToMaturity(ttm), quantity(q) {}
};

class FinancialModel {
private:
    double mu;  // Annual mean return
    double sigma;  // Annual standard deviation
    std::default_random_engine generator;
    std::lognormal_distribution<double> distribution;
    std::vector<double> logReturns;
    std::vector<Option> optionsList;  // List of options

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

        // Load options from input.csv
        loadOptionsFromCSV("input.csv");
    }

    // Function to split a string based on a delimiter
    std::vector<std::string> splitString(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }

    // Function to load options from a CSV file
    void loadOptionsFromCSV(const std::string& filename) {
        std::ifstream inputFile(filename);
        if (!inputFile.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            return;
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(inputFile, line)) {
            ++lineNumber;

            // Add this line for additional debugging
            std::cout << "Processing line " << lineNumber << ": " << line << std::endl;

            std::vector<std::string> optionData = splitString(line, ',');
            if (optionData.size() == 5) {
                // Modify the first element of the first option
                if (lineNumber == 1 && !optionData[0].empty()) {
                    optionData[0] = optionData[0].back();
                }

                Option option(optionData[0], optionData[1], std::stod(optionData[2]),
                            std::stod(optionData[3]), std::stoi(optionData[4]));
                optionsList.push_back(option);
            } else {
                std::cerr << "Error reading line " << lineNumber << ": " << line << std::endl;
            }
        }
        inputFile.close();

        // Debugging statement
        std::cout << "Number of options loaded: " << optionsList.size() << std::endl;
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
            double logReturn = std::log(annualReturn);
            logReturns.push_back(logReturn);
        }
    }

    // Function to calculate prices from log returns and initial spot price
    std::vector<double> calculatePrices(double spotPrice) const {
        std::vector<double> prices;
        prices.reserve(logReturns.size());

        for (const auto& logReturn : logReturns) {
            double scenarioPrice = spotPrice * (1 + logReturn);
            prices.push_back(scenarioPrice);
        }

        return prices;
    }

    // Function to calculate option payouts
    std::vector<std::vector<double>> calculateOptionPayouts(const std::vector<double>& prices) const {
        // Debugging print statement
        std::cout << "Number of options: " << optionsList.size() << std::endl;

        std::vector<std::vector<double>> payouts;
        payouts.reserve(prices.size());

        for (size_t i = 0; i < prices.size(); ++i) {
            std::vector<double> rowPayouts(optionsList.size(), 0.0);

            for (size_t j = 0; j < optionsList.size(); ++j) {
                const auto& option = optionsList[j];
                double optionPayout = 0.0;

                if (option.type == "C") {
                    optionPayout = std::max(0.0, prices[i] - option.strike);
                } else if (option.type == "P") {
                    optionPayout = std::max(0.0, option.strike - prices[i]);
                }

                if (option.direction == "S") {
                    optionPayout *= -1;
                }

                optionPayout *= option.quantity;
                rowPayouts[j] = optionPayout;
            }

            // Add rowPayouts to the payouts vector
            payouts.push_back(rowPayouts);
        }

        return payouts;
    }

    // Function to calculate the sum of values per price (row)
    std::vector<double> calculateRowSum(const std::vector<std::vector<double>>& payouts) const {
        std::vector<double> rowSums;
        rowSums.reserve(payouts.size());

        for (const auto& row : payouts) {
            double rowSum = 0.0;
            for (const auto& value : row) {
                rowSum += value;
            }
            rowSums.push_back(rowSum);
        }

        return rowSums;
    }

    // Function to save data to a CSV file
    void saveToCSV(const std::vector<double>& data, const std::string& filename) {
        std::ofstream outputFile(filename);

        if (outputFile.is_open()) {
            for (const auto& value : data) {
                outputFile << std::fixed << std::setprecision(2) << value << "\n";
            }
            outputFile.close();
            std::cout << "Data has been saved to: " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    // Function to save matrix data to a CSV file
    void saveToCSV(const std::vector<std::vector<double>>& data, const std::string& filename) {
        std::ofstream outputFile(filename);

        if (outputFile.is_open()) {
            for (const auto& row : data) {
                for (const auto& value : row) {
                    outputFile << std::fixed << std::setprecision(2) << value << ",";
                }
                outputFile << "\n";
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

        // Calculate option payouts
        std::vector<std::vector<double>> optionPayouts = calculateOptionPayouts(priceList);

        // Save option payouts to a CSV file
        saveToCSV(optionPayouts, "option_payouts.csv");

        // Calculate row sums
        std::vector<double> rowSums = calculateRowSum(optionPayouts);

        // Save row sums to a CSV file
        saveToCSV(rowSums, "row_sums.csv");

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
    financialModel.runModel(100, spotPrice);

    return 0;
}
