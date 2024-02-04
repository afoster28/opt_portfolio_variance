// Option Portfolio Variance Calculator
// Authors: Adam Foster, Sercan Akbay
//
// The program feeds in a portfolio of options, samples underlying asset prices at maturity, 
// calculates payoffs, portfolio payoffs and summary statistics per portfolio and persists
// final data and reference data

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

class Option
{
public:
    std::string type;      // Call or put
    std::string direction; // Long or short
    double strike;         // Strike price
    double timeToMaturity; // Time to maturity in days
    int quantity;          // Volume traded

    Option(const std::string &t, const std::string &d, double s, double ttm, int q)
        : type(t), direction(d), strike(s), timeToMaturity(ttm), quantity(q) {}
};

class Model
{
private:
    double spot;                                   // Spot price
    double mu;                                     // Annual mean return
    double sigma;                                  // Annual standard deviation
    std::default_random_engine generator;          // Random number generator for price sampling
    std::normal_distribution<double> distribution; // Normal distribution for log prices
    std::vector<Option> optionsList;               // List of options
    double ttm_annual;                             // Annualised time to maturity

public:
    Model(double spotPrice, double meanReturn, double stddev, const std::string &filename)
        : spot(spotPrice), mu(meanReturn), sigma(stddev)
    {
        // Load options from input csv
        if (!loadOptionsFromCSV(filename)) {
            throw std::runtime_error("Unable to open file: " + filename);
        }

        // Calculate mu and sigma sq of log price for the normal distribution
        double muLog = std::log(spot) + (mu - 0.5 * sigma * sigma) * ttm_annual;
        double sigmaLog = std::sqrt(sigma * sigma * ttm_annual);

        // Set up the normal distribution
        distribution = std::normal_distribution<double>(muLog, sigmaLog);

        // Seed the random number generator
        std::random_device rd;
        generator.seed(rd());
    }

    // Split a string based on a delimiter for CSV input routine
    std::vector<std::string> splitString(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }

    // Load options from a CSV file
    bool loadOptionsFromCSV(const std::string& filename) {
        std::ifstream inputFile(filename);
        if (!inputFile.is_open()) {
            throw std::runtime_error("Unable to open file: " + filename);
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(inputFile, line)) {
            ++lineNumber;
            std::cout << "Processing line " << lineNumber << ": " << line << std::endl;

            std::vector<std::string> optionData = splitString(line, ',');
            if (optionData.size() == 5) {
                // Modify the first element of the first option
                if (lineNumber == 1 && !optionData[0].empty()) {
                    optionData[0] = optionData[0].back();
                    ttm_annual = std::stod(optionData[3]) / 365;
                    std::cout << "Time to maturity: " << ttm_annual << std::endl;
                }

                Option option(optionData[0], optionData[1], std::stod(optionData[2]),
                            std::stod(optionData[3]), std::stoi(optionData[4]));
                optionsList.push_back(option);
            } else {
                std::cerr << "Error reading line " << lineNumber << ": " << line << std::endl;
            }
        }
        inputFile.close();

        std::cout << "Number of options loaded: " << optionsList.size() << std::endl;

        return true;
    }

    // Generate a list of prices
    double generatePrice() {
        return distribution(generator);
    }

    std::vector<double> generatePrices(int numPrices) {
        std::vector<double> sampledPrices;
        sampledPrices.reserve(numPrices);

        for (int i = 0; i < numPrices; ++i) {
            double sampledPrice = std::exp(generatePrice());
            sampledPrices.push_back(sampledPrice);
        }

        return sampledPrices;
    }

    // Calculate option payouts matrix
    std::vector<std::vector<double>> calculateOptionPayouts(const std::vector<double>& prices, size_t numOptions) const {
        std::vector<std::vector<double>> payouts;
        payouts.reserve(prices.size());

        for (size_t i = 0; i < prices.size(); ++i) { // Run for all prices
            std::vector<double> rowPayouts(numOptions, 0.0);

            for (size_t j = 0; j < numOptions; ++j) { // Consider all options up to the a specified index, needed for sequencing
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

            payouts.push_back(rowPayouts);
        }

        return payouts;
    }

    // Save combined data (price sample and option payouts) to a CSV file for reference
    void writeCombined(const std::vector<double>& assetPrices,
                const std::vector<std::vector<double>>& optionPayouts,
                const std::vector<Option>& optionsList,
                const std::string& filename, size_t numOptions) {
        std::ofstream outputFile(filename);

        if (outputFile.is_open()) {
            // Column headers: price and multiple options
            outputFile << "asset_prices";

            for (size_t j = 0; j < numOptions; ++j) {
                const auto& option = optionsList[j];
                outputFile << ',' << option.type << '|' << option.direction << '|' << option.strike
                        << '|' << option.timeToMaturity << '|' << option.quantity;
            }
            outputFile << "\n";

            // Write data rows
            for (size_t i = 0; i < assetPrices.size(); ++i) {
                // Write asset prices
                outputFile << std::fixed << std::setprecision(2) << assetPrices[i];

                // Write option payouts for the relevant options
                for (size_t j = 0; j < numOptions; ++j) {
                    outputFile << ',' << optionPayouts[i][j];
                }

                outputFile << "\n";
            }

            outputFile.close();
            std::cout << "Data has been saved to: " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    // Save result matrix to a CSV file to summarise final trade count, mean, variance, standard deviation
    void writeStats(const std::vector<std::vector<double>> &resultMatrix, const std::string &filename)
    {
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            // Write column headers
            outputFile << "trades,mean,variance,standard_deviation\n";

            // Write data rows
            for (size_t i = 0; i < resultMatrix.size(); ++i)
            {
                for (size_t j = 0; j < resultMatrix[i].size(); ++j)
                {
                    outputFile << std::fixed << std::setprecision(6) << resultMatrix[i][j];
                    if (j < resultMatrix[i].size() - 1)
                    {
                        outputFile << ',';
                    }
                }
                outputFile << "\n";
            }

            outputFile.close();
            std::cout << "Result matrix has been saved to: " << filename << std::endl;
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    // Add initial zeroes row to result matrix ahead of any trade processing
    std::vector<std::vector<double>> resultMatrixAddZeroes(std::vector<std::vector<double>> resultMatrixInput) {
        std::vector<double> initialRow{0, 0, 0, 0};
        resultMatrixInput.push_back(initialRow);

        return resultMatrixInput;
    }

    void runModel(int numPrices, double spotPrice)
    {
        // Price list
        std::vector<double> priceList = generatePrices(numPrices);

        // Calculate and persist trade count, mean, variance, standard deviation and reference table
        std::vector<std::vector<double>> resultMatrix = calculatePortfolioStatistics(optionsList.size(), priceList);
    }

private:
    // Validate user input as 'complete' or 'sequential'
    bool isValidRunType(const std::string& runType) const {
        return (runType == "complete" || runType == "sequential");
    }

    // Validate user input as a non-negative double
    double isValidVarianceCutoff(const std::string& prompt) const {
        double value;
        bool validInput = false;

        do {
            std::cout << prompt;
            std::cin >> value;

            if (std::cin.fail() || value < 0.0) {
                std::cout << "Invalid input. Please enter a non-negative number." << std::endl;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            } else {
                validInput = true;
            }
        } while (!validInput);

        return value;
    }
    
    // Calculate and persist trade count, mean, variance, standard deviation and reference table
    std::vector<std::vector<double>> calculatePortfolioStatistics(size_t numOptions, const std::vector<double> &priceList) {
        std::vector<std::vector<double>> resultMatrix;

        // Prompt user for 'complete' or 'sequential' run
        std::string runType;
        do {
            std::cout << "Enter 'complete' for complete run or 'sequential' for sequential run: ";
            std::cin >> runType;
        } while (!isValidRunType(runType));

        double varianceCutoff = 0.0;
        bool resetTrigger = false;

        resultMatrix = resultMatrixAddZeroes(resultMatrix);

        // Print the initial row of statistics
        std::cout << "Option Count: " << 0
                << ", Mean: " << 0
                << ", Variance: " << 0
                << ", Standard Deviation: " << 0 << std::endl;

        for (size_t i = 0; i < numOptions; ++i) {
            // Reset feature
            if (resetTrigger == true) {
                i = 0;
                resultMatrix.clear();
                resultMatrix = resultMatrixAddZeroes(resultMatrix);
                resetTrigger = false;
            }
            
            // Prompt the user for a variance cutoff in the sequential run
            if (runType == "sequential") {
                if (i == 0) {
                    varianceCutoff = isValidVarianceCutoff("Please input variance cutoff: ");
                }
                
                char continueOption;
                std::cout << "Do you want to continue with the next option? (y/n/r): ";
                std::cin >> continueOption;

                if (continueOption == 'r') {
                    resetTrigger = true;
                }
                else if (continueOption == 'n') {
                    // Break the loop if the user chooses not to continue
                    std::vector<Option> partialOptionsOut(optionsList.begin(), optionsList.begin() + i);
                    std::vector<std::vector<double>> optionPayoutsOut = calculateOptionPayouts(priceList, numOptions);
                    
                    // Save combined data to CSV if the user chooses not to continue
                    writeCombined(priceList, optionPayoutsOut, partialOptionsOut, "combined_table.csv", partialOptionsOut.size());
                    break;
                }
            }

            if (resetTrigger == false ) {
                std::vector<Option> partialOptions(optionsList.begin(), optionsList.begin() + i + 1);

                // Calculate option payouts with the specified number of options
                std::vector<std::vector<double>> optionPayouts = calculateOptionPayouts(priceList, numOptions);

                // Accumulate payouts across options up to and including the current option
                std::vector<double> portfolioPayouts(priceList.size(), 0.0);
                for (size_t j = 0; j <= i; ++j) {
                    for (size_t k = 0; k < priceList.size(); ++k) {
                        portfolioPayouts[k] += optionPayouts[k][j]; // Accumulate payouts for the current portfolio
                    }
                }

                // Calculate mean
                double mean = 0.0;
                for (size_t k = 0; k < priceList.size(); ++k) {
                    mean += portfolioPayouts[k] / priceList.size();
                }

                // Calculate variance
                double varianceSum = 0.0;
                for (size_t k = 0; k < priceList.size(); ++k) {
                    double sum = portfolioPayouts[k];
                    varianceSum += (sum - mean) * (sum - mean);
                }
                double variance = varianceSum / priceList.size();

                // Calculate standard deviation
                double stdDeviation = std::sqrt(variance);

                // Store the results in a row
                std::vector<double> resultRow{static_cast<double>(partialOptions.size()), mean, variance, stdDeviation};
                resultMatrix.push_back(resultRow);

                // Print the row of statistics
                std::cout << "Option Count: " << partialOptions.size()
                        << ", Mean: " << mean
                        << ", Variance: " << variance
                        << ", Standard Deviation: " << stdDeviation << std::endl;

                if (runType == "sequential" && variance > varianceCutoff){
                    writeCombined(priceList, optionPayouts, partialOptions, "combined_table.csv", partialOptions.size());
                    std::cout << "Variance cutoff exceed, ending program" << std::endl;
                    break;
                }

                // Save combined data to CSV in the last iteration
                if (i == numOptions - 1) {
                    writeCombined(priceList, optionPayouts, partialOptions, "combined_table.csv", partialOptions.size());
                }
            }
        }

        // Persist the statistics matrix to CSV
        writeStats(resultMatrix, "result_matrix.csv");

        return resultMatrix;
    }
};

int main()
{
    try {
        // Input filename from the user
        std::string inputFilename;
        std::cout << "Enter the name of the input CSV file: ";
        std::cin >> inputFilename;

        // Asset parameters
        double meanReturn = 0.05; // 5% annual mean return
        double stdDev = 0.2;      // 20% annual standard deviation
        double spotPrice = 5000;  // Typical S&P500 spot price

        // Create Model object, loading input data and defining inputs for price sampling
        Model portfolioModel(spotPrice, meanReturn, stdDev, inputFilename);

        // Run the model, generating prices, option payouts and summary statistics per incremental portfolio
        portfolioModel.runModel(10000, spotPrice);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1; // Terminate the program with an error code, e.g. if cannot load input file
    }

    return 0;
}
