#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>
#include <algorithm> // Include for std::remove_if
#include <cctype>    // Include for std::isspace

class Option
{
public:
    std::string type;      // Call or Put
    std::string direction; // Long or Short
    double strike;
    double timeToMaturity;
    int quantity;

    Option(const std::string &t, const std::string &d, double s, double ttm, int q)
        : type(t), direction(d), strike(s), timeToMaturity(ttm), quantity(q) {}
};

class FinancialModel
{
private:
    double spot;
    double mu; // Annual mean return
    double sigma; // Annual standard deviation
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    std::vector<Option> optionsList; // List of options

public:
    FinancialModel(double spotPrice, double meanReturn, double stddev)
        : spot(spotPrice), mu(meanReturn), sigma(stddev)
    {
        // Calculate mu and sigma sq of log price for the normal distribution
        double muLog = std::log(spot) + (mu - 0.5 * sigma * sigma) / 6;
        double sigmaLog = std::sqrt(sigma * sigma / 6);

        // Set up the normal distribution
        distribution = std::normal_distribution<double>(muLog, sigmaLog);

        // Seed the random number generator
        std::random_device rd;
        generator.seed(rd());

        // Load options from input csv
        loadOptionsFromCSV("sample_small.csv");
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

    // Function to generate price
    double generatePrice() {
        return distribution(generator);
    }

    // Function to generate a list of prices
    std::vector<double> generatePrices(int numPrices) {
        std::vector<double> sampledPrices;
        sampledPrices.reserve(numPrices);

        for (int i = 0; i < numPrices; ++i) {
            double sampledPrice = std::exp(generatePrice());
            sampledPrices.push_back(sampledPrice);
        }

        return sampledPrices;
    }

    // Function to calculate option payouts
    std::vector<std::vector<double>> calculateOptionPayouts(const std::vector<double>& prices) const {
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
    void writeCombined(const std::vector<double>& assetPrices,
                const std::vector<std::vector<double>>& optionPayouts,
                const std::vector<Option>& optionsList,
                const std::string& filename) {
        std::ofstream outputFile(filename);

        if (outputFile.is_open()) {
            // Write column headers
            outputFile << "asset_prices";

            // Write option_payouts headers
            for (const auto& option : optionsList) {
                outputFile << ',' << option.type << '|' << option.direction << '|' << option.strike
                        << '|' << option.timeToMaturity << '|' << option.quantity;
            }
            outputFile << "\n";

            // Write data rows
            for (size_t i = 0; i < assetPrices.size(); ++i) {
                // Write asset_prices
                outputFile << std::fixed << std::setprecision(2) << assetPrices[i];

                // Write option_payouts columns
                for (const auto& optionPayout : optionPayouts[i]) {
                    outputFile << ',' << optionPayout;
                }

                outputFile << "\n";
            }

            outputFile.close();
            std::cout << "Data has been saved to: " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    // Function to save result matrix to a CSV file
    void writeStats(const std::vector<std::vector<double>> &resultMatrix, const std::string &filename)
    {
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            // Write column headers
            outputFile << "OptionCount,Mean,Variance,StdDeviation\n";

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

    void runModel(int numPrices, double spotPrice)
    {
        // Generate prices (moved outside the loop)
        std::vector<double> priceList = generatePrices(numPrices);

        // Save combined data to CSV
        std::vector<std::vector<double>> optionPayouts = calculateOptionPayouts(priceList);
        writeCombined(priceList, optionPayouts, optionsList, "combined_table.csv");

        // Calculate and store option count, mean, variance, and standard deviation
        std::vector<std::vector<double>> resultMatrix = calculateOptionStatistics(optionsList, priceList);

        // Process and save the resultMatrix
        processResultMatrix(resultMatrix);
    }

private:
    // Function to calculate and store option count, mean, variance, and standard deviation
    std::vector<std::vector<double>> calculateOptionStatistics(const std::vector<Option> &optionsList,
                                                              const std::vector<double> &priceList)
    {
        std::vector<std::vector<double>> resultMatrix;

        for (size_t i = 0; i < optionsList.size(); ++i)
        {
            std::vector<Option> partialOptions(optionsList.begin(), optionsList.begin() + i + 1);

            // Calculate option payouts
            std::vector<std::vector<double>> optionPayouts = calculateOptionPayouts(priceList);

            // Accumulate payoffs across options up to and including the current option
            std::vector<double> portfolioPayoffs(priceList.size(), 0.0);
            for (size_t j = 0; j <= i; ++j)
            {
                for (size_t k = 0; k < priceList.size(); ++k)
                {
                    portfolioPayoffs[k] += optionPayouts[k][j]; // Accumulate payoffs for the current option
                }
            }

            // Calculate mean
            double mean = 0.0;
            for (size_t k = 0; k < priceList.size(); ++k)
            {
                mean += portfolioPayoffs[k] / priceList.size();
            }

            // Calculate variance
            double varianceSum = 0.0;
            for (size_t k = 0; k < priceList.size(); ++k)
            {
                double sum = portfolioPayoffs[k];
                varianceSum += (sum - mean) * (sum - mean);
            }
            double variance = varianceSum / priceList.size();

            // Calculate standard deviation
            double stdDeviation = std::sqrt(variance);

            // Store the results in a row
            std::vector<double> resultRow{static_cast<double>(partialOptions.size()), mean, variance, stdDeviation};
            resultMatrix.push_back(resultRow);
        }

        return resultMatrix;
    }

    // Function to process and print statistics for resultMatrix
    void processResultMatrix(const std::vector<std::vector<double>> &resultMatrix)
    {
        // Debug print for the resultMatrix
        std::cout << "Result Matrix:" << std::endl;
        for (size_t i = 0; i < resultMatrix.size(); ++i)
        {
            // Extract values for readability
            size_t optionCount = static_cast<size_t>(resultMatrix[i][0]);
            double mean = resultMatrix[i][1];
            double variance = resultMatrix[i][2];
            double stdDeviation = resultMatrix[i][3];

            // Print option count, mean, variance, and standard deviation on one line
            std::cout << "Option Count: " << optionCount
                    << ", Mean: " << mean
                    << ", Variance: " << variance
                    << ", Standard Deviation: " << stdDeviation << std::endl;
        }

        // Save results to CSV
        writeStats(resultMatrix, "result_matrix.csv");
    }
};

int main()
{
    // Set asset parameters
    double meanReturn = 0.05; // 5% annual mean return
    double stdDev = 0.2;      // 20% annual standard deviation
    double spotPrice = 5000;  // typical S&P500 spot price

    // Create FinancialModel object
    FinancialModel financialModel(spotPrice, meanReturn, stdDev);

    // Run the financial model
    financialModel.runModel(10000, spotPrice);

    return 0;
}