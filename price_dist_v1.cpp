#include <iostream>
#include <fstream>
#include <vector>
#include <random>

class AssetReturnsSampler {
private:
    double mu;  // Annual mean return
    double sigma;  // Annual standard deviation
    std::default_random_engine generator;
    std::lognormal_distribution<double> distribution;

public:
    AssetReturnsSampler(double meanReturn, double stddev)
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

    // Function to save the list to a CSV file
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
};

int main() {
    // Set asset return parameters
    double meanReturn = 0.05;      // 5% annual mean return
    double stdDev = 0.2;           // 20% annual standard deviation

    // Create AssetReturnsSampler object
    AssetReturnsSampler assetReturnsSampler(meanReturn, stdDev);

    // Generate a list of 100 annual returns
    std::vector<double> returnList = assetReturnsSampler.generateAnnualReturnList(10000);

    // Save the generated annual returns to a CSV file
    assetReturnsSampler.saveToCSV(returnList, "annual_returns.csv");

    return 0;
}
