#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip> // For std::setprecision
#include <regex> // For std::regex

struct Matrix {
    std::vector<std::vector<double>> data;
};

// Function to read matrix from CSV file
Matrix readMatrixFromCSV(const std::string& filename) {
    Matrix matrix;
    std::ifstream file(filename);
    std::string line;

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }

        matrix.data.push_back(row);
    }

    return matrix;
}

// Adjusted function to handle potential parsing issues for TXT files
Matrix convertTXTToMatrix(const std::string& filename) {
    Matrix matrix;
    std::ifstream file(filename);
    std::string line;
    // Skip the header
    std::getline(file, line);

    std::regex numberRegex(R"([\d\.]+)"); // Regular expression to match numbers

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream ss(line);
        std::string token;
        while (ss >> token) {
            // Use regex to check if token is a valid number
            if (std::regex_match(token, numberRegex)) {
                try {
                    row.push_back(std::stod(token));
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Conversion error: " << e.what() << " for token: " << token << '\n';
                    // Optionally handle the error, e.g., by continuing to the next token
                }
            }
        }

        if (!row.empty()) {
            matrix.data.push_back(row);
        }
    }

    return matrix;
}

// Function to calculate percentage difference matrix
Matrix calculatePercentageDifference(const Matrix& modelOutput, const Matrix& expectedOutput) {
    Matrix differenceMatrix;

    for (size_t i = 0; i < modelOutput.data.size(); ++i) {
        std::vector<double> row;
        for (size_t j = 0; j < modelOutput.data[i].size(); ++j) {
            double difference = 0.0;
            if (expectedOutput.data[i][j] != 0) { // Avoid division by zero
                difference = (modelOutput.data[i][j] - expectedOutput.data[i][j]) / expectedOutput.data[i][j];
            }
            row.push_back(difference);
        }
        differenceMatrix.data.push_back(row);
    }

    return differenceMatrix;
}

// Function to save matrix to CSV file
void saveMatrixToCSV(const Matrix& matrix, const std::string& filename) {
    std::ofstream file(filename);

    for (const auto& row : matrix.data) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << std::fixed << std::setprecision(2) << row[j];
            if (j < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
}

int main() {
    std::string modelOutputFilename, expectedOutputFilename;

    // Prompt user for filenames
    std::cout << "Enter model output file name: ";
    std::cin >> modelOutputFilename;
    std::cout << "Enter expected output file name: ";
    std::cin >> expectedOutputFilename;

    Matrix modelOutput = readMatrixFromCSV(modelOutputFilename);
    Matrix expectedOutput = convertTXTToMatrix(expectedOutputFilename);
    Matrix differenceMatrix = calculatePercentageDifference(modelOutput, expectedOutput);

    std::string outputFilename = "difference_matrix.csv";
    saveMatrixToCSV(differenceMatrix, outputFilename);

    std::cout << "Percentage difference matrix saved to " << outputFilename << std::endl;

return 0;
}