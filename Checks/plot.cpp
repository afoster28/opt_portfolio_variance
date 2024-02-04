#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <map>

// Gnuplot headers
#include "gnuplot-iostream.h"

int main() {
    // Open the CSV file
    std::ifstream inputFile("your_file.csv");
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return 1;
    }

    // Skip the header
    std::string header;
    std::getline(inputFile, header);

    // Read values from the first column
    std::vector<double> values;
    double value;
    while (inputFile >> value) {
        values.push_back(value);
        char comma;
        if (inputFile >> comma && comma != ',') {
            std::cerr << "Error: CSV format mismatch." << std::endl;
            return 1;
        }
    }

    // Gnuplot
    Gnuplot gp;
    gp << "set boxwidth 0.5 relative\n";
    gp << "set style fill solid\n";
    gp << "plot '-' using 1:2 with boxes title 'Histogram'\n";

    // Count occurrences of each value
    std::map<double, int> histogram;
    for (const auto &val : values) {
        histogram[val]++;
    }

    // Send data to gnuplot
    for (const auto &entry : histogram) {
        gp << entry.first << " " << entry.second << "\n";
    }

    gp << "e\n";

    return 0;
}
