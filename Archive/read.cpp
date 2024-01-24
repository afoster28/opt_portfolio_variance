#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// Function to split a string based on a delimiter
vector<string> splitString(const string& s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

int main() {
    // Open the CSV file
    ifstream inputFile("input.csv");
    if (!inputFile.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    // Read and display header
    string header;
    getline(inputFile, header);
    cout << "Header: " << header << endl;

    // Read and display options
    vector<vector<string>> options;
    string line;
    while (getline(inputFile, line)) {
        options.push_back(splitString(line, ','));
    }

    // Print available options
    cout << "Available Options:" << endl;
    for (int i = 0; i < options.size(); ++i) {
        cout << i + 1 << ": ";
        for (const string& element : options[i]) {
            cout << element << " ";
        }
        cout << endl;
    }

    // Get user input for the row to display
    int selectedRow;
    cout << "Enter the row number to display: ";
    cin >> selectedRow;

    // Validate user input
    if (selectedRow < 1 || selectedRow > options.size()) {
        cerr << "Invalid row number!" << endl;
        return 1;
    }

    // Display the selected row
    cout << "Selected Row: ";
    for (const string& element : options[selectedRow - 1]) {
        cout << element << " ";
    }
    cout << endl;

    return 0;
}
