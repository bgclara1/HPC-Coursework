#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

string readLastLine(const string &filename) {
    ifstream infile(filename.c_str());
    if (!infile) {
        return "File not found: " + filename;
    }
    string line, lastLine;
    while (getline(infile, line)) {
        if (!line.empty())
            lastLine = line;
    }
    return lastLine;
}

void printPositions(const string &line) {
    
    stringstream ss(line);// split  line into tokens.
    vector<string> tokens;
    string token;
    while(ss >> token) {
        tokens.push_back(token);
    }
    
    cout << "Final Positions:" << endl;
    if (tokens.empty()) {
        cout << "  No data available" << endl;
        return;
    }
    
    cout << "  Timestamp: " << tokens[0] << endl;
    
    if (tokens.size() == 3) {
        cout << "  Particle 0: X = " << tokens[1] << ", Y = " << tokens[2] << endl;
    } else if (tokens.size() == 5) {
        cout << "  Particle 0: X = " << tokens[1] << ", Y = " << tokens[2] << endl;
        cout << "  Particle 1: X = " << tokens[3] << ", Y = " << tokens[4] << endl;
    } else {
        cout << "  ";
        for (const auto &tok : tokens) {
            cout << tok << " ";
        }
        cout << endl;
    }
}

void printEnergies(const string &line) {
    stringstream ss(line);
    vector<string> tokens;
    string token;
    while(ss >> token) {
        tokens.push_back(token);
    }
    
    cout << "Final Kinetic Energies:" << endl;
    if (tokens.empty()) {
        cout << "  No data available" << endl;
        return;
    }
    
    cout << "  Timestamp: " << tokens[0] << endl;
    
    if (tokens.size() == 2) {
        cout << "  Particle 0: KE = " << tokens[1] << endl;
    } else if (tokens.size() == 3) {
        cout << "  Particle 0: KE = " << tokens[1] << endl;
        cout << "  Particle 1: KE = " << tokens[2] << endl;
    } else {
        cout << "  ";
        for (const auto &tok : tokens) {
            cout << tok << " ";
        }
        cout << endl;
    }
}

int main() {
    cout << "Compiling serialSim.cpp..." << endl;

    vector<string> testCases = {
        "./serialSim --ic-one",
        "./serialSim --ic-one-vel",
        "./serialSim --ic-two",
        "./serialSim --ic-two-pass1",
        "./serialSim --ic-two-pass2",
        "./serialSim --ic-two-pass3"
    };
    
    for (size_t i = 0; i < testCases.size(); i++) {
        cout << "\nRunning test case " << i + 1 << ": " << testCases[i] << endl;
        int ret = system(testCases[i].c_str());
        if (ret != 0) {
            cout << "Test case " << i + 1 << " returned error code " << ret << endl;
        } else {
            cout << "Test case " << i + 1 << " executed successfully." << endl;
        }
        
        string finalPositions = readLastLine("positions.txt");
        string finalEnergies = readLastLine("energy.txt");
        
        cout << "\n----------------------------------" << endl;
        printPositions(finalPositions);
        cout << endl;
        printEnergies(finalEnergies);
        cout << "----------------------------------" << endl;
    }
    
    return 0;
}
