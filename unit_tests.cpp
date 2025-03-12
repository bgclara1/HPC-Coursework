#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

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

int main() {
    cout << "Compiling simulation.cpp..." << endl;

    vector<string> testCases = {
        "./simulation --ic-one",
        "./simulation --ic-one-vel",
        "./simulation --ic-two",
        "./simulation --ic-two-pass1",
        "./simulation --ic-two-pass2",
        "./simulation --ic-two-pass3"
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
        
        cout << "Final Positions: " << endl;
        cout << "Time Stamp, X position, Y position ... " << endl;
        cout << finalPositions << endl;

        cout << "Final Kinetic Energy: " << endl;
        cout << "Time Stamp, Kinetic Energy ... " << endl;
        cout << finalEnergies << endl;
    }
    
    return 0;
}



