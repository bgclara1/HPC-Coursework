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
    string line, lastLine;
    while (getline(infile, line)) {
        if (!line.empty())
            lastLine = line;
    }
    return lastLine;
}

void printPositions(const string &line) {
    
    stringstream ss(line);// split  line into elements.
    vector<string> elements;
    string element;
    while(ss >> element) {
        elements.push_back(element);
    }
    
    cout << "Final Positions:" << endl;
    cout << "  Timestamp: " << elements[0] << endl;
    
    if (elements.size() == 3) {   // as in time , x , y so one particle case
        cout << "  Particle 0: X = " << elements[1] << ", Y = " << elements[2] << endl;
    } else if (elements.size() == 5) { // as in time x1 y1 x2 y2 , two partivcle case
        cout << "  Particle 0: X = " << elements[1] << ", Y = " << elements[2] << endl;
        cout << "  Particle 1: X = " << elements[3] << ", Y = " << elements[4] << endl;
    } else {
        cout << "  ";
        for (const auto &i : elements) {
            cout << i << " ";
        }
        cout << endl;
    }
}

void printEnergies(const string &line) {
    stringstream ss(line);
    vector<string> elements;
    string element;
    while(ss >> element) {
        elements.push_back(element);
    }
    
    cout << "Final Kinetic Energies:" << endl;
    cout << "  Timestamp: " << elements[0] << endl;
    
    if (elements.size() == 2) {
        cout << "  Particle 0: KE = " << elements[1] << endl;
    } else if (elements.size() == 3) {
        cout << "  Particle 0: KE = " << elements[1] << endl;
        cout << "  Particle 1: KE = " << elements[2] << endl;
    } else {
        cout << "  ";
        for (const auto &i : elements) {
            cout << i << " ";
        }
        cout << endl;
    }
}

int main() {

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
