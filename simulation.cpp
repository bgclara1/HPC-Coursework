#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// cmd + / to comment out on mass
// g++ simulation.cpp -std=c++11 -o simulation



int main(int argc, char *argv[]) {              //read cmd args w main params.
    int i = 0;
    double Lx = 20;
    double Ly = 20;
    double Lz = 20;
    double dt = 0.001;
    bool testCase = false;
    bool timeProvided = false;
    bool nProvided = false;
    bool icRandomChosen = false;
    

    map<string, map<string, vector<double> > > testCaseDict;
    testCaseDict["--ic-one"] = {
        {"time", {1.0}},
        {"numParticles", {1}},
        {"x", {10.0}},
        {"y", {10.0}},
        {"z", {10.0}},
        {"u", {0.0}},
        {"v", {0.0}},
        {"w", {0.0}},
        {"type", {0}}
    };
    testCaseDict["--ic-one-vel"] = {
        {"time", {20.0}},
        {"numParticles", {1}},
        {"x", {10.0}},
        {"y", {10.0}},
        {"z", {10.0}},
        {"u", {5.0}},
        {"v", {2.0}},
        {"w", {1.0}},
        {"type", {0}}
    };
    testCaseDict["--ic-two"] = {
        {"time", {50.0}},
        {"numParticles", {2}},
        {"x", {8.5, 11.5}},
        {"y", {10.0, 10.0}},
        {"z", {10.0, 10.0}},
        {"u", {0.0, 0.0}},
        {"v", {0.0, 0.0}},
        {"w", {0.0, 0.0}},
        {"type", {0, 0}}
    };
    testCaseDict["--ic-two-pass1"] = {
        {"time", {50.0}},
        {"numParticles", {2}},
        {"x", {8.5, 11.5}},
        {"y", {11.5, 8.5}},
        {"z", {10.0, 10.0}},
        {"u", {0.5, -0.5}},
        {"v", {0.0, 0.0}},
        {"w", {0.0, 0.0}},
        {"type", {0, 0}}
    };
    testCaseDict["--ic-two-pass2"] = {
        {"time", {50.0}},
        {"numParticles", {2}},
        {"x", {8.5, 11.5}},
        {"y", {11.3, 8.7}},
        {"z", {10.0, 10.0}},
        {"u", {0.5, -0.5}},
        {"v", {0.0, 0.0}},
        {"w", {0.0, 0.0}},
        {"type", {0, 0}}
    };
    testCaseDict["--ic-two-pass3"] = {
        {"time", {50.0}},
        {"numParticles", {2}},
        {"x", {8.5, 11.5}},
        {"y", {11.3, 8.7}},
        {"z", {10.0, 10.0}},
        {"u", {0.5, -0.5}},
        {"v", {0.0, 0.0}},
        {"w", {0.0, 0.0}},
        {"type", {1, 1}}
    };
    

    double time;
    int numParticles;
    vector<double> x, y, z, u, v, w, type;

    // test reading in args properly
    while (i < argc) {  //argc num args provided
        cout << "Argument " << i + 1 << ": " << argv[i]
             << endl;

        if (string(argv[i]) == "--Lx") {
            cout << "hey " << endl;
            Lx = stod(argv[i+1]);
            cout << argv[i+1] << endl;
        } else if (string(argv[i]) == "--Ly") {
            Ly = stod(argv[i+1]);
        } else if (string(argv[i]) == "--Lz") {
            Lz = stod(argv[i+1]);
        } else if (string(argv[i]) == "--T") {
            time = stod(argv[i+1]);
            timeProvided = true;
        } else if (string(argv[i]) == "--N") {
            numParticles = stod(argv[i+1]);
            nProvided = true;
        } else if (string(argv[i]) == "--ic-random") {
            icRandomChosen = true;
        } else if (testCaseDict.find(string(argv[i])) != testCaseDict.end()) {
            string key(argv[i]);
            time = testCaseDict[key]["time"][0];
            numParticles = testCaseDict[key]["numParticles"][0];
            x = testCaseDict[key]["x"];
            y = testCaseDict[key]["y"];
            z = testCaseDict[key]["z"];
            u = testCaseDict[key]["u"];
            v = testCaseDict[key]["v"];
            w = testCaseDict[key]["w"];
            type = testCaseDict[key]["type"];
            testCase = true;
        } else if (string(argv[i]) == "--help") {
            cout << "--help" << "                            " << "Print available options" << endl;
            cout << "--Lx arg (=20)" << "                    " << "x length (Angstroms)" << endl;
            cout << "--Ly arg (=20)" << "                    " << "y length (Angstroms)" << endl;
            cout << "--Lz arg (=20)" << "                    " << "z length (Angstroms)" << endl;
            cout << "--dt arg (=0.001)" << "                 " << "Time-step" << endl;
            cout << "--T" << "                               " << "Final Time" << endl;
            cout << "--ic-one" << "                          " << "Initial condition: one stationary particle" << endl;
            cout << "--ic-one-vel" << "                      " << "Initial condition: one moving particle" << endl;
            cout << "--ic-two" << "                          " << "Initial condition: two bouncing particles" << endl;
            cout << "--ic-two-pass1" << "                    " << "Initial condition: two passing particles" << endl;
            cout << "--ic-two-pass2" << "                    " << "Initial condition: two passing particles close" << endl;
            cout << "--ic-two-pass3" << "                    " << "Initial condition: two passing particles close, heavy" << endl;
            cout << "--ic-random" << "                       " << "Initial condition: N random particles" << endl;
            cout << "--percent-type1 arg (=10)" << "         " << "Percentage of type 1 particles with random IC" << endl;
            cout << "--N arg" << "                           " << "Number of particles to spawn with random IC" << endl;
            cout << "--temp arg" << "                        " << "Temperature (degrees Kelvin)" << endl;

        }

        i++;
    }

    if ((testCase == true) || (icRandomChosen == true && nProvided == true && timeProvided == true)) {
        cout << "Command Line input well-formatted, carrying on... " << endl;
    } else {
        cout << testCase << endl;
        cout << icRandomChosen << endl;
        cout << nProvided << endl;
        cout << timeProvided << endl;
        cout << "Command line input formatted incorrectly, exiting program." << endl;
        exit(1);
    }

    

    cout << time << endl;

    if (testCase == true) {
        cout << numParticles << endl;
        cout << x[0] << endl;
        cout << y[0] << endl;
        cout << z[0] << endl;
        cout << u[0] << endl;
        cout << v[0] << endl;
        cout << w[0] << endl;
        cout << type[0] << endl;
    }
    cout << Lx << endl;
    cout << Ly << endl;
    cout << Lz << endl;





    return 0;

}