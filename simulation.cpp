#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream> 
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
    testCaseDict["--Claras-test-case"] = {
        {"time", {50.0}},
        {"numParticles", {3}},
        {"x", {8.5, 11.5, 10}},
        {"y", {11.3, 8.7,20}},
        {"z", {10.0, 10.0,13}},
        {"u", {0.5, -0.5,-0.5}},
        {"v", {0.0, 0.0,0}},
        {"w", {0.0, 0.0,0}},
        {"type", {1, 1,0}}
    };
    
    double time;
    int numParticles;
    vector<double> x, y, z, u, v, w, type;
    
    // test reading in args properly
    while (i < argc) {  //argc num args provided
        cout << "Argument " << i + 1 << ": " << argv[i] << endl;
    
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
            exit(1);
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
    
    // GENERATE TYPE IF NOT --IC-...
    
    ////////////////// SIM ALGO /////////////////////////
    
    double steps = (time / dt) + 1;
    vector<double> timestamps;
    for (int i = 0; i < steps; i++) {
        timestamps.push_back(i * dt);
    }
    
    int epsilon[2][2] = {
        {3,15},
        {15,60}
    };
    
    int sigma[2][2] = {
        {1,2},
        {2,3}
    };
    
    vector<vector<double>> xij(numParticles, vector<double>(numParticles, 0.0)); // initialise w zeros
    vector<vector<double>> yij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> zij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> rij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dx(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dy(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dz(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> F(numParticles, vector<double>(numParticles, 0.0));
    
    for (int i = 0; i < numParticles - 1; i++) {            // only calculating upper triangle (no diag ) avoids double counting
        for (int j = i + 1; j < numParticles; j++) {
            xij[i][j] = x[i] - x[j];
            yij[i][j] = y[i] - y[j];
            zij[i][j] = z[i] - z[j];
            rij[i][j] = sqrt(xij[i][j]*xij[i][j] + yij[i][j]*yij[i][j] + zij[i][j]*zij[i][j]);
            int t1 = type[i];
            int t2 = type[j];
            int e = epsilon[t1][t2];
            int s = sigma[t1][t2];
            dPhi_dx[i][j] = -24 * e * xij[i][j]*((2*pow(s,12)/pow(rij[i][j],14)) - (pow(s,6)/pow(rij[i][j],8)));
            dPhi_dy[i][j] = -24 * e * yij[i][j]*((2*pow(s,12)/pow(rij[i][j],14)) - (pow(s,6)/pow(rij[i][j],8)));
            dPhi_dz[i][j] = -24 * e * zij[i][j]*((2*pow(s,12)/pow(rij[i][j],14)) - (pow(s,6)/pow(rij[i][j],8)));
        }
    }
    
    cout << "Force matrix F:" << endl;
    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < numParticles; j++) {
            cout << F[i][j] << "\t";
        }
        cout << endl;
    }
    
    vector<double> Fx(numParticles, 0.0), Fy(numParticles, 0.0), Fz(numParticles, 0.0);
    for (int i = 0; i < numParticles - 1; i++) {
        for (int j = i + 1; j < numParticles; j++) {
            Fx[i] += dPhi_dx[i][j];
            Fy[i] += dPhi_dy[i][j];
            Fz[i] += dPhi_dz[i][j];
            Fx[j] -= dPhi_dx[i][j];
            Fy[j] -= dPhi_dy[i][j];
            Fz[j] -= dPhi_dz[i][j];
        }
    }
    
    int m;
    for (int i = 0; i < numParticles - 1; i++) {
        if (type[i] == 0) {
            m = 1;
        } else {
            m = 10;
        }
        u[i+1] = u[i] + dt * Fx[i]/m;
        v[i+1] = v[i] + dt * Fy[i]/m;
        w[i+1] = w[i] + dt * Fz[i]/m;
        x[i+1] = x[i] + dt * u[i];
        y[i+1] = y[i] + dt * v[i];
        z[i+1] = z[i] + dt * w[i];
    }
    
    // Open file once outside the simulation loop
    ofstream outfile("output.txt");
    for (int t = 0; t < time; t++) {
        outfile << "Time step " << t << endl;
        outfile << "Net force arrays:" << endl;
        outfile << "Fx: ";
        for (int i = 0; i < numParticles; i++) {
            outfile << Fx[i] << "\t";
        }
        outfile << "\nFy: ";
        for (int i = 0; i < numParticles; i++) {
            outfile << Fy[i] << "\t";
        }
        outfile << "\nFz: ";
        for (int i = 0; i < numParticles; i++) {
            outfile << Fz[i] << "\t";
        }
        outfile << "\nUpdated positions and velocities:" << endl;
        for (int i = 0; i < numParticles; i++) {
            outfile << "Particle " << i << ": "
                    << " x = " << x[i]
                    << ", y = " << y[i]
                    << ", z = " << z[i]
                    << " ; u = " << u[i]
                    << ", v = " << v[i]
                    << ", w = " << w[i]
                    << endl;
        }
        outfile << "\n";
    }
    outfile.close();
    
    return 0;
}
