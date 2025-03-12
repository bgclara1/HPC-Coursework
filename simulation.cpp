#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

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
        {"runtime", {0.2}},
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
        {"runtime", {20.0}},
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
        {"runtime", {50}},
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
        {"runtime", {50.0}},
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
        {"runtime", {50.0}},
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
        {"runtime", {50.0}},
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
        {"runtime", {5}},
        {"numParticles", {3}},
        {"x", {8.5, 11.5, 10}},
        {"y", {11.3, 8.7, 20}},
        {"z", {10.0, 10.0, 13}},
        {"u", {0.5, -0.5, -0.5}},
        {"v", {0.0, 0.0, 0}},
        {"w", {0.0, 0.0, 0}},
        {"type", {1, 1, 0}}
    };
    
    double runtime, percent_type1;
    int numParticles;
    vector<double> x, y, z, u, v, w, type;
    while (i < argc) {
        cout << "Argument " << i + 1 << ": " << argv[i] << endl;
        if (string(argv[i]) == "--Lx") {
            Lx = stod(argv[i+1]);
        } else if (string(argv[i]) == "--Ly") {
            Ly = stod(argv[i+1]);
        } else if (string(argv[i]) == "--Lz") {
            Lz = stod(argv[i+1]);
        } else if (string(argv[i]) == "--T") {
            runtime = stod(argv[i+1]);
            timeProvided = true;
        } else if (string(argv[i]) == "--N") {
            numParticles = stod(argv[i+1]);
            nProvided = true;

        } else if (string(argv[i]) == "--percent-type1") {
            percent_type1 = stod(argv[i+1]);
        } else if (string(argv[i]) == "--ic-random") {
            icRandomChosen = true;
        } else if (testCaseDict.find(string(argv[i])) != testCaseDict.end()) {
            string key(argv[i]);
            runtime = testCaseDict[key]["runtime"][0];
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
            cout << "--help" << endl;
            exit(1);
        }
        i++;
    }
    
    if ((testCase == true) || (icRandomChosen == true && nProvided == true && timeProvided == true)) {
        cout << "Command Line input well-formatted, carrying on..." << endl;
    } else {
        cout << "Command line input formatted incorrectly, exiting program." << endl;
        exit(1);
    }
    
    double totalSteps = (runtime / dt) + 1;
    vector<double> timestamps(totalSteps);
    for (int i = 0; i < totalSteps; i++) {
        timestamps[i] = i * dt;
    }

    
    
    vector<vector<double>> X(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> Y(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> Z(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> U(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> V(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> W(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> E_history(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> speed_history(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> xij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> yij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> zij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> rij(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dx(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dy(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> dPhi_dz(numParticles, vector<double>(numParticles, 0.0));
    vector<vector<double>> Fx(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> Fy(totalSteps, vector<double>(numParticles, 0.0));
    vector<vector<double>> Fz(totalSteps, vector<double>(numParticles, 0.0));

    

    if (icRandomChosen) {
        srand(time(0));  
        for (int i = 0; i < numParticles; i++) {
            double cx, cy, cz;
            while (true) {
                cx = ((double)rand() / RAND_MAX) * Lx;
                cy = ((double)rand() / RAND_MAX) * Ly;
                cz = ((double)rand() / RAND_MAX) * Lz;
                bool valid = true;
                for (int j = 0; j < i; j++) {
                    double dx = cx - X[0][j];
                    double dy = cy - Y[0][j];
                    double dz = cz - Z[0][j];
                    if (dx*dx + dy*dy + dz*dz < 0.25) {  // 0.5^2 = 0.25
                        valid = false;
                        break;
                    }
                }
                if (valid)
                    break;
            }
            X[0][i] = cx;
            Y[0][i] = cy;
            Z[0][i] = cz;
            U[0][i] = ((double)rand() / RAND_MAX) - 0.5;
            V[0][i] = ((double)rand() / RAND_MAX) - 0.5;
            W[0][i] = ((double)rand() / RAND_MAX) - 0.5;

            int numType1 = (int)ceil(numParticles * (percent_type1 / 100.0));
            vector<int> particleTypes(numParticles, 0);
            for (int i = 0; i < numType1; i++) {
                particleTypes[i] = 1;
            }
            for (int i = 0; i < numParticles; i++) {
                int j = rand() % numParticles;
                int temp = particleTypes[i];
                particleTypes[i] = particleTypes[j];
                particleTypes[j] = temp;
            }
            type.resize(numParticles);
            for (int i = 0; i < numParticles; i++) {
                type[i] = particleTypes[i];
            }


        }

    } else {
        for (int i = 0; i < numParticles; i++) {
            // X at time 0 at position i = 0, 1 ...
            X[0][i] = x[i];
            Y[0][i] = y[i];
            Z[0][i] = z[i];
            U[0][i] = u[i];
            V[0][i] = v[i];
            W[0][i] = w[i];
        }
    }

    int epsilon[2][2] = { {3,15}, {15,60} };
    int sigma[2][2] = { {1,2}, {2,3} };

    int m;
    for (int t = 0; t < timestamps.size()-1; t++) {
        for (int i = 0; i < numParticles ; i++) {
            for (int j = i + 1; j < numParticles; j++) {
                xij[i][j] = X[t][i] - X[t][j];
                yij[i][j] = Y[t][i] - Y[t][j];
                zij[i][j] = Z[t][i] - Z[t][j];
                rij[i][j] = (xij[i][j]*xij[i][j] + yij[i][j]*yij[i][j] + zij[i][j]*zij[i][j]);   // r squared

                int t1 = type[i];
                int t2 = type[j];
                int e = epsilon[t1][t2];
                int s = sigma[t1][t2];

                double dPhi_coeff = -24*e*(  (2*pow(s,12)/pow(rij[i][j],7))   -   (pow(s,6)/pow(rij[i][j],4))   );

                dPhi_dx[i][j] = xij[i][j] * dPhi_coeff;
                dPhi_dy[i][j] = yij[i][j] * dPhi_coeff;
                dPhi_dz[i][j] = zij[i][j] * dPhi_coeff;
            }
            
        }
        for (int i = 0; i < numParticles ; i++) {
            for (int j = i + 1; j < numParticles; j++) {
                if (i!=j) {
                    Fx[t][i] -= dPhi_dx[i][j] ;
                    Fy[t][i] -= dPhi_dy[i][j];
                    Fz[t][i] -= dPhi_dz[i][j];
                    Fx[t][j] += dPhi_dx[i][j] ;
                    Fy[t][j] += dPhi_dy[i][j];
                    Fz[t][j] += dPhi_dz[i][j];
                }
            }
        }        

        

        for (int i = 0; i < numParticles; i++) {
            if (type[i] == 0) {
                m = 1;
            } else {
                m = 10;
            }

            U[t+1][i] = U[t][i] + dt * Fx[t][i] / m;
            V[t+1][i] = V[t][i] + dt * Fy[t][i] / m;
            W[t+1][i] = W[t][i] + dt * Fz[t][i] / m;

            X[t+1][i] = X[t][i] + dt * U[t][i];
            Y[t+1][i] = Y[t][i] + dt * V[t][i];
            Z[t+1][i] = Z[t][i] + dt * W[t][i];



            if (X[t+1][i] > Lx) {
                X[t+1][i] = 2*Lx - X[t+1][i];
                U[t+1][i] = -abs(U[t+1][i]);
                cout << "boundary" << endl;
            }
            if (Y[t+1][i] > Ly) {
                Y[t+1][i] = 2*Ly - Y[t+1][i];
                V[t+1][i] = -abs(V[t+1][i]);
                cout << "boundary" << endl;
            }
            if (Z[t+1][i] > Lz) {
                Z[t+1][i] = 2*Lz - Z[t+1][i];
                W[t+1][i] = -abs(W[t+1][i]);
                cout << "boundary" << endl;
            }

            if (X[t+1][i] < 0) {
                X[t+1][i] = - X[t+1][i];
                U[t+1][i] = abs(U[t+1][i]);
                cout << "boundary" << endl;
            }
            if (Y[t+1][i] < 0) {
                Y[t+1][i] = - Y[t+1][i];
                V[t+1][i] = abs(V[t+1][i]);
                cout << "boundary" << endl;
            }
            if (Z[t+1][i] < 0) {
                Z[t+1][i] = - Z[t+1][i];
                W[t+1][i] = abs(W[t+1][i]);
                cout << "boundary" << endl;
            }

            }
            for (int i = 0; i < numParticles; i++) {
                if (type[i] == 0) {
                    m = 1;
                } else {
                    m = 10;
                }
                speed_history[t][i] = sqrt(U[t][i]*U[t][i] + V[t][i]*V[t][i] + W[t][i]*W[t][i]);
                E_history[t][i] = 0.5 * m * speed_history[t][i] * speed_history[t][i];
            }
    }
    
    ofstream outfile("output.txt");
    for (int t = 0; t < totalSteps; t++) {
        outfile << "Time step " << t << "\n";
        for (int i = 0; i < numParticles; i++) {
            outfile << "Particle " << i << ": x = " << X[t][i] 
                    << " y = " << Y[t][i] 
                    << " z = " << Z[t][i]
                    << " u = " << U[t][i] 
                    << " v = " << V[t][i] 
                    << " w = " << W[t][i] 
                    << " E = " << E_history[t][i] << "\n";
        }
        outfile << "\n";
    }
    outfile.close();
    
    ofstream energyfile("energy.txt");
    energyfile << "runtime";
    for (int i = 0; i < numParticles; i++) {
        energyfile << " E" << i;
    }
    energyfile << "\n";
    for (int t = 0; t < totalSteps; t++) {
        energyfile << timestamps[t];
        for (int i = 0; i < numParticles; i++) {
            energyfile << " " << E_history[t][i];
        }
        energyfile << "\n";
    }
    energyfile.close();

    ofstream posfile("positions.txt");
     posfile << fixed << setprecision(10);
    posfile << "runtime";
    for (int i = 0; i < numParticles; i++) {
        posfile << " x" << i << " y" << i;
    }
    posfile << "\n";
    for (int i = 0; i < totalSteps; i++) {
        posfile << timestamps[i];
        for (int j = 0; j < numParticles; j++) {
           

            posfile << " " << X[i][j] << " " << Y[i][j];
        }
        posfile << "\n";
    }
    posfile.close();


    return 0;
}










