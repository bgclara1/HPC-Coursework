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


void variableInitialisation(int totalSteps, int numParticles,vector<vector<double>>& X,vector<vector<double>>& Y,
    vector<vector<double>>& Z,vector<vector<double>>& U,vector<vector<double>>& V,vector<vector<double>>& W,
    vector<vector<double>>& E,vector<vector<double>>& speed,vector<vector<double>>& xij,vector<vector<double>>& yij,
    vector<vector<double>>& zij,vector<vector<double>>& rij,vector<vector<double>>& dPhi_dx,vector<vector<double>>& dPhi_dy,
    vector<vector<double>>& dPhi_dz,vector<vector<double>>& Fx,vector<vector<double>>& Fy,vector<vector<double>>& Fz)
{
    X.resize(totalSteps, vector<double>(numParticles, 0.0));
    Y.resize(totalSteps, vector<double>(numParticles, 0.0));
    Z.resize(totalSteps, vector<double>(numParticles, 0.0));
    U.resize(totalSteps, vector<double>(numParticles, 0.0));
    V.resize(totalSteps, vector<double>(numParticles, 0.0));
    W.resize(totalSteps, vector<double>(numParticles, 0.0));
    E.resize(totalSteps, vector<double>(numParticles, 0.0));
    speed.resize(totalSteps, vector<double>(numParticles, 0.0));

    xij.resize(numParticles, vector<double>(numParticles, 0.0));
    yij.resize(numParticles, vector<double>(numParticles, 0.0));
    zij.resize(numParticles, vector<double>(numParticles, 0.0));
    rij.resize(numParticles, vector<double>(numParticles, 0.0));
    dPhi_dx.resize(numParticles, vector<double>(numParticles, 0.0));
    dPhi_dy.resize(numParticles, vector<double>(numParticles, 0.0));
    dPhi_dz.resize(numParticles, vector<double>(numParticles, 0.0));

    Fx.resize(totalSteps, vector<double>(numParticles, 0.0));
    Fy.resize(totalSteps, vector<double>(numParticles, 0.0));
    Fz.resize(totalSteps, vector<double>(numParticles, 0.0));
}

void icRandom(int numParticles, double Lx, double Ly, double Lz, double percent_type1,
    vector<vector<double>>& X,vector<vector<double>>& Y,vector<vector<double>>& Z,
    vector<vector<double>>& U,vector<vector<double>>& V,vector<vector<double>>& W,
    vector<double>& type)
{
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
                if (dx * dx + dy * dy + dz * dz < 0.25) {  // 0.5^2 = 0.25
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
    }
    int numType1 = (int)ceil(numParticles * (percent_type1 / 100.0));
    vector<double> particleTypes(numParticles, 0);
    for (int i = 0; i < numType1; i++) {
        particleTypes[i] = 1;
    }
    for (int i = 0; i < numParticles; i++) {
        int j = rand() % numParticles;
        int temp = particleTypes[i];
        particleTypes[i] = particleTypes[j];
        particleTypes[j] = temp;
    }
    type = particleTypes;
}

map<string, map<string, vector<double>>> getTestCases() {
    map<string, map<string, vector<double>>> testCaseDict;
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
    return testCaseDict;
}

void updateVars(const vector<double>& timestamps, int numParticles, double dt,double Lx, double Ly, double Lz,
                vector<double>& type,double temperature, bool tempProvided, double kb,const int epsilon[2][2], const int sigma[2][2],
                vector<vector<double>>& X,vector<vector<double>>& Y,vector<vector<double>>& Z,
                vector<vector<double>>& U,vector<vector<double>>& V,vector<vector<double>>& W,
                vector<vector<double>>& E,vector<vector<double>>& speed,vector<vector<double>>& xij,
                vector<vector<double>>& yij,vector<vector<double>>& zij,vector<vector<double>>& rij,
                vector<vector<double>>& dPhi_dx,vector<vector<double>>& dPhi_dy,
                vector<vector<double>>& dPhi_dz,vector<vector<double>>& Fx,
                vector<vector<double>>& Fy,vector<vector<double>>& Fz)
{
    int totalSteps = timestamps.size();
    for (int t = 0; t < totalSteps - 1; t++) {
        if (t % 5000 == 0) {
            cout << "Time: " << t * dt << endl;
        }
        for (int i = 0; i < numParticles; i++) {
            for (int j = i + 1; j < numParticles; j++) {
                xij[i][j] = X[t][i] - X[t][j];
                yij[i][j] = Y[t][i] - Y[t][j];
                zij[i][j] = Z[t][i] - Z[t][j];
                rij[i][j] = xij[i][j]*xij[i][j] + yij[i][j]*yij[i][j] + zij[i][j]*zij[i][j]; // r squared

                int t1 = type[i];
                int t2 = type[j];
                int e = epsilon[t1][t2];
                int s = sigma[t1][t2];

                double dPhi_coeff = -24 * e * ((2 * pow(s, 12) / pow(rij[i][j], 7)) - (pow(s, 6) / pow(rij[i][j], 4)));
                dPhi_dx[i][j] = xij[i][j] * dPhi_coeff;
                dPhi_dy[i][j] = yij[i][j] * dPhi_coeff;
                dPhi_dz[i][j] = zij[i][j] * dPhi_coeff;
            }
        }
        for (int i = 0; i < numParticles; i++) {
            for (int j = i + 1; j < numParticles; j++) {
                Fx[t][i] -= dPhi_dx[i][j];
                Fy[t][i] -= dPhi_dy[i][j];
                Fz[t][i] -= dPhi_dz[i][j];
                Fx[t][j] += dPhi_dx[i][j];
                Fy[t][j] += dPhi_dy[i][j];
                Fz[t][j] += dPhi_dz[i][j];
            }
        }
        for (int i = 0; i < numParticles; i++) {
            int m = (type[i] == 0) ? 1 : 10;
            U[t+1][i] = U[t][i] + dt * Fx[t][i] / m;
            V[t+1][i] = V[t][i] + dt * Fy[t][i] / m;
            W[t+1][i] = W[t][i] + dt * Fz[t][i] / m;
        }
        for (int i = 0; i < numParticles; i++) {
            int m = (type[i] == 0) ? 1 : 10;
            speed[t+1][i] = sqrt(U[t+1][i]*U[t+1][i] + V[t+1][i]*V[t+1][i] + W[t+1][i]*W[t+1][i]);
            E[t+1][i] = 0.5 * m * speed[t+1][i] * speed[t+1][i];
            if (tempProvided) {
                double BoltzTemp = (2.0 / (3.0 * kb)) * E[t+1][i];
                double lambda = sqrt(temperature / BoltzTemp);
                U[t+1][i] *= lambda;
                V[t+1][i] *= lambda;
                W[t+1][i] *= lambda;
            }
        }
        for (int i = 0; i < numParticles; i++) {
            X[t+1][i] = X[t][i] + dt * U[t][i];
            Y[t+1][i] = Y[t][i] + dt * V[t][i];
            Z[t+1][i] = Z[t][i] + dt * W[t][i];
            // Apply BCs
            if (X[t+1][i] > Lx) {
                X[t+1][i] = 2*Lx - X[t+1][i];
                U[t+1][i] = -abs(U[t+1][i]);
            }
            if (Y[t+1][i] > Ly) {
                Y[t+1][i] = 2*Ly - Y[t+1][i];
                V[t+1][i] = -abs(V[t+1][i]);
            }
            if (Z[t+1][i] > Lz) {
                Z[t+1][i] = 2*Lz - Z[t+1][i];
                W[t+1][i] = -abs(W[t+1][i]);
            }
            if (X[t+1][i] < 0) {
                X[t+1][i] = -X[t+1][i];
                U[t+1][i] = abs(U[t+1][i]);
            }
            if (Y[t+1][i] < 0) {
                Y[t+1][i] = -Y[t+1][i];
                V[t+1][i] = abs(V[t+1][i]);
            }
            if (Z[t+1][i] < 0) {
                Z[t+1][i] = -Z[t+1][i];
                W[t+1][i] = abs(W[t+1][i]);
            }
        }
    }
}

void writeToFiles(int totalSteps, int numParticles, const vector<double>& timestamps,const vector<vector<double>>& X, const vector<vector<double>>& Y,
    const vector<vector<double>>& Z, const vector<vector<double>>& U,const vector<vector<double>>& V, const vector<vector<double>>& W,
    const vector<vector<double>>& E)
    {
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
        << " E = " << E[t][i] << "\n";
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
    energyfile << " " << E[t][i];
    }
    energyfile << "\n";
    }
    energyfile.close();

    ofstream posfile("positions.txt");
    posfile << "runtime";
    for (int i = 0; i < numParticles; i++) {
    posfile << " x" << i << " y" << i;
    }
    posfile << "\n";
    for (int t = 0; t < totalSteps; t++) {
    posfile << defaultfloat << timestamps[t];
    for (int i = 0; i < numParticles; i++) {
    posfile << " " << fixed << setprecision(6) << X[t][i]
        << " " << fixed << setprecision(6) << Y[t][i];
    }
    posfile << "\n";
    }
    posfile.close();
    }



/////////////////////////////////////////////////////
////////////////////////////////////////////////////
//
//          M A I N   P R O G R A M
//
////////////////////////////////////////////////////
////////////////////////////////////////////////////


    
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
    bool tempProvided = false;

    vector<vector<double>> X, Y, Z, U, V, W, E, speed,xij, yij, zij, rij,dPhi_dx, dPhi_dy, dPhi_dz,Fx, Fy, Fz;

    map<string, map<string, vector<double>>> testCaseDict = getTestCases();
    
    double runtime, percent_type1, temperature;
    double kb = 0.8314459920816467;
    int numParticles;
    vector<double> x, y, z, u, v, w, type;
    while (i < argc) {
      //  cout << "Argument " << i + 1 << ": " << argv[i] << endl;
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
        
        } else if (string(argv[i]) == "--temp")  {
            temperature = stod(argv[i+1]);
            tempProvided = true;
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
            cout << "Allowed options:\n"
                 << "--help                Print available options.\n"
                 << "--Lx arg (=20)        x length (Angstroms)\n"
                 << "--Ly arg (=20)        y length (Angstroms)\n"
                 << "--Lz arg (=20)        z length (Angstroms)\n"
                 << "--dt arg (=0.001)     Time-step\n"
                 << "--T arg               Final time\n"
                 << "--ic-one              Initial condition: one stationary particle\n"
                 << "--ic-one-vel          Initial condition: one moving particle\n"
                 << "--ic-two              Initial condition: two bouncing particles\n"
                 << "--ic-two-pass1        Initial condition: two passing particles close\n"
                 << "--ic-two-pass2        Initial condition: two passing particles close\n"
                 << "--ic-two-pass3        Initial condition: two passing particles close, heavy\n"
                 << "--percent-type1 arg (=10)  Percentage of type 1 particles with random IC\n"
                 << "--ic-random           Number of particles to spawn with random IC\n"
                 << "--temp arg            Temperature (degrees Kelvin)\n";
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

    
    variableInitialisation(totalSteps, numParticles,X, Y, Z,U, V, W,E, speed, xij, yij, zij, rij,
        dPhi_dx, dPhi_dy, dPhi_dz, Fx, Fy, Fz);

    if (icRandomChosen) {
        icRandom(numParticles, Lx, Ly, Lz, percent_type1,
                                          X, Y, Z, U, V, W, type);
    } else {
        for (int i = 0; i < numParticles; i++) {
            X[0][i] = x[i];             // the lower case x,y,z etc are in testCaseDict
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

    updateVars(timestamps, numParticles, dt, Lx, Ly, Lz,type, temperature, tempProvided, kb,
        epsilon, sigma,X, Y, Z,U, V, W,E, speed,xij, yij, zij, rij,dPhi_dx, dPhi_dy, dPhi_dz,Fx, Fy, Fz);            
    

    writeToFiles(totalSteps, numParticles, timestamps, X, Y, Z, U, V, W, E);



    return 0;
}










