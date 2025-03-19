#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace std;

// Global constant lookup table for s^6 (index 0 unused)
const double s6_table[4] = {0.0, 1.0, 64.0, 729.0};

void variableInitialisation(int totalSteps, int numParticles,
    vector<double>& X, vector<double>& Y, vector<double>& Z,
    vector<double>& U, vector<double>& V, vector<double>& W,
    vector<double>& E, vector<double>& speed,
    vector<double>& Fx, vector<double>& Fy, vector<double>& Fz)
{
    int size = totalSteps * numParticles;
    X.resize(size, 0.0);
    Y.resize(size, 0.0);
    Z.resize(size, 0.0);
    U.resize(size, 0.0);
    V.resize(size, 0.0);
    W.resize(size, 0.0);
    E.resize(size, 0.0);
    speed.resize(size, 0.0);
    Fx.resize(size, 0.0);
    Fy.resize(size, 0.0);
    Fz.resize(size, 0.0);
}

void icRandom(int numParticles, double Lx, double Ly, double Lz, double percent_type1,
    vector<double>& X, vector<double>& Y, vector<double>& Z,
    vector<double>& U, vector<double>& V, vector<double>& W,
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
                double dx = cx - X[j];
                double dy = cy - Y[j];
                double dz = cz - Z[j];
                if (dx * dx + dy * dy + dz * dz < 0.25) { // 0.5^2 = 0.25
                    valid = false;
                    break;
                }
            }
            if (valid)
                break;
        }
        X[i] = cx;
        Y[i] = cy;
        Z[i] = cz;
        U[i] = ((double)rand() / RAND_MAX) - 0.5;
        V[i] = ((double)rand() / RAND_MAX) - 0.5;
        W[i] = ((double)rand() / RAND_MAX) - 0.5;
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
        {"runtime", {1}},
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

//-----------------------------------------------------------------
// Serial flattened force updater
//-----------------------------------------------------------------
// This function uses a single loop over the total number of unique pairs.
void updateVarsSerialFlattened(int numParticles, double dt, double Lx, double Ly, double Lz,
    const vector<double>& type, bool tempProvided, double temperature, double kb,
    const int epsilon[2][2], const int sigma[2][2],
    vector<double>& X, vector<double>& Y, vector<double>& Z,
    vector<double>& U, vector<double>& V, vector<double>& W,
    vector<double>& E, vector<double>& speed,
    vector<double>& Fx, vector<double>& Fy, vector<double>& Fz)
{
    // Flattened force loop
    long long totalPairs = (long long)numParticles * (numParticles - 1) / 2;
    for (long long p = 0; p < totalPairs; p++) {
        double tmp = sqrt(8.0 * p + 1.0);
        int i = (int)((tmp - 1.0) / 2.0);
        int j = p - (i * (i + 1) / 2) + i + 1;

        double xij = X[i] - X[j];
        double yij = Y[i] - Y[j];
        double zij = Z[i] - Z[j];
        double rij = xij*xij + yij*yij + zij*zij;  // squared distance

        int t1 = static_cast<int>(type[i]);
        int t2 = static_cast<int>(type[j]);
        int e = epsilon[t1][t2];
        int s = sigma[t1][t2];

        double inv_r4 = 1.0 / (rij * rij * rij * rij);
        double sigma6_val = s6_table[s] * inv_r4;
        double sigma12_val = sigma6_val * sigma6_val * rij;
        double coeff = -24.0 * e * (2.0 * sigma12_val - sigma6_val);

        Fx[i] -= xij * coeff;
        Fy[i] -= yij * coeff;
        Fz[i] -= zij * coeff;
        Fx[j] += xij * coeff;
        Fy[j] += yij * coeff;
        Fz[j] += zij * coeff;
    }

    // Update velocities.
    for (int i = 0; i < numParticles; i++) {
        int m = (type[i] == 0) ? 1 : 10;
        U[i] += dt * Fx[i] / m;
        V[i] += dt * Fy[i] / m;
        W[i] += dt * Fz[i] / m;
    }

    // Compute energies.
    double E_total = 0.0;
    for (int i = 0; i < numParticles; i++) {
        int m = (type[i] == 0) ? 1 : 10;
        double speed2 = U[i] * U[i] + V[i] * V[i] + W[i] * W[i];
        E[i] = 0.5 * m * speed2;
        E_total += E[i];
    }

    // Temperature scaling.
    if (tempProvided) {
        double currentTemp = (2.0 / (3.0 * numParticles * kb)) * E_total;
        double lambda = sqrt(temperature / currentTemp);
        for (int i = 0; i < numParticles; i++) {
            U[i] *= lambda;
            V[i] *= lambda;
            W[i] *= lambda;
        }
    }

    // Update positions and apply boundary conditions.
    for (int i = 0; i < numParticles; i++) {
        X[i] += dt * U[i];
        Y[i] += dt * V[i];
        Z[i] += dt * W[i];
        if (X[i] > Lx) {
            X[i] = 2 * Lx - X[i];
            U[i] = -abs(U[i]);
        }
        if (Y[i] > Ly) {
            Y[i] = 2 * Ly - Y[i];
            V[i] = -abs(V[i]);
        }
        if (Z[i] > Lz) {
            Z[i] = 2 * Lz - Z[i];
            W[i] = -abs(W[i]);
        }
        if (X[i] < 0) {
            X[i] = -X[i];
            U[i] = abs(U[i]);
        }
        if (Y[i] < 0) {
            Y[i] = -Y[i];
            V[i] = abs(V[i]);
        }
        if (Z[i] < 0) {
            Z[i] = -Z[i];
            W[i] = abs(W[i]);
        }
    }
}

void writeToFiles(int t, int numParticles, const vector<double>& timestamps,
    const vector<double>& X, const vector<double>& Y,
    const vector<double>& Z, const vector<double>& U,
    const vector<double>& V, const vector<double>& W,
    const vector<double>& E)
{
    ofstream outfile("output.txt", ios::app);
    outfile << "Time step " << t << "\n";
    for (int i = 0; i < numParticles; i++) {
        outfile << "Particle " << i << ": x = " << X[i]
                << " y = " << Y[i]
                << " z = " << Z[i]
                << " u = " << U[i]
                << " v = " << V[i]
                << " w = " << W[i]
                << " E = " << E[i] << "\n";
    }
    outfile << "\n";
    outfile.close();

    // Append energy data.
    ofstream energyfile("energy.txt", ios::app);
    energyfile << "runtime";
    for (int i = 0; i < numParticles; i++) {
        energyfile << " E" << i;
    }
    energyfile << "\n";
    energyfile << timestamps[t];
    for (int i = 0; i < numParticles; i++) {
        energyfile << " " << E[i];
    }
    energyfile << "\n";
    energyfile.close();

    // Append positions.
    ofstream posfile("positions.txt", ios::app);
    posfile << "runtime";
    for (int i = 0; i < numParticles; i++) {
        posfile << " x" << i << " y" << i;
    }
    posfile << "\n";
    posfile << defaultfloat << timestamps[t];
    for (int i = 0; i < numParticles; i++) {
        posfile << " " << fixed << setprecision(6) << X[i]
                << " " << fixed << setprecision(6) << Y[i];
    }
    posfile << "\n";
    posfile.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//                              MAIN PROGRAM
///////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    auto start = chrono::high_resolution_clock::now();
    int i = 0;
    double Lx = 20, Ly = 20, Lz = 20;
    double dt = 0.001;
    bool testCase = false, timeProvided = false, nProvided = false, icRandomChosen = false, tempProvided = false;

    // Clean up old output files.
    ifstream file1("output.txt");
    if (file1) { file1.close(); remove("output.txt"); }
    ifstream file2("energy.txt");
    if (file2) { file2.close(); remove("energy.txt"); }
    ifstream file3("positions.txt");
    if (file3) { file3.close(); remove("positions.txt"); }

    vector<double> X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz;
    double xij, yij, zij, rij; // temporary variables used in force updater
    map<string, map<string, vector<double>>> testCaseDict = getTestCases();
    double runtime, percent_type1, temperature;
    double kb = 0.8314459920816467;
    int numParticles;
    vector<double> x, y, z, u, v, w, type;

    while (i < argc) {
        string arg = argv[i];
        if (arg == "--Lx") {
            Lx = stod(argv[i+1]);
        } else if (arg == "--Ly") {
            Ly = stod(argv[i+1]);
        } else if (arg == "--Lz") {
            Lz = stod(argv[i+1]);
        } else if (arg == "--T") {
            runtime = stod(argv[i+1]);
            timeProvided = true;
        } else if (arg == "--N") {
            numParticles = stoi(argv[i+1]);
            nProvided = true;
        } else if (arg == "--temp") {
            temperature = stod(argv[i+1]);
            tempProvided = true;
        } else if (arg == "--percent-type1") {
            percent_type1 = stod(argv[i+1]);
        } else if (arg == "--ic-random") {
            icRandomChosen = true;
        } else if (testCaseDict.find(arg) != testCaseDict.end()) {
            string key = arg;
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
        } else if (arg == "--help") {
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

    if ((testCase == true) || (icRandomChosen && nProvided && timeProvided)) {
        cout << "Command Line input well-formatted, carrying on..." << endl;
    } else {
        cout << "Command line input formatted incorrectly, exiting program." << endl;
        exit(1);
    }

    int totalSteps = (runtime / dt) + 1;
    vector<double> timestamps(totalSteps);
    for (int i = 0; i < totalSteps; i++) {
        timestamps[i] = i * dt;
    }

    variableInitialisation(totalSteps, numParticles, X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz);
    //cout << "hi 1" << endl;
    if (icRandomChosen) {
        icRandom(numParticles, Lx, Ly, Lz, percent_type1, X, Y, Z, U, V, W, type);
    } else {
        for (int i = 0; i < numParticles; i++) {
            X[i] = x[i];
            Y[i] = y[i];
            Z[i] = z[i];
            U[i] = u[i];
            V[i] = v[i];
            W[i] = w[i];
        }
    }

    int epsilon[2][2] = { {3,15}, {15,60} };
    int sigma[2][2] = { {1,2}, {2,3} };

    double min_dist = ((X[1] - X[0]) * (X[1] - X[0]) +
                       (Y[1] - Y[0]) * (Y[1] - Y[0]) +
                       (Z[1] - Z[0]) * (Z[1] - Z[0]));

    // Main simulation loop.
    for (int t = 0; t < totalSteps; t++) {
        // Reset forces.
        for (int i = 0; i < numParticles; i++) {
            Fx[i] = 0.0;
            Fy[i] = 0.0;
            Fz[i] = 0.0;
        }
        // Use the flattened serial updater.
        updateVarsSerialFlattened(numParticles, dt, Lx, Ly, Lz, type, tempProvided, temperature, kb,
                                  epsilon, sigma, X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz);
        if (t % 100 == 0) {
            writeToFiles(t, numParticles, timestamps, X, Y, Z, U, V, W, E);
        }
    }
    cout << "minimum distance: " << sqrt(min_dist) << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Runtime: " << duration.count() << " seconds" << endl;

    return 0;
}
