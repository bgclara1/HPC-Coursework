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
#include <cuda_runtime.h>

using namespace std;

__constant__ double s6_table[4] = {0.0, 1.0, 64.0, 729.0};

void variableInitialisation(int totalSteps, int numParticles,
    double** X, double** Y, double** Z,
    double** U, double** V, double** W,
    double** E, double** speed,
    double** Fx, double** Fy, double** Fz,
    double** type)
{
    size_t size = totalSteps * numParticles * sizeof(double);
    cudaMallocManaged((void**)X, size);
    cudaMallocManaged((void**)Y, size);
    cudaMallocManaged((void**)Z, size);
    cudaMallocManaged((void**)U, size);
    cudaMallocManaged((void**)V, size);
    cudaMallocManaged((void**)W, size);
    cudaMallocManaged((void**)E, size);
    cudaMallocManaged((void**)speed, size);
    cudaMallocManaged((void**)Fx, size);
    cudaMallocManaged((void**)Fy, size);
    cudaMallocManaged((void**)Fz, size);

    cudaMemset(*X, 0, size);
    cudaMemset(*Y, 0, size);
    cudaMemset(*Z, 0, size);
    cudaMemset(*U, 0, size);
    cudaMemset(*V, 0, size);
    cudaMemset(*W, 0, size);
    cudaMemset(*E, 0, size);
    cudaMemset(*speed, 0, size);
    cudaMemset(*Fx, 0, size);
    cudaMemset(*Fy, 0, size);
    cudaMemset(*Fz, 0, size);

    size_t typeSize = numParticles * sizeof(double);
    cudaMallocManaged((void**)type, typeSize);
    cudaMemset(*type, 0, typeSize);
}

void icRandom(int numParticles, double Lx, double Ly, double Lz, double percent_type1,
    double* X, double* Y, double* Z,
    double* U, double* V, double* W,
    double* type)
{
    srand(time(0));

    // Initialize positions and velocities.
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
                if (dx * dx + dy * dy + dz * dz < 0.25) {  // 0.5^2 = 0.25
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
    for (int i = 0; i < numParticles; i++) {
        type[i] = 0;
    }
    for (int i = 0; i < numType1; i++) {
        type[i] = 1;
    }
    for (int i = 0; i < numParticles; i++) {
        int j = rand() % numParticles;
        double temp = type[i];
        type[i] = type[j];
        type[j] = temp;
    }
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

__global__
void updateVars(int numParticles, double dt, double Lx, double Ly, double Lz,
    double* type, double temperature, bool tempProvided, double kb,
    const int epsilon[2][2], const int sigma[2][2],
    double* X, double* Y, double* Z,
    double* U, double* V, double* W,
    double* E, double* speed, double* Fx, double* Fy, double* Fz)
{

    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    if (tid >= numParticles) return;
    for (int i = 0; i < numParticles; i++) {
        for (int j = i + 1; j < numParticles; j++) {
            double xij = X[i] - X[j];
            double yij = Y[i] - Y[j];
            double zij = Z[i] - Z[j];
            double rij = xij*xij + yij*yij + zij*zij; // r squared
            int t1 = static_cast<int>(type[i]);
            int t2 = static_cast<int>(type[j]);
            int e = epsilon[t1][t2];
            int s = sigma[t1][t2];

            double inv_r4 = 1.0 / (rij * rij * rij * rij);
            double sigma6_val = s6_table[s] * inv_r4;
            double sigma12 = sigma6_val * sigma6_val * rij;
            double coeff = -24.0 * e * (2.0 * sigma12 - sigma6_val);

            Fx[i] -= xij * coeff;
            Fy[i] -= yij * coeff;
            Fz[i] -= zij * coeff;
            Fx[j] += xij * coeff;
            Fy[j] += yij * coeff;
            Fz[j] += zij * coeff;
        }
    }
    
    for (int i = 0; i < numParticles; i++) {
        int m = (type[i] == 0) ? 1 : 10;
        U[i] += dt * Fx[i] / m;
        V[i] += dt * Fy[i] / m;
        W[i] += dt * Fz[i] / m;
    }
    
    double E_total = 0.0;
    for (int i = 0; i < numParticles; i++) {
        int m = (type[i] == 0) ? 1 : 10;
        double speed2 = U[i]*U[i] + V[i]*V[i] + W[i]*W[i];
        E[i] = 0.5 * m * speed2;
        E_total += E[i];
    }
    
    if (tempProvided) {
        double currentTemp = (2.0 / (3.0 * numParticles * kb)) * E_total;
        double lambda = sqrt(temperature / currentTemp);
        for (int i = 0; i < numParticles; i++) {
            U[i] *= lambda;
            V[i] *= lambda;
            W[i] *= lambda;
        }
    }
    
    for (int i = 0; i < numParticles; i++) {
        X[i] += dt * U[i];
        Y[i] += dt * V[i];
        Z[i] += dt * W[i];
        if (X[i] > Lx) {
            X[i] = 2*Lx - X[i];
            U[i] = -fabs(U[i]);
        }
        if (Y[i] > Ly) {
            Y[i] = 2*Ly - Y[i];
            V[i] = -fabs(V[i]);
        }
        if (Z[i] > Lz) {
            Z[i] = 2*Lz - Z[i];
            W[i] = -fabs(W[i]);
        }
        if (X[i] < 0) {
            X[i] = -X[i];
            U[i] = fabs(U[i]);
        }
        if (Y[i] < 0) {
            Y[i] = -Y[i];
            V[i] = fabs(V[i]);
        }
        if (Z[i] < 0) {
            Z[i] = -Z[i];
            W[i] = fabs(W[i]);
        }
    }
}

void writeToFiles(int t, int numParticles, const vector<double>& timestamps,
                  const double* X, const double* Y, const double* Z,
                  const double* U, const double* V, const double* W,
                  const double* E)
{

    {
        ofstream outfile("output.txt", ios::app);
        outfile << "Time step " << t << "\n";
    }

    {
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
    }

    {
        ofstream posfile("positions.txt", ios::app);
        posfile << "runtime";
        for (int i = 0; i < numParticles; i++) {
            posfile << " x" << i << " y" << i;
        }
        posfile << "\n";
        posfile << std::defaultfloat << timestamps[t];
        for (int i = 0; i < numParticles; i++) {
            posfile << " " << std::fixed << std::setprecision(6) << X[i]
                    << " " << std::fixed << std::setprecision(6) << Y[i];
        }
        posfile << "\n";
    }
}

int main(int argc, char *argv[]) {
    auto start = chrono::high_resolution_clock::now();
    int i = 0;
    double Lx = 20, Ly = 20, Lz = 20, dt = 0.001;
    bool testCase = false, timeProvided = false, nProvided = false, icRandomChosen = false, tempProvided = false;

    ifstream file1("output.txt"); if (file1) { file1.close(); remove("output.txt"); }
    ifstream file2("energy.txt"); if (file2) { file2.close(); remove("energy.txt"); }
    ifstream file3("positions.txt"); if (file3) { file3.close(); remove("positions.txt"); }

    double *X, *Y, *Z, *U, *V, *W, *E, *speed, *Fx, *Fy, *Fz;
    double xij, yij, zij, rij;

    map<string, map<string, vector<double>>> testCaseDict = getTestCases();
    
    double runtime, percent_type1, temperature;
    double kb = 0.8314459920816467;
    int numParticles;
    vector<double> x, y, z, u, v, w;
    double* type;  // This will be allocated in variableInitialisation

    while (i < argc) {
        if (string(argv[i]) == "--Lx") { Lx = stod(argv[i+1]); }
        else if (string(argv[i]) == "--Ly") { Ly = stod(argv[i+1]); }
        else if (string(argv[i]) == "--Lz") { Lz = stod(argv[i+1]); }
        else if (string(argv[i]) == "--T") { runtime = stod(argv[i+1]); timeProvided = true; }
        else if (string(argv[i]) == "--N") { numParticles = stoi(argv[i+1]); nProvided = true; }
        else if (string(argv[i]) == "--temp") { temperature = stod(argv[i+1]); tempProvided = true; }
        else if (string(argv[i]) == "--percent-type1") { percent_type1 = stod(argv[i+1]); }
        else if (string(argv[i]) == "--dt") { dt = stod(argv[++i]); }
        else if (string(argv[i]) == "--ic-random") { icRandomChosen = true; }
        else if (testCaseDict.find(string(argv[i])) != testCaseDict.end()) {
            string key(argv[i]);
            runtime = testCaseDict[key]["runtime"][0];
            numParticles = testCaseDict[key]["numParticles"][0];
            x = testCaseDict[key]["x"];
            y = testCaseDict[key]["y"];
            z = testCaseDict[key]["z"];
            u = testCaseDict[key]["u"];
            v = testCaseDict[key]["v"];
            w = testCaseDict[key]["w"];
            // Instead of assigning a pointer from a temporary vector, we'll copy later.
            testCase = true;
        }
        else if (string(argv[i]) == "--help") {
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

    // Allocate all arrays (including type) in managed memory.
    variableInitialisation(totalSteps, numParticles, &X, &Y, &Z, &U, &V, &W, &E, &speed, &Fx, &Fy, &Fz, &type);

    if (icRandomChosen) {
        icRandom(numParticles, Lx, Ly, Lz, percent_type1, X, Y, Z, U, V, W, type);
    } else if (testCase == true) {
        for (int i = 0; i < numParticles; i++) {
            X[i] = x[i];
            Y[i] = y[i];
            Z[i] = z[i];
            U[i] = u[i];
            V[i] = v[i];
            W[i] = w[i];
            //type[i] = type[i];
        }
    }
    
    
    int epsilon[2][2] = { {3,15}, {15,60} };
    int sigma[2][2] = { {1,2}, {2,3} };

    constexpr int n = 2048;
    int threads = min(256, n);
    int blocks = max(n/256, 1);

    for (int t = 0; t < totalSteps; t++) {
        for (int i = 0; i < numParticles; i++) {
            Fx[i] = 0.0;
            Fy[i] = 0.0;
            Fz[i] = 0.0;           
        }

        updateVars<<<blocks, threads>>>(numParticles, dt, Lx, Ly, Lz, type, temperature, tempProvided, kb,
            epsilon, sigma, X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz);
            
        if (t % 10 == 0) {
            writeToFiles(t, numParticles, timestamps,X,Y,Z,U,V, W,E);
        }
        cudaDeviceSynchronize();
    }

    cudaFree(X);
    cudaFree(Y);
    cudaFree(Z);
    cudaFree(U);
    cudaFree(V);
    cudaFree(W);
    cudaFree(E);
    cudaFree(speed);
    cudaFree(Fx);
    cudaFree(Fy);
    cudaFree(Fz);
    cudaFree(type);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Runtime: " << duration.count() << " seconds" << endl;
    
    return 0;
}
