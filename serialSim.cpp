/**
 * @file serialSim.cpp
 * @brief A serial solver for particle interaction in an enclosed space modelled using Lennard-Jones potential.
 *
 * This program simulates the interaction of particles within an enclosed space. 
 * It uses Lennard-Jones potential to model the interactions. 
 * This programme uses a serial algorithm to solve the problem. 
 * The programme allows the user to specify the dimensions of the box, the number of particles, the ratio of heavy to light particles,
 * the temperature, the duration of the simulation and the timestep. 
 * The user may also chose to run 6 example simulations.
 */

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

const double s6_table[4] = {0.0, 1.0, 64.0, 729.0}; //The values of sigma ^ 6 precalculated and set as a global constant as 
                                                    //the variable never changes and is computationally intensive at a high num of particles

/**
 * @brief Initializes particle variables.
 *
 * This function initializes the variables used for particle properties by setting them
 * to totalSteps * numParticles with initial values of 0.0. The variables that get updated (X, Y etc.)
 * are stored as vector doubles such that they contain the information of each particle in the system.
 *
 * @param totalSteps Total number of steps in the simulation
 * @param numParticles Number of particles 
 * @param X x coordinate of each particle
 * @param Y y coordinate of each particle
 * @param Z z coordinate of each particle
 * @param U x velocity of each particle
 * @param V y velocity of each particle
 * @param W z velocity of each particle
 * @param E Kinetic Energy of each particle
 * @param speed speed of each particle
 * @param Fx x direction force component of each particle
 * @param Fy y direction force component of each particle
 * @param Fz z direction force component of each particle
 */
void variableInitialisation(int totalSteps, int numParticles,
    vector<double>& X, vector<double>& Y, vector<double>& Z,
    vector<double>& U, vector<double>& V, vector<double>& W,
    vector<double>& E, vector<double>& speed,
    vector<double>& Fx, vector<double>& Fy, vector<double>& Fz)
{
    X.resize(totalSteps * numParticles, 0.0);
    Y.resize(totalSteps * numParticles, 0.0);
    Z.resize(totalSteps * numParticles, 0.0);
    U.resize(totalSteps * numParticles, 0.0);
    V.resize(totalSteps * numParticles, 0.0);
    W.resize(totalSteps * numParticles, 0.0);
    E.resize(totalSteps * numParticles, 0.0);
    speed.resize(totalSteps * numParticles, 0.0);
    Fx.resize(totalSteps * numParticles, 0.0);
    Fy.resize(totalSteps * numParticles, 0.0);
    Fz.resize(totalSteps * numParticles, 0.0);
}

/**
 * @brief Initializes particles with random positions and velocities for when the user inputs --ic-random as a command line parameter.
 *
 * This function assigns random initial positions and velocities to particles. The posititions must be within the default (20) or requested dimensions
 * of the container. The velocities are generated between -0.5 and 0.5.
 * Additionally, it randomly assigns particle types according to the default (10%) percentage of heavy particles or the percentage requested by the user.
 *
 * @param numParticles Number of particles
 * @param Lx Length of the container in the x direction
 * @param Ly Length of the container in the y direction
 * @param Lz Length of the container in the z direction
 * @param percent_type1 Percentage of particles of type 1 (heavy)
 * @param X x coordinate of each particle
 * @param Y y coordinate of each particle
 * @param Z z coordinate of each particle
 * @param U x velocity of each particle
 * @param V y velocity of each particle
 * @param W z velocity of each particle
 * @param type type of each particle (0 or 1 ie. light or heavy)
 */
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
            for (int j = 0; j < i; j++) {   // check the particles aren't initialised too close together
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

/**
 * @brief Fetches predefined test cases.
 *
 * The brief specifies six ecample cases. This function generates a map of the paramaters for each of the examples.
 *
 * @return A map where keys are test case names (eg. --ic-one-vel) and their values are maps containing their respective simulation parameters.
 */
map<string, map<string, vector<double>>> getTestCases() {   // nested dictionaries/maps to store the test cases
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


/**
 * @brief Updates particle position and velocity variables for each time step.
 *
 * This function updates the positions, velocities, energies, and forces on particles based on
 * Lennard-Jones potential equations. Boundary conditions are applied
 * to keep particles within the simulation box. If a temperature is set by the user it is enforced at this point.
 *
 * @param min_dist Minimum distance between any two particles
 * @param dt Time step
 * @param numParticles Number of particles
 * @param Lx Length of the container in the x direction
 * @param Ly Length of the container in the y direction
 * @param Lz Length of the container in the z direction
 * @param type Type of each particle (0 or 1 ie. light or heavy)
 * @param temperature Chosen simulation temperature
 * @param tempProvided Boolean indicating if the temperature is provided
 * @param kb Boltzmann constant
 * @param epsilon Lennard-Jones Potential epsilon values
 * @param sigma Lennard-Jones Potential sigma values
 * @param X x coordinate of each particle
 * @param Y y coordinate of each particle
 * @param Z z coordinate of each particle
 * @param U x velocity of each particle
 * @param V y velocity of each particle
 * @param W z velocity of each particle
 * @param E Kinetic Energy of each particle
 * @param speed Speed of each particleß
 * @param xij Difference in x position between two particles 
 * @param yij Difference in y position between two particles
 * @param zij Difference in z position between two particles 
 * @param rij Distance between particles squared (!)
 * @param Fx x direction force component of each particle
 * @param Fy y direction force component of each particle
 * @param Fz z direction force component of each particle
 */
void updateVars(double min_dist, int numParticles, double dt, double Lx, double Ly, double Lz,
    vector<double>& type, double temperature, bool tempProvided, double kb,
    const int epsilon[2][2], const int sigma[2][2],
    vector<double>& X, vector<double>& Y, vector<double>& Z,
    vector<double>& U, vector<double>& V, vector<double>& W,
    vector<double>& E, vector<double>& speed, double& xij,
    double& yij, double& zij, double& rij, vector<double>& Fx,
    vector<double>& Fy, vector<double>& Fz)
{
    double s6 = 0;
    for (int i = 0; i < numParticles; i++) {
        for (int j = i + 1; j < numParticles; j++) {        // builds upper triangle of matrix to avoid the double calculation
            xij = X[i] - X[j];
            yij = Y[i] - Y[j];
            zij = Z[i] - Z[j];
            rij = xij * xij + yij * yij + zij * zij; // distance squared

            int t1 = type[i];
            int t2 = type[j];
            int e = epsilon[t1][t2];        // finds e and s of the particular particle pair
            int s = sigma[t1][t2];

            double inv_r4 = 1.0 / (rij * rij * rij * rij); //calculation split up and not including pow for optimisation reasons
            double sigma6 = s6_table[s] * inv_r4;
            double sigma12 = sigma6 * sigma6 * rij;
            double coeff = -24.0 * e * (2.0 * sigma12 - sigma6);

            Fx[i] -= xij * coeff;       //calculate the net forces
            Fy[i] -= yij * coeff;
            Fz[i] -= zij * coeff;
            Fx[j] += xij * coeff;   // the opposite sign is applied to the j indices as this represents the missing lower triangle of the matrix
            Fy[j] += yij * coeff;
            Fz[j] += zij * coeff;
        }
    }
    for (int i = 0; i < numParticles; i++) {
        int m = (type[i] == 0) ? 1 : 10; // if true pick 1 else 10
        U[i] += dt * Fx[i] / m; //update velocities
        V[i] += dt * Fy[i] / m;
        W[i] += dt * Fz[i] / m;
    }
    double E_total = 0.0;
    for (int i = 0; i < numParticles; i++) {                // calculate kinetic energy
        int m = (type[i] == 0) ? 1 : 10;
        double speed2 = (U[i] * U[i] + V[i] * V[i] + W[i] * W[i]);
        E[i] = 0.5 * m * speed2;
        E_total += E[i];
    }
    if (tempProvided) {                 //update velocity if temperature is defined by the user
        double currentTemp = (2.0 / (3.0 * numParticles * kb)) * E_total;
        double lambda = sqrt(temperature / currentTemp);
        for (int i = 0; i < numParticles; i++) {
            U[i] *= lambda;
            V[i] *= lambda;
            W[i] *= lambda;
        }
    }
    for (int i = 0; i < numParticles; i++) {
        X[i] = X[i] + dt * U[i];                // update position
        Y[i] = Y[i] + dt * V[i];
        Z[i] = Z[i] + dt * W[i];

        if (X[i] > Lx) {                            // Apply Boundary conditions
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

/**
 * @brief Writes simulation data to output files.
 *
 * This function writes particle positions, velocities, kinetic energy, and timestamps to the files energy.txt and positions.txt.
 * energy.txt containes timestamp and kinetic energy of each particle.alignas. eg. Time E1 E2 E3 ... 
 * positions.txt containes the timestamp and x and y position of each particle. eg Time X1 Y1 X2 Y2 X3 Y3...
 * 
 *
 * @param t Current time step index
 * @param numParticles Number of particles
 * @param timestamps timestamps from 0 to the time set by the user in increments dt set by the user
 * @param X x coordinate of each particle
 * @param Y y coordinate of each particle
 * @param Z z coordinate of each particle
 * @param U x velocity of each particle
 * @param V y velocity of each particle
 * @param W z velocity of each particle
 * @param E Kinetic Energy of each particle
 */
void writeToFiles(int t, int numParticles, const vector<double>& timestamps,
    const vector<double>& X, const vector<double>& Y,
    const vector<double>& Z, const vector<double>& U,
    const vector<double>& V, const vector<double>& W,
    const vector<double>& E)
{

    {
        ofstream energyfile("energy.txt", ios::app);        // write time stamp and KE to kinetic energy file
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
        ofstream posfile("positions.txt", ios::app);        // write time stamp, x and y position to position file
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

/**
 * @brief Main simulation program.
 *
 * This function reads command line arguments, initializes the variables, runs the simulation for each timestep,
 * and writes to output files.
 *
 * @param argc Number of arguments provided
 * @param argv Arguments provided stored as strings
 * @return Exit value
 */
int main(int argc, char *argv[]) { // read cmd args w main params.
    auto start = chrono::high_resolution_clock::now();              // start runtime clock
    int i = 0;
    double Lx = 20;         // initialise default params
    double Ly = 20;
    double Lz = 20;
    double dt = 0.001;
    bool testCase = false;
    bool timeProvided = false;
    bool nProvided = false;
    bool icRandomChosen = false;
    bool tempProvided = false;

    ifstream file1("output.txt"); // close files in case make clean isn't run. Function write to file appends so it's worth doing just in case.
    if (file1) {
        file1.close();
        remove("output.txt");
    }
    ifstream file2("energy.txt");
    if (file2) {
        file2.close();
        remove("energy.txt");
    }
    ifstream file3("positions.txt");
    if (file3) {
        file3.close();
        remove("positions.txt");
    }

    vector<double> X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz;
    double xij, yij, zij, rij, dPhi_dx, dPhi_dy, dPhi_dz;

    map<string, map<string, vector<double>>> testCaseDict = getTestCases();        // get the params for the test cases 1 to 6
    double runtime, percent_type1, temperature;
    double kb = 0.8314459920816467;         // Boltzman constant
    int numParticles;
    vector<double> x, y, z, u, v, w, type;
    while (i < argc) {                              // save args given by user into relevant variables
        if (string(argv[i]) == "--Lx") {
            Lx = stod(argv[i + 1]);
        } else if (string(argv[i]) == "--Ly") {
            Ly = stod(argv[i + 1]);
        } else if (string(argv[i]) == "--Lz") {
            Lz = stod(argv[i + 1]);
        } else if (string(argv[i]) == "--T") {
            runtime = stod(argv[i + 1]);
            timeProvided = true;
        } else if (string(argv[i]) == "--N") {
            numParticles = stoi(argv[i + 1]);
            nProvided = true;
        } else if (string(argv[i]) == "--temp") {
            temperature = stod(argv[i + 1]);
            tempProvided = true;
        } else if (string(argv[i]) == "--percent-type1") {
            percent_type1 = stod(argv[i + 1]);
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

    if ((testCase == true) || (icRandomChosen == true && nProvided == true && timeProvided == true)) {      // check if args are valid
        cout << "Command Line input well-formatted, carrying on..." << endl;
    } else {
        cout << "Command line input formatted incorrectly, exiting program." << endl;
        exit(1);
    }

    int totalSteps = (runtime / dt) + 1;
    vector<double> timestamps(totalSteps);
    for (int i = 0; i < totalSteps; i++) {
        timestamps[i] = i * dt;             // vector going from 0 to time T in increments dt
    }

    variableInitialisation(totalSteps, numParticles, X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz);       //initialise variables - given a separate function for main function readability

    if (icRandomChosen) {
        icRandom(numParticles, Lx, Ly, Lz, percent_type1, X, Y, Z, U, V, W, type);
    } else {
        for (int i = 0; i < numParticles; i++) {
            X[i] = x[i]; // the lower case x,y,z etc are how the vars are stored in testCaseDict
            Y[i] = y[i];
            Z[i] = z[i];
            U[i] = u[i];
            V[i] = v[i];
            W[i] = w[i];
        }
    }

    int epsilon[2][2] = { {3,15}, {15,60} };    //initialise epsilon and sigma as stated in brief
    int sigma[2][2] = { {1,2}, {2,3} };

    double min_dist = ((X[1] - X[0]) * (X[1] - X[0]) +      // calculate an initial minimum distance to compare against using first 2 particles at t=0
                       (Y[1] - Y[0]) * (Y[1] - Y[0]) +
                       (Z[1] - Z[0]) * (Z[1] - Z[0]));

    for (int t = 0; t < totalSteps; t++) {
        for (int i = 0; i < numParticles; i++) {
            Fx[i] = 0.0;
            Fy[i] = 0.0;
            Fz[i] = 0.0;
        }
        updateVars(min_dist, numParticles, dt, Lx, Ly, Lz, type, temperature,
                   tempProvided, kb, epsilon, sigma, X, Y, Z, U, V, W, E, speed,
                   xij, yij, zij, rij, Fx, Fy, Fz);
        if (t % 100 == 0) {     // save to file every 100 timestamps to balance performance and fineness of file data
            writeToFiles(t, numParticles, timestamps, X, Y, Z, U, V, W, E);
        }
    }
    cout << "minimum distance: " << sqrt(min_dist) << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Runtime: " << duration.count() << " seconds" << endl;

    return 0;
}
