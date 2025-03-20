/**
 * @file cuda.cu
 * @brief  A solver for particle interaction in an enclosed space modelled using Lennard-Jones potential using parallelisation with CUDA.
 *
 * This program simulates the interaction of particles within an enclosed space. 
 * It uses Lennard-Jones potential to model the interactions. 
 * This programme uses a parallel OMP algorithm to solve the problem. 
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
 #include <cuda_runtime.h>
 
 using namespace std;
 
 __constant__ double s6_table[4] = {0.0, 1.0, 64.0, 729.0}; //The values of sigma ^ 6 precalculated and set as a global constant as 
 //the variable never changes and is computationally intensive.
 
 
 /**
  * @brief Allocates and initializes unified memory for simulation arrays.
  *
  * This function allocates managed memory for the updating variables X,Y etc...,
  * Then the memory block contents is set to zero.
  *
  * @param totalSteps Number of simulation steps         
  * @param numParticles Number of particles
  * @param X Pointer to the pointer for x coordinates
  * @param Y Pointer to the pointer for y coordinate
  * @param Z Pointer to the pointer for z coordinates
  * @param U Pointer to the pointer for x velocities
  * @param V Pointer to the pointer for y velocities
  * @param W Pointer to the pointer for z velocities
  * @param E Pointer to the pointer that will hold the kinetic energy
  * @param speed Pointer to the pointer that will hold the speeds
  * @param Fx Pointer to the pointer of the x direction force component of each particle
  * @param Fy Pointer to the pointer of the y direction force component of each particle
  * @param Fz Pointer to the pointer of the z direction force component of each particle
  * @param type Pointer to the pointer that will hold the particle types.
  */
 void variableInitialisation(int totalSteps, int numParticles,
     double** X, double** Y, double** Z,
     double** U, double** V, double** W,
     double** E, double** speed,
     double** Fx, double** Fy, double** Fz,
     double** type)
 {
     double size = totalSteps * numParticles * sizeof(double);
     cudaMallocManaged((void**)X, size);  // basically allocate memory and put address in the pointer that X points to
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
 
     double typeSize = numParticles * sizeof(double);
     cudaMallocManaged((void**)type, typeSize);
     cudaMemset(*type, 0, typeSize);
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
  * @param X Pointer to x coordinates of each particle
  * @param Y Pointer to y coordinates of each particle
  * @param Z Pointer to z coordinates of each particle
  * @param U Pointer to x velocities of each particle
  * @param V Pointer to y velocities of each particle
  * @param W Pointer to z velocities of each particle
  * @param type Pointer to an array of the type of each particle (0 or 1 ie. light or heavy)
  */
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
 
 
 /**
  * @brief Fetches predefined test cases.
  *
  * The brief specifies six ecample cases. This function generates a map of the paramaters for each of the examples.
  *
  * @return A map where keys are test case names (eg. --ic-one-vel) and their values are maps containing their respective simulation parameters.
  */
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
 
 
 
 /**
  * @brief CUDA kernel to update variables for particles.
  *
  * This kernel updates the positions, velocities, energies, and forces on particles based on
  * Lennard-Jones potential equations. Boundary conditions are applied
  * to keep particles within the simulation box. If a temperature is set by the user it is enforced at this point.
  * 
  * It uses CUDA GPU parallelisation to improve the runtime of complex particle simulations.
  *
  * @param min_dist Minimum distance between any two particles
  * @param dt Time step
  * @param numParticles Number of particles
  * @param Lx Length of the container in the x direction
  * @param Ly Length of the container in the y direction
  * @param Lz Length of the container in the z direction
  * @param type Pointer to an array of the type of each particle (0 or 1 ie. light or heavy)
  * @param temperature Chosen simulation temperature
  * @param tempProvided Boolean indicating if the temperature is provided
  * @param kb Boltzmann constant
  * @param epsilon Lennard-Jones Potential epsilon values
  * @param sigma Lennard-Jones Potential sigma values
  * @param X Pointer to x coordinates of each particle
  * @param Y Pointer to y coordinates of each particle
  * @param Z Pointer to z coordinates of each particle
  * @param U Pointer to x velocities of each particle
  * @param V Pointer to y velocities of each particle
  * @param W Pointer to z velocities of each particle
  * @param E Pointer to the pointer that will hold the kinetic energies of each particle
  * @param speed Pointer to the speeds
  * @param Fx Pointer to the x direction force component of each particle
  * @param Fy Pointer to the y direction force component of each particle
  * @param Fz Pointer to the z direction force component of each particle
  */
 __global__   //CUDA kernal updateVars
 void updateVars(int numParticles, double dt, double Lx, double Ly, double Lz,
     double* type, double temperature, bool tempProvided, double kb,
     const int epsilon[2][2], const int sigma[2][2],
     double* X, double* Y, double* Z,
     double* U, double* V, double* W,
     double* E, double* speed, double* Fx, double* Fy, double* Fz)
 {
 
     
     //global thread index = thread id + num threads per block + current block index
     int tid = threadIdx.x + blockDim.x * blockIdx.x; 
     if (tid >= numParticles) return; // exit condition if thread index pointing to a particle that doesn't exist
     for (int i = 0; i < numParticles; i++) {
         for (int j = i + 1; j < numParticles; j++) {
             double xij = X[i] - X[j];
             double yij = Y[i] - Y[j];
             double zij = Z[i] - Z[j];
             double rij = xij*xij + yij*yij + zij*zij; // r squared
             int t1 = type[i];
             int t2 = type[j];
             int e = epsilon[t1][t2];  // finds e and s of the particular particle pair
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
     for (int i = 0; i < numParticles; i++) {        // calculate kinetic energy
         int m = (type[i] == 0) ? 1 : 10;
         double speed2 = U[i]*U[i] + V[i]*V[i] + W[i]*W[i];
         E[i] = 0.5 * m * speed2;
         E_total += E[i];
     }
 
     if (tempProvided) {                                                      //update velocity if temperature is defined by the user
         double currentTemp = (2.0 / (3.0 * numParticles * kb)) * E_total;
         double lambda = sqrt(temperature / currentTemp);
         for (int i = 0; i < numParticles; i++) {
             U[i] *= lambda;
             V[i] *= lambda;
             W[i] *= lambda;
         }
     }
     
     for (int i = 0; i < numParticles; i++) { //apply boundary conditions
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
  * @param X Pointer to the pointer for x coordinates
  * @param Y Pointer to the pointer for y coordinate
  * @param Z Pointer to the pointer for z coordinates
  * @param U Pointer to the pointer for x velocities
  * @param V Pointer to the pointer for y velocities
  * @param W Pointer to the pointer for z velocities
  * @param E Pointer to the pointer that will hold the kinetic energy
  */
 void writeToFiles(int t, int numParticles, const vector<double>& timestamps,
                   const double* X, const double* Y, const double* Z,
                   const double* U, const double* V, const double* W,
                   const double* E)
 {
 
     {
         ofstream energyfile("energy.txt", ios::app);  // write time stamp and KE to kinetic energy file
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
         ofstream posfile("positions.txt", ios::app);   // write time stamp, x and y position to position file
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
     double *X, *Y, *Z, *U, *V, *W, *E, *speed, *Fx, *Fy, *Fz;
     double xij, yij, zij, rij;
 
     map<string, map<string, vector<double>>> testCaseDict = getTestCases();
     
     double runtime, percent_type1, temperature;
     double kb = 0.8314459920816467;
     int numParticles;
     vector<double> x, y, z, u, v, w;
     double* type;  // This will be allocated in variableInitialisation
 
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
     
     if ((testCase == true) || (icRandomChosen && nProvided && timeProvided)) {  // check if args are valid
         cout << "Command Line input well-formatted, carrying on..." << endl;
     } else {
         cout << "Command line input formatted incorrectly, exiting program." << endl;
         exit(1);
     }
     
     int totalSteps = (runtime / dt) + 1;
     vector<double> timestamps(totalSteps);
     for (int i = 0; i < totalSteps; i++) {
         timestamps[i] = i * dt;  // vector going from 0 to time T in increments dt
     }
 
     // allocate vars in managed memory
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
         }
     }
     
     
     int epsilon[2][2] = { {3,15}, {15,60} }; //initialise epsilon and sigma as stated in brief
     int sigma[2][2] = { {1,2}, {2,3} };
 
     constexpr int n = 2048; //compile time constant. n value found to be most optimal through trial and error
     int threads = min(256, n);
     int blocks = max(n/256, 1);
 
     for (int t = 0; t < totalSteps; t++) {
         for (int i = 0; i < numParticles; i++) {
             Fx[i] = 0.0;
             Fy[i] = 0.0;
             Fz[i] = 0.0;           
         }
 
         updateVars<<<blocks, threads>>>(numParticles, dt, Lx, Ly, Lz, type, temperature, tempProvided, kb,
             epsilon, sigma, X, Y, Z, U, V, W, E, speed, Fx, Fy, Fz);  //run CUDA kernal updateVars 
             
         if (t % 100 == 0) {
             writeToFiles(t, numParticles, timestamps,X,Y,Z,U,V, W,E);
         }
         cudaDeviceSynchronize(); //wait for GPU ops to complete before next time step iteration
     }
 
     cudaFree(X); // release managed memory
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
 