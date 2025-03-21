
 #include "md.cuh"

 // CUDA kernel to compute forces
 _global_
 void compute_forces_kernel(Particle *d_particles, int numParticles, double *d_epsilon, double *d_sigma) {
     int i = blockIdx.x * blockDim.x + threadIdx.x;
     if (i >= numParticles) return;
 
     Particle &pi = d_particles[i];
     pi.fx = 0.0;
     pi.fy = 0.0;
     pi.fz = 0.0;
 
     for (int j = i + 1; j < numParticles; ++j) {
         Particle &pj = d_particles[j];
 
         double dx = pj.x - pi.x;
         double dy = pj.y - pi.y;
         double dz = pj.z - pi.z;
 
         double r2 = dx * dx + dy * dy + dz * dz;
         if (r2 < 0.0001) continue;
 
         int type1 = pi.type;
         int type2 = pj.type;
         double eps = d_epsilon[type1 * 2 + type2];
         double sig = d_sigma[type1 * 2 + type2];
 
         double sig2 = sig * sig;
         double r2_inv = 1.0 / r2;
         double sig_r2_inv = sig2 * r2_inv;
         double r6 = sig_r2_inv * sig_r2_inv * sig_r2_inv;
         double r12 = r6 * r6;
 
         double f = 24.0 * eps * (2.0 * r12 - r6) * r2_inv;
 
         double fx = f * dx;
         double fy = f * dy;
         double fz = f * dz;
 
         pi.fx -= fx;
         pi.fy -= fy;
         pi.fz -= fz;
 
         pj.fx+= fx;
         pj.fy+= fy;
         pj.fz+= fz;
     }
 }
 
 // Host function to manage CUDA force computation
 void compute_forces_gpu() {
     Particle *d_particles;
     double *d_epsilon, *d_sigma;
 
     int numParticles = particles.size();
     size_t particleSize = numParticles * sizeof(Particle);
     size_t matrixSize = 4 * sizeof(double);
 
     // Allocate GPU memory
     cudaMalloc((void **)&d_particles, particleSize);
     cudaMalloc((void **)&d_epsilon, matrixSize);
     cudaMalloc((void **)&d_sigma, matrixSize);
 
     // Copy data to GPU
     cudaMemcpy(d_particles, particles.data(), particleSize, cudaMemcpyHostToDevice);
     cudaMemcpy(d_epsilon, epsilon, matrixSize, cudaMemcpyHostToDevice);
     cudaMemcpy(d_sigma, sigma, matrixSize, cudaMemcpyHostToDevice);
 
     // Define CUDA grid dimensions
     int threadsPerBlock = 256;
     int blocksPerGrid = (numParticles + threadsPerBlock - 1) / threadsPerBlock;
 
     // Launch kernel
     compute_forces_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles, numParticles, d_epsilon, d_sigma);
     cudaDeviceSynchronize();
 
     // Copy results back to CPU
     cudaMemcpy(particles.data(), d_particles, particleSize, cudaMemcpyDeviceToHost);
 
     // Free GPU memory
     cudaFree(d_particles);
     cudaFree(d_epsilon);
     cudaFree(d_sigma);
 }
 
 using namespace std;
 void simulate_cuda() {
     const int steps = static_cast<int>(final_time / dt);
     const int output_interval = static_cast<int>(0.1 / dt);
     double global_min_dist = std::numeric_limits<double>::max();
 
     for (int step = 0; step <= steps; ++step) {
         double current_time = step * dt;
         if (step % output_interval == 0) {
             write_output(current_time);
         }
 
         if (current_time != final_time) {
             compute_forces_gpu();
             update_velocities();
             update_positions();
             if (test_case!=0){
                 compute_min_distance(global_min_dist);
 
             }
         }
     }
     if (test_case!=0) {
         ofstream min_dist_file("minimum_distance.txt");
         min_dist_file << "Smallest Distance between any 2 particles across All Time: " << global_min_dist << "\n";
         min_dist_file.close();
         data_file.close();
     }
 }