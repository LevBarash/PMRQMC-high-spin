-----------------------------------------------------------------------------------------------------------

This program is introduced in the paper: Arman Babakhani, Lev Barash, Itay Hen, A quantum Monte Carlo algorithm for arbitrary high-spin Hamiltonians, preprint arXiv:2503.08039 (2025).

-----------------------------------------------------------------------------------------------------------

Instructions for the Permutation Matrix Representation Quantum Monte Carlo:

1. Check the parameter values of the simulation in the header file "parameters.hpp" such as the number of Monte-Carlo updates and the inverse temperature.

2. Check the list of standard observables in the header file "parameters.hpp"

3. Compile PMRQMC:

		g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp

4. Prepare the Hamiltonian input text file "H.txt".
   Each line corresponds to a summand of the Hamiltonian and contains: "J q_1 sigma_1 power_1 (2s_1+1) q_2 sigma_2 power_2 (2s_2+1) ...", where J is a constant, q_i is a particle index, sigma_i = X, Y, or Z corresponds to the spin matrices, power_i is the exponent to which the spin matrix is raised, and (2s_i+1) defines the particle type. It is also possible to use 1, 2, and 3 instead of X, Y, and Z.

5. Run PMRQMC:

		./PMRQMC.bin

-----------------------------------------------------------------------------------------------------------

Instructions for parallel computing with MPI (optional):

1. Perform the steps 1, 2, and 4 above.
   Note that the parameters "Tsteps" and "steps" are numbers of Monte Carlo updates performed by each MPI process rather than total numbers of Monte Carlo updates.

2. Compile "PMRQMC_mpi.cpp":

		mpicxx -O3 -std=c++11 -o PMRQMC_mpi.bin PMRQMC_mpi.cpp

3. Run the MPI application "PMRQMC_mpi.bin" using mpirun or any MPI-compatible job scheduler.

-----------------------------------------------------------------------------------------------------------
