//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//

//
// Below are the parameter values:
//

#define Tsteps 10000000              // number of Monte-Carlo initial equilibration updates
#define steps  100000000             // number of Monte-Carlo updates
#define stepsPerMeasurement 10       // number of Monte-Carlo updates per measurement
#define beta   1.0                   // inverse temperature
#define alpha  3.0

//
// Below is the list of standard observables:
//

#define MEASURE_H                    // <H>               is measured when this line is not commented
#define MEASURE_H2                   // <H^2>             is measured when this line is not commented
// #define MEASURE_HDIAG                // <H_{diag}>        is measured when this line is not commented
// #define MEASURE_HDIAG2               // <H_{diag}^2>      is measured when this line is not commented
// #define MEASURE_HOFFDIAG             // <H_{offdiag}>     is measured when this line is not commented
// #define MEASURE_HOFFDIAG2            // <H_{offdiag}^2>   is measured when this line is not commented
#define MEASURE_Z_MAGNETIZATION      // Z-magnetization   is measured when this line is not commented
#define MEASURE_Z_MAGNETIZATION2     // Z-magnetization^2 is measured when this line is not commented

//
// Below are the implementation parameters:
//

#define qmax     1000                // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
// #define USE_ABS_WEIGHTS              // uncomment this line to employ absolute values of weights rather than real parts of weights
#define N_CYCLE_TESTS           2000 // number of random tests in test_cycle()
#define N_CYCLE_FIX_ATTEMPTS    5000 // number of attempts to correct cycle in fix_cycle()
#define MAX_FIX_TIME           36000 // maximum time (in seconds) for fixing the cycles. To disable fixing the cycles, set it to zero.
#define GAPS_GEOMETRIC_PARAMETER 0.8 // parameter of geometric distribution for the length of gaps in the cycle completion update
#define EXHAUSTIVE_CYCLE_SEARCH      // comment this line for a more restrictive cycle search
// #define COMPOSITE_UPDATE_BREAK_PROBABILITY  0   // exit composite update at each step with this probability (if this line is uncommented)
#define WORM_UPDATE                  // uncomment this line to employ universal worm update in addition to other updates
// #define WORM_BREAK_PROBABILITY     0 // exit worm update at each step with this probability (if this line is uncommented)
#define MAX_WORM_LENGTH       1000000   // worm exits if its length exceeds MAX_WORM_LENGTH (if this line is uncommented)
// #define COMPOSITE_INSIDE_WORM     // employing composite update inside worm update

//
// Use macros below to enable or disable possibility to checkpoint and restart
//

#define SAVE_UNFINISHED_CALCULATION  // save calculation state to the files "qmc_data_*.dat" prior to exiting when SIGTERM signal is detected
#define SAVE_COMPLETED_CALCULATION   // save detailed data to "qmc_data_*.dat" when calculaiton is completed
#define RESUME_CALCULATION           // attempt to read data from "qmc_data_*.dat" to resume the previous calculation
// #define EXACTLY_REPRODUCIBLE      // uncomment this to always employ the same RNG seeds and reproduce exactly the same results
