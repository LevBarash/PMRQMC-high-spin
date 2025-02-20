//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include<iostream>
#include<iomanip>
#include<complex>
#include<random>
#include<cstdlib>
#include<algorithm>
#include<csignal>
#include"DatatoPMR.hpp"
#include"ModnLA.hpp"
#include"divdiff.hpp"
#include"parameters.hpp"  // parameters of the simulation such as the number of Monte-Carlo updates

#define measurements (steps/stepsPerMeasurement)

#ifdef EXHAUSTIVE_CYCLE_SEARCH
	#define rmin 0                        // length r of sub-sequence is chosen randomly between rmin and rmax
	#define rmax cycle_max_len
	#define lmin r                        // cycle lengths must be between lmin and lmax
	#define lmax cycle_max_len
#else
	#define rmin (cycle_min_len-1)/2
	#define rmax (cycle_max_len+1)/2
	#define lmin 2*r-1
	#define lmax 2*r+1
#endif

int N;
int Nop;
int Ncycles;
PDdata CDPdata;
vector<vector<int>> cycles;

static std::random_device rd;
static std::mt19937 rng;
static std::uniform_int_distribution<> dice2(0,1);
static std::uniform_int_distribution<> diceS;
static std::uniform_int_distribution<> diceSminus1;
static std::uniform_int_distribution<> diceN;
static std::uniform_int_distribution<> diceNop;
static std::uniform_real_distribution<> val(0.0,1.0);
static std::geometric_distribution<> geometric_int(GAPS_GEOMETRIC_PARAMETER);

ExExFloat beta_pow_factorial[qmax]; // contains the values (-beta)^q / q!
double factorial[qmax]; // contains the values q!
vector<int> cycle_len;
vector<int> cycles_used;
vector<int> cycles_used_backup;
int cycle_min_len, cycle_max_len, found_cycles, min_index, max_index;

#ifndef MEASURE_CUSTOM_OBSERVABLES
#define Nobservables 0
#endif

const int N_all_observables = Nobservables + 8;
int valid_observable[N_all_observables];

unsigned long long bin_length = measurements / Nbins, bin_length_old;
double in_bin_sum[N_all_observables];
double bin_mean[N_all_observables][Nbins];
double in_bin_sum_sgn;
double bin_mean_sgn[Nbins];

int q;
int qmax_achieved=0;

divdiff* d;

std::vector<int> lattice;
std::vector<int> z;
std::vector<int> P;

int Sq[qmax];	// list of q operator indices
int Sq_backup[qmax];
int Sq_subseq[qmax];
int Sq_gaps[qmax];
double Energies[qmax+1];
double Energies_backup[qmax+1];
int eoccupied[qmax+1];
double currEnergy;
std::complex<double> currD, old_currD;
std::complex<double> currD_partial[qmax];

int abs_weights = 
#ifdef USE_ABS_WEIGHTS
	1;
#else
	0;
#endif

#define REALW(x)  ( abs_weights ? std::abs(x) : std::real(x) )

#ifdef WORM_UPDATE
int q_saved;
std::vector<int> lattice_backup;
int Sq_saved[qmax];
int distance_from_identity = 0; int worm_length;
int saved_distance_from_identity;
#endif

ExExFloat one, currWeight; int TstepsFinished = 0;
unsigned long long step = 0;
unsigned long long measurement_step = 0;

double Svalue; int twoSplus1;
vector<int> operator_count; // contains numbers of permutation operators in Sq_subseq of length r (or in its inverse).
double permutation_count; // multinomial coefficient (the number of distinct permutations) corresponding to operator_count
vector<int> InvPermutationIndex;
vector<int> MainPermutationIndex;

unsigned int rng_seed;
double meanq = 0;
double maxq = 0;
double start_time;
int save_data_flag = 0, mpi_rank = 0, mpi_size = 0, resume_calc = 0;

#define fout(obj) { foutput.write((char *)& obj, sizeof(obj)); }
#define fout_vec(obj) { std::size_t tt = obj.size(); fout(tt); for(auto & val : obj) fout(val); }
#define fout_vec_vec(obj) { std::size_t tt = obj.size(); fout(tt); for(auto & val : obj) fout_vec(val); }
#define fout_vec_vec_vec(obj) { std::size_t tt = obj.size(); fout(tt); for(auto & val : obj) fout_vec_vec(val); }
#define fin(obj)  { finput.read((char *)& obj, sizeof(obj)); }
#define fin_vec(obj) { std::size_t tt; fin(tt); obj.resize(tt); for(auto & val : obj) fin(val); }
#define fin_vec_vec(obj) { std::size_t tt; fin(tt); obj.resize(tt); for(auto & val : obj) fin_vec(val); }
#define fin_vec_vec_vec(obj) { std::size_t tt; fin(tt); obj.resize(tt); for(auto & val : obj) fin_vec_vec(val); }

void save_QMC_data(int printout = 1){
	if(printout) cout<<"SIGTERM signal detected. Saving unfinished calculation...";
	char fname[100]; if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
	ofstream foutput(fname,ios::binary); double elapsed;
	fout(CDPdata.twoSplusOne); fout(CDPdata.NumOfParticles);
	fout_vec_vec(CDPdata.Permutations); fout_vec_vec_vec(CDPdata.Diagonals); fout_vec_vec(CDPdata.Coefficients);
	fout_vec_vec(CDPdata.D0); fout_vec(CDPdata.D0Coeffs);
	fout(N); fout(Nop); fout(Ncycles); fout_vec_vec(cycles); fout(rng); fout(bin_length);
	fout_vec(cycle_len); fout_vec(cycles_used); fout_vec(cycles_used_backup); fout(cycle_min_len); fout(cycle_max_len);
	fout(found_cycles); fout(min_index); fout(max_index); fout(TstepsFinished);
	fout(in_bin_sum); fout(bin_mean); fout(in_bin_sum_sgn); fout(bin_mean_sgn); 
	fout(q); fout(qmax_achieved); fout_vec(lattice); fout_vec(z); fout_vec(P); fout(Sq); fout(Sq_backup); fout(Sq_subseq); fout(Sq_gaps);
	fout(Energies); fout(Energies_backup); fout(eoccupied); fout(currEnergy); fout(valid_observable);
	fout(currD); fout(old_currD); fout(currD_partial); fout(one); fout(currWeight); fout(step); fout(measurement_step); 
	fout(Svalue); fout(twoSplus1); fout_vec(operator_count); fout(permutation_count);
	fout_vec(InvPermutationIndex); fout_vec(MainPermutationIndex); fout(rng_seed); fout(meanq); fout(maxq); fout(abs_weights);
#ifdef WORM_UPDATE
	fout(q_saved); fout_vec(lattice_backup); fout(Sq_saved); fout(distance_from_identity); fout(saved_distance_from_identity);
#endif
#ifdef MPI_VERSION
	elapsed = MPI_Wtime() - start_time;
#else
	elapsed = (double)clock() / CLOCKS_PER_SEC - start_time;
#endif
	fout(elapsed);
	foutput.close(); if(printout) cout<<"done"<<endl; fflush(stdout);
}

int check_QMC_data(){
#ifdef RESUME_CALCULATION
	int i,r,g; char fname[100];
	if(mpi_rank==0 && mpi_size>0){
		r = 1;
		for(i=0;i<mpi_size;i++){
			sprintf(fname,"qmc_data_%d.dat",mpi_rank);
			std::ifstream finput(fname,std::ios::binary);
			g = finput.good(); if(g) finput.close(); r = r && g;
		}
	} else{
		if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
		std::ifstream finput(fname,std::ios::binary);
		r = finput.good(); if(r) finput.close();
	}
#else
	int r = 0;
#endif
	return r;
}

void load_QMC_data(){
	char fname[100]; if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
	std::ifstream finput(fname,std::ios::binary); double elapsed;
	if(finput.good()){
		fin(CDPdata.twoSplusOne); fin(CDPdata.NumOfParticles);
		fin_vec_vec(CDPdata.Permutations); fin_vec_vec_vec(CDPdata.Diagonals); fin_vec_vec(CDPdata.Coefficients);
		fin_vec_vec(CDPdata.D0); fin_vec(CDPdata.D0Coeffs);
		fin(N); fin(Nop); fin(Ncycles); fin_vec_vec(cycles); fin(rng); fin(bin_length_old);
		fin_vec(cycle_len); fin_vec(cycles_used); fin_vec(cycles_used_backup); fin(cycle_min_len); fin(cycle_max_len);
		fin(found_cycles); fin(min_index); fin(max_index); fin(TstepsFinished);
		fin(in_bin_sum); fin(bin_mean); fin(in_bin_sum_sgn); fin(bin_mean_sgn); 
		fin(q); fin(qmax_achieved); fin_vec(lattice); fin_vec(z); fin_vec(P); fin(Sq); fin(Sq_backup); fin(Sq_subseq); fin(Sq_gaps);
		fin(Energies); fin(Energies_backup); fin(eoccupied); fin(currEnergy); fin(valid_observable);
		fin(currD); fin(old_currD); fin(currD_partial); fin(one); fin(currWeight); fin(step); fin(measurement_step); 
		fin(Svalue); fin(twoSplus1); fin_vec(operator_count); fin(permutation_count);
		fin_vec(InvPermutationIndex); fin_vec(MainPermutationIndex); fin(rng_seed); fin(meanq); fin(maxq); fin(abs_weights);
#ifdef WORM_UPDATE
		fin(q_saved); fin_vec(lattice_backup); fin(Sq_saved); fin(distance_from_identity); fin(saved_distance_from_identity);
#endif
		fin(elapsed); if(finput.gcount()==0) elapsed = 0;
		finput.close(); start_time -= elapsed;
		if(mpi_size > 0){
			if(TstepsFinished){
				cout<<"MPI process No. "<< mpi_rank <<" loaded unfinished calculation: all Tsteps completed, main steps completed: "
				    <<measurement_step*stepsPerMeasurement<<" of "<<steps<<"."<<endl;
			} else cout <<"MPI process No. "<< mpi_rank <<" loaded unfinished calculation: Tsteps completed: "<<step<<" of "<<Tsteps<<"."<<endl;
		} else{
			if(TstepsFinished){
				cout<<"Loaded unfinished calculation: all Tsteps completed, main steps completed: "
				    <<measurement_step*stepsPerMeasurement<<" of "<<steps<<"."<<endl;
			} else cout <<"Loaded unfinished calculation: Tsteps completed: "<<step<<" of "<<Tsteps<<"."<<endl;
		}
	} else{
		if(mpi_size>0) cout<<"MPI process No. "<<mpi_rank<<": error opening file "<<fname<<endl;
		else cout<<"Error opening file "<<fname<<endl; fflush(stdout);
#ifdef MPI_VERSION
		MPI_Abort(MPI_COMM_WORLD,1);
#else
		exit(1);
#endif
	}
	if(bin_length != bin_length_old){ // It is allowed to increase steps by an integer number of times for a completed calculation
		if(bin_length > 0 && bin_length % bin_length_old == 0){ // All other parameters should remain unchanged
			double sum; int i, j, o, m = bin_length / bin_length_old, curr_bins = Nbins/m; // merging bins together
			for(o=0;o<N_all_observables;o++) for(i=0;i<curr_bins;i++){
				sum = 0; for(j=0;j<m;j++) sum += bin_mean[o][m*i+j]; sum /= m;
				bin_mean[o][i] = sum;
			}
			for(i=0;i<curr_bins;i++){
				sum = 0; for(j=0;j<m;j++) sum += bin_mean_sgn[m*i+j]; sum /= m;
				bin_mean_sgn[i] = sum;
			}
			for(o=0;o<N_all_observables;o++){
				sum = 0; for(i=curr_bins*m;i<Nbins;i++) sum += bin_mean[o][i]*bin_length_old;
				in_bin_sum[o] += sum;
			}
			sum = 0; for(i=curr_bins*m;i<Nbins;i++) sum += bin_mean_sgn[i]*bin_length_old;
			in_bin_sum_sgn += sum;
		} else{
			cout << "Error: bin_length = " << bin_length <<" is not divisible by bin_length_old = " << bin_length_old << endl; fflush(stdout);
#ifdef MPI_VERSION
			MPI_Abort(MPI_COMM_WORLD,1);
#else
			exit(1);
#endif
		}
	}
}

double Dz(int j){ return Svalue-j;}
double Dplus(int j){ return sqrt((twoSplus1-j-1)*(j+1))/2; }
double Dk(int k, int j){ return Dplus((j+k+twoSplus1) % twoSplus1); }
double Dzk(int k, int j){ return Dz((j+k+twoSplus1) % twoSplus1); }

double CalcEnergy(){ // calculate the energy <z | D_0 | z> of a given configuration of spins
	int i, j, index, D0size, Dlength; TotalDiag t;
	std::complex<double> sum = 0; double summand;
	D0size = CDPdata.D0.size();
        for(i=0;i<D0size;i++){
		summand = 1; Dlength = CDPdata.D0[i].size();
		for(j=0;j<Dlength;j++){
			t = CDPdata.D0[i][j];
			index = lattice[t.particle-1];
			summand *= t.ztype ? Dzk(t.k,index) : Dk(t.k,index);
		}
		sum += CDPdata.D0Coeffs[i] * summand;
	}
	return sum.real();
}

std::complex<double> calc_d(int k){ // calculate d_k = <z | D_k | z> for the current configuration of spins
	int i, j, index, Dsize, Dlength; TotalDiag t;
	std::complex<double> sum = 0; double summand;
	Dsize = CDPdata.Diagonals[k].size();
        for(i=0;i<Dsize;i++){
		summand = 1; Dlength = CDPdata.Diagonals[k][i].size();
		for(j=0;j<Dlength;j++){
			t = CDPdata.Diagonals[k][i][j];
			index = lattice[t.particle-1];
			summand *= t.ztype ? Dzk(t.k,index) : Dk(t.k,index);
		}
		sum += CDPdata.Coefficients[k][i] * summand;
	}
	return sum;
}

void ApplyOperator(int k){
	pair<int,int> t; int i; int Perm_size = CDPdata.Permutations[k].size();
	for(i=0;i<Perm_size;i++){
		t = CDPdata.Permutations[k][i];
		lattice[t.first-1] = (lattice[t.first-1] - t.second + twoSplus1) % twoSplus1;
	}
}

void GetEnergies(){
#ifdef WORM_UPDATE
	z = lattice;
	int dist=0;
#endif
	currD = currD_partial[0] = 1;
	for(int i=0;i<q;i++){
		Energies[i] = CalcEnergy();
		ApplyOperator(Sq[i]);
		currD *= calc_d(Sq[i]);
		currD_partial[i+1] = currD;
	}
	Energies[q] = CalcEnergy(); currEnergy = Energies[0];
#ifdef WORM_UPDATE
	for(int i=0;i<N;i++) if(lattice[i] != z[i]) dist++;
	distance_from_identity = dist; currD /= exp(dist*alpha);
	lattice = z;
#endif
}

#include"cycles.hpp"

ExExFloat GetWeight(){
	d->CurrentLength=0; GetEnergies();
	for(int i=0;i<=q;i++) d->AddElement(-beta*Energies[i]);
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);
}

ExExFloat UpdateWeight(){
	int i, j, notfound, n=d->CurrentLength; double value;
	GetEnergies(); memset(eoccupied,0,(q+1)*sizeof(int));
	for(i=0;i<n;i++){
	        notfound = 1; value = d->z[i];
		for(j=0;j<=q;j++) if(eoccupied[j]==0 && value == -beta*Energies[j]){ notfound = 0; break; }
		if(notfound) break; eoccupied[j] = 1;
	}
	if(i==0) d->CurrentLength=0; else while(n>i){ d->RemoveElement(); n--; }
	j=0; while(i<=q){ while(eoccupied[j]) j++; d->AddElement(-beta*Energies[j++]); i++; }
        return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);
}

ExExFloat UpdateWeightReplace(double removeEnergy, double addEnergy){
	if(removeEnergy != addEnergy){
		if(d->RemoveValue(-beta*removeEnergy)) d->AddElement(-beta*addEnergy); else{
			std::cout << "Error: energy not found" << std::endl; exit(1);
		}
	}
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD); // use this value only when the values of q and currD are correct
}

ExExFloat UpdateWeightDel(double removeEnergy1, double removeEnergy2){
	if(d->RemoveValue(-beta*removeEnergy1) && d->RemoveValue(-beta*removeEnergy2))
		return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);  // use this value only when the values of q and currD are correct
	else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
}

ExExFloat UpdateWeightIns(double addEnergy1, double addEnergy2){
	d->AddElement(-beta*addEnergy1); d->AddElement(-beta*addEnergy2);
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);    // use this value only when the values of q and currD are correct
}

int NoRepetitionCheck(int* sequence, int r){ // check for absence of repetitions in a sequence of length r
	int i,j,rep = 1;
	for(i=0;i<r && rep;i++) for(j=0;j<i;j++) if(sequence[j]==sequence[i]){ rep = 0; break;}
	return rep;
}

void PickSubsequence(int r){ // randomly picks a sequential sub-sequence of length r from Sq
	int i,j,m; m = int(val(rng)*(q-r+1)); // m is random integer between 0 and q-r
	for(i=0;i<r;i++) Sq_subseq[i] = Sq[i+m];
	min_index = m; max_index = m+r-1;
}

int vector_sum(vector<int> v){
	int i, sum = 0;
	for(i=0;i<v.size();i++) sum += v[i];
	return sum;
}

void lattice_flip(int p){
	int spin_new = twoSplus1==2 ? 0 : diceSminus1(rng);
	if(spin_new >= lattice[p]) spin_new++;
	lattice[p] = spin_new;
}

int FindCycles(int r, int inv){  // find all cycles of length between lmin and lmax, each containing all operators of Sq_subseq of length r.
	int i,j,k,l,not_contained;  // if inv == 1, then containing all operators of the inverse of Sq_subseq of length r.
	vector<int> found_cycle_list(Ncycles);
	fill(operator_count.begin(), operator_count.end(), 0);
	if(inv)
		for(j=0;j<r;j++) operator_count[InvPermutationIndex[Sq_subseq[j]]]++;
	else    for(j=0;j<r;j++) operator_count[Sq_subseq[j]]++;
	permutation_count = factorial[r];
	for(i=0;i<Nop;i++) if(operator_count[i]>1) permutation_count /= factorial[operator_count[i]];
	found_cycles = 0;                                  // found_cycles contains the number of cycles found. it is global variable.
	for(i=0;i<Ncycles;i++){
		if(cycle_len[i]<lmin || cycle_len[i]>lmax) continue;
		not_contained = 0;
		if(inv){
			for(j=0;j<r;j++) if(cycles[i][InvPermutationIndex[Sq_subseq[j]]] < operator_count[InvPermutationIndex[Sq_subseq[j]]]){ not_contained = 1; break;}
		} else{
			for(j=0;j<r;j++) if(cycles[i][Sq_subseq[j]] < operator_count[Sq_subseq[j]]){ not_contained = 1; break;}
		}
		if(not_contained) continue;
		found_cycle_list[found_cycles++] = i;
	}
	return found_cycles>0 ? found_cycle_list[int(val(rng)*found_cycles)] : -1; // returns one of the found cycles chosen randomly.
}

int is_inverse(int k, int l){ // checks whether permutation # k and permutation # l are mutually inverse permutations
	int i, size1, size2, r = 1; std::vector<int> P1(N,0); std::vector<int> P2(N,0);
	size1 = CDPdata.Permutations[k].size(); size2 = CDPdata.Permutations[l].size();
	if(size1 == size2){
		for(i=0;i<size1;i++) P1[CDPdata.Permutations[k][i].first-1] = CDPdata.Permutations[k][i].second;
		for(i=0;i<size1;i++) P2[CDPdata.Permutations[l][i].first-1] = CDPdata.Permutations[l][i].second;
		for(i=0;i<N;i++) if((P1[i] + P2[i]) % twoSplus1 != 0){ r = 0; break; }
	} else r = 0;
	return r;
}

void init_rng(){
#ifdef EXACTLY_REPRODUCIBLE
	rng_seed = 1000 + mpi_rank;
#else
	rng_seed = rd();
#endif
	rng.seed(rng_seed);
	diceN.param(std::uniform_int_distribution<int>::param_type(0,N-1));
	diceNop.param(std::uniform_int_distribution<int>::param_type(0,Nop-1));
	diceS.param(std::uniform_int_distribution<int>::param_type(0,twoSplus1-1));
	diceSminus1.param(std::uniform_int_distribution<int>::param_type(0,twoSplus1-2));
}

void init_basic(){
	double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	for(int k=1;k<qmax;k++){ curr1*=(-double(beta))/k; curr2*=k; beta_pow_factorial[k] = curr1; factorial[k] = curr2;}
	currWeight = GetWeight();
}

void init(){
	operator_count.resize(Nop);
	int i,j,k,l; double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	for(q=1;q<qmax;q++){ curr1*=(-double(beta))/q; curr2*=q; beta_pow_factorial[q] = curr1; factorial[q] = curr2;}
	lattice.resize(N); for(i=N-1;i>=0;i--) lattice[i] = diceS(rng); z = lattice; q = 0;
	currWeight = GetWeight();
	for(i=0;i<Nop;i++){
		cycles.push_back(vector<int>(Nop,0));
		cycles[Ncycles++][i] = twoSplus1;
	}
	for(i=0;i<Ncycles;i++) cycle_len.push_back(vector_sum(cycles[i]));
	cycle_min_len = 64; for(i=0;i<Ncycles;i++) cycle_min_len = min(cycle_min_len,cycle_len[i]);
	cycle_max_len = 0; for(i=0;i<Ncycles;i++) cycle_max_len = max(cycle_max_len,cycle_len[i]);
#ifdef WORM_UPDATE
	cycles_used.resize(Ncycles,1);
#else
	cycles_used.resize(Ncycles,0);
#endif
	InvPermutationIndex.resize(Nop,-1); MainPermutationIndex.resize(Nop,0);
	for(k=0;k<Nop;k++) for(l=k;l<Nop;l++) if(is_inverse(k,l)){
		InvPermutationIndex[k] = l; InvPermutationIndex[l] = k;
		MainPermutationIndex[k] = 1; MainPermutationIndex[l] = 0;
	}
	for(i=0;i<Nop;i++) if(InvPermutationIndex[i] < 0){
		std::cout << "Error: some of the inverse permutations are not found" << std::endl; exit(1);
	} else if(InvPermutationIndex[InvPermutationIndex[i]] != i){
		std::cout << "Error: inverse(inverse(P_"<<i<<")) is not equal to " << "P_"<<i<< std::endl; exit(1);
	}
	for(i=0;i<N_all_observables;i++) in_bin_sum[i] = 0; in_bin_sum_sgn = 0;
	for(i=0;i<N_all_observables;i++) valid_observable[i] = 0;
#ifdef MEASURE_CUSTOM_OBSERVABLES
	for(i=0;i<Nobservables;i++) valid_observable[i] = 1;
#endif
#ifdef MEASURE_H
	valid_observable[Nobservables] = 1;
#endif
#ifdef MEASURE_H2
	valid_observable[Nobservables + 1] = 1;
#endif
#ifdef MEASURE_HDIAG
	valid_observable[Nobservables + 2] = 1;
#endif
#ifdef MEASURE_HDIAG2
	valid_observable[Nobservables + 3] = 1;
#endif
#ifdef MEASURE_HOFFDIAG
	valid_observable[Nobservables + 4] = 1;
#endif
#ifdef MEASURE_HOFFDIAG2
	valid_observable[Nobservables + 5] = 1;
#endif
#ifdef MEASURE_Z_MAGNETIZATION
	valid_observable[Nobservables + 6] = 1;
#endif
#ifdef MEASURE_Z_MAGNETIZATION2
	valid_observable[Nobservables + 7] = 1;
#endif
}

double Metropolis(ExExFloat newWeight){
	return min(1.0,fabs((newWeight/currWeight).get_double()));
}

#ifdef WORM_UPDATE

int worm_exit_check(){
	int r = save_data_flag;
#ifdef MAX_WORM_LENGTH
	if(worm_length > MAX_WORM_LENGTH) r = 1;
#endif
#ifdef WORM_BREAK_PROBABILITY
	if(val(rng) < WORM_BREAK_PROBABILITY) r = 1;
#endif
	return r;
}

void worm_common_part(double v){ // assuming v >= 0.6
	int i,m,p,r,u; double oldE; ExExFloat newWeight;
	if(v < 0.7){ // attempt to insert a randomly chosen single permutation operator
		if(q+1<qmax){
			saved_distance_from_identity = distance_from_identity;
			m = int(val(rng)*(q+1)); // m is between 0 and q
			memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
			old_currD = currD; p = diceNop(rng);
			for(i=q-1;i>=m;i--) Sq[i+1] = Sq[i]; q++; Sq[m] = p; newWeight = UpdateWeight();
			if(saved_distance_from_identity && distance_from_identity && worm_exit_check()){
				q=q_saved; memcpy(Sq,Sq_saved,q_saved*sizeof(int));
				currWeight = UpdateWeight();
			} else if(val(rng) < Metropolis(newWeight*Nop)){
				currWeight = newWeight;
				if(saved_distance_from_identity == 0){
					memcpy(Sq_saved,Sq_backup,(q-1)*sizeof(int)); q_saved = q-1;
				}				
			} else{
				q--; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
				currD = old_currD; UpdateWeight();
			}
		} else qmax_achieved = 1;
	} else if(v < 0.8){ // attempt to remove a randomly chosen single permutation operator
		if(q>0){
			saved_distance_from_identity = distance_from_identity;
			m = int(val(rng)*q); p = Sq[m]; // m is between 0 and (q-1)
			old_currD = currD; memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
			for(i=m;i<q-1;i++) Sq[i] = Sq[i+1]; q--; newWeight = UpdateWeight();
			if(saved_distance_from_identity && distance_from_identity && worm_exit_check()){
				q=q_saved; memcpy(Sq,Sq_saved,q_saved*sizeof(int));
				currWeight = UpdateWeight();
			} else if(val(rng) < Metropolis(newWeight/Nop)){
				currWeight = newWeight;
				if(saved_distance_from_identity == 0){
					memcpy(Sq_saved,Sq_backup,(q+1)*sizeof(int)); q_saved = q+1;
				}
			} else{
				q++; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
				currD = old_currD; UpdateWeight();
			}
		}
	} else if(v < 0.9){ // attempting a block swap
		if(q>=2 && distance_from_identity==0){
			m = q==2 ? 0 : int(val(rng)*(q-1)); // m is between 0 and (q-2)
			lattice_backup = lattice;
			oldE = currEnergy; for(i=0;i<=m;i++) ApplyOperator(Sq[i]);
			for(i=0;i<=m;i++)  Sq_backup[q-1-m+i] = Sq[i];
			for(i=m+1;i<q;i++) Sq_backup[i-m-1] = Sq[i];
			for(i=0;i<q;i++) { p = Sq[i]; Sq[i] = Sq_backup[i]; Sq_backup[i] = p;}
			memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
			GetEnergies(); newWeight = UpdateWeightReplace(oldE,currEnergy);
			if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else{
				UpdateWeightReplace(currEnergy,oldE); currEnergy = oldE;
				lattice = lattice_backup; memcpy(Sq,Sq_backup,q*sizeof(int));
				memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
			}
		}
	} else if(distance_from_identity==0){ // flip of a random spin
		p = diceN(rng); u = lattice[p]; lattice_flip(p); newWeight = UpdateWeight();
		if(val(rng) < Metropolis(newWeight)) { z = lattice; currWeight = newWeight;}
			else { lattice[p] = u; currWeight = UpdateWeight();}
	}
}

#ifdef COMPOSITE_INSIDE_WORM

void worm_update(){  // employing worm update in addition to other updates + employing composite update inside worm update
	int i,m,p,r,u,oldq,cont; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight; double Rfactor;
	if(v < 0.6){  // composite update
		Rfactor = 1; oldq = q; memcpy(Sq_backup,Sq,q*sizeof(int));
		cycles_used_backup = cycles_used; newWeight = currWeight;
		do{
			cont = 0; v = val(rng);
			if(v < 0.25){ // attempt to swap Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(Sq[m]!=Sq[m+1]){
						oldE = Energies[m+1];
						p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p;
						GetEnergies();
						newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.5){ // attempt to delete Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(InvPermutationIndex[Sq[m]]==Sq[m+1]){
						oldE = Energies[m]; oldE2 = Energies[m+1];
						for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
						GetEnergies(); Rfactor /= Nop;
						newWeight = UpdateWeightDel(oldE,oldE2);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.75){ // attempt to insert Sq[m] and Sq[m+1]
				if(q+2<qmax){
					m = int(val(rng)*(q+1)); // m is between 0 and q
					p = diceNop(rng);
					for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = p; Sq[m+1] = InvPermutationIndex[p];
					GetEnergies(); Rfactor *= Nop;
					newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
					if(REALW(currD) == 0) cont = 1;
				} else qmax_achieved = 1;
			} else{ // attempting a fundamental cycle completion
				int j = 0, s, t, inv, inv_pr; double wfactor;
				u = geometric_int(rng); // a random integer u is picked according to geometric distribution
				if(q >= u+rmin){
					inv = dice2(rng); inv_pr = min(rmax,q-u)-(rmin)+1;
					r = int(val(rng)*inv_pr) + (rmin);  // r is random integer between rmin and min(rmax,q-u)
					PickSubsequence(r+u); // indexes of the subsequence are min_index, min_index+1,..., max_index=min_index+r+u-1
					std::shuffle(Sq_subseq,Sq_subseq+r+u,rng);
					for(i=0;i<u;i++) Sq_gaps[i] = Sq_subseq[i+r];
					m = FindCycles(r,inv);
					if(found_cycles > 0){ // cycles[m] is one of the found cycles, containing all the operators of Sq_subseq
						P = cycles[m]; for(i=0;i<Nop;i++) if(P[i]!=0) P[i] -= operator_count[i]; // P is remaining part of the cycle
						if(inv==0){ // if inv = 0, then P is inverse of the remaining part of the cycle
							for(i=0;i<Nop;i++) if(MainPermutationIndex[i]){
								t = P[i];
								P[i] = P[InvPermutationIndex[i]];
								P[InvPermutationIndex[i]] = t;
							}
						}
						p = vector_sum(P); // here, p is length of the sequence S' to which S will be replaced
						inv = 1-inv;
						if(q+p-r >= qmax) qmax_achieved = 1; else if(p >= rmin && p <= rmax){
							if(r<p)	     for(i=q-1;i>max_index;i--) Sq[i+p-r] = Sq[i]; // shift the values to the right
							else if(r>p) for(i=max_index+1;i<q;i++) Sq[i+p-r] = Sq[i]; // shift the values to the left
							i=0; while(i<p){ 
								while(P[j]==0) j++; t = i + P[j];
								for(s=i; s<t; s++) Sq_subseq[s] = Sq[min_index+s] = j; j++; i = t;
							}
							for(i=0;i<u;i++) Sq[min_index+p+i] = Sq_gaps[i]; // S' contains the remaining operators
							std::shuffle(Sq+min_index,Sq+min_index+p+u,rng);
							q += p-r; // the length q may have changed
							newWeight = UpdateWeight();
							wfactor = found_cycles/permutation_count; FindCycles(p,inv); wfactor /= found_cycles/permutation_count;
							wfactor *= inv_pr; inv_pr = min(rmax,q-u)-(rmin)+1; wfactor /= inv_pr;
							Rfactor *= wfactor;
							cycles_used[m] = 1;
							if(REALW(currD) == 0) cont = 1;
						}
					}
				}
			}
		} while(cont && !save_data_flag 
#ifdef COMPOSITE_UPDATE_BREAK_PROBABILITY
			&& val(rng) > COMPOSITE_UPDATE_BREAK_PROBABILITY
#endif
			);
		if(REALW(currD) != 0 && val(rng) < Metropolis(newWeight*Rfactor)){
			currWeight = newWeight;
		} else{
			q = oldq;
			memcpy(Sq,Sq_backup,q*sizeof(int));
			cycles_used = cycles_used_backup;
			currWeight = UpdateWeight();
		}
	} else worm_common_part(v);
}

#else

void worm_update(){  // employing worm update in addition to other updates
	int i,m,p,r,u; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight;
	if(v < 0.3){ //local move
		if(v<0.1){ // attempt to swap Sq[m] and Sq[m+1]
			if(q>=2){
				m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
				if(Sq[m]!=Sq[m+1]){
					oldE = Energies[m+1]; old_currD = currD;
					p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p; GetEnergies(); newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
					if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else {
						Sq[m+1] = Sq[m]; Sq[m] = p; currD = old_currD;
						UpdateWeightReplace(Energies[m+1],oldE); Energies[m+1] = oldE;
					}
				}
			}
		} else if(v<0.2){ // attempt to delete Sq[m] and Sq[m+1]
			if(q>=2){
				m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
				if(InvPermutationIndex[Sq[m]]==Sq[m+1]){
					oldE = Energies[m]; oldE2 = Energies[m+1]; old_currD = currD;
					memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
					for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
					GetEnergies(); newWeight = UpdateWeightDel(oldE,oldE2);
					if(val(rng) < Metropolis(newWeight/Nop)) currWeight = newWeight; else{
						q+=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
						currD = old_currD; UpdateWeightIns(oldE,oldE2);
					}
				}
			}
		} else{
			if(q+2<qmax){ // attempt to insert Sq[m] and Sq[m+1]
				m = int(val(rng)*(q+1)); // m is between 0 and q
				memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
				old_currD = currD; p = diceNop(rng);
				for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = p; Sq[m+1] = InvPermutationIndex[p];
				GetEnergies(); newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
				if(val(rng) < Metropolis(newWeight*Nop)) currWeight = newWeight; else{
					q-=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
					currD = old_currD; d->RemoveElement(); d->RemoveElement();
				}
			} else qmax_achieved = 1;
		}
	} else if(v < 0.6){ // attempting a fundamental cycle completion
		int oldq, j = 0, s, t, inv, inv_pr; double wfactor;
		u = geometric_int(rng); // a random integer u is picked according to geometric distribution
		if(q >= u+rmin){
			inv = dice2(rng); inv_pr = min(rmax,q-u)-(rmin)+1;
			r = int(val(rng)*inv_pr) + (rmin);  // r is random integer between rmin and min(rmax,q-u)
			PickSubsequence(r+u); // indexes of the subsequence are min_index, min_index+1,..., max_index=min_index+r+u-1
			std::shuffle(Sq_subseq,Sq_subseq+r+u,rng);
			for(i=0;i<u;i++) Sq_gaps[i] = Sq_subseq[i+r];
			m = FindCycles(r,inv);
			if(found_cycles > 0){ // cycles[m] is one of the found cycles, containing all the operators of Sq_subseq
				P = cycles[m]; for(i=0;i<Nop;i++) if(P[i]!=0) P[i] -= operator_count[i]; // P is remaining part of the cycle
				if(inv==0){ // if inv = 0, then P is inverse of the remaining part of the cycle
					for(i=0;i<Nop;i++) if(MainPermutationIndex[i]){
						t = P[i];
						P[i] = P[InvPermutationIndex[i]];
						P[InvPermutationIndex[i]] = t;
					}
				}
				p = vector_sum(P); // here, p is length of the sequence S' to which S will be replaced
				inv = 1-inv;
				if(q+p-r >= qmax) qmax_achieved = 1; else if(p >= rmin && p <= rmax){
					memcpy(Sq_backup,Sq,q*sizeof(int)); oldq = q;
					if(r<p)	     for(i=q-1;i>max_index;i--) Sq[i+p-r] = Sq[i]; // shift the values to the right
					else if(r>p) for(i=max_index+1;i<q;i++) Sq[i+p-r] = Sq[i]; // shift the values to the left
					i=0; while(i<p){ 
						while(P[j]==0) j++; t = i + P[j];
						for(s=i; s<t; s++) Sq_subseq[s] = Sq[min_index+s] = j; j++; i = t;
					}
					for(i=0;i<u;i++) Sq[min_index+p+i] = Sq_gaps[i]; // S' contains the remaining operators
					std::shuffle(Sq+min_index,Sq+min_index+p+u,rng);
					q += p-r; // the length q may have changed
					newWeight = UpdateWeight();
					wfactor = found_cycles/permutation_count; FindCycles(p,inv); wfactor /= found_cycles/permutation_count;
					wfactor *= inv_pr; inv_pr = min(rmax,q-u)-(rmin)+1; wfactor /= inv_pr;
					if(val(rng) < Metropolis(newWeight*wfactor)){
						currWeight = newWeight; cycles_used[m] = 1;
					} else{
						q = oldq; memcpy(Sq,Sq_backup,q*sizeof(int)); currWeight = UpdateWeight();
					}
				}
			}
		}
	} else worm_common_part(v);
}

#endif

void update(){
	if(save_data_flag){ save_QMC_data(); exit(0); }
	worm_length = 0; do{ worm_update(); worm_length++; } while(distance_from_identity);
}

#else

void update(){  // employing composite update rather than worm update
	if(save_data_flag){ save_QMC_data(); exit(0); };
	int i,m,p,r,u,oldq,cont; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight; double Rfactor;
	if(v < 0.8){
		Rfactor = 1; oldq = q; memcpy(Sq_backup,Sq,q*sizeof(int));
		cycles_used_backup = cycles_used; newWeight = currWeight;
		do{
			cont = 0; v = val(rng);
			if(v < 0.25){ // attempt to swap Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(Sq[m]!=Sq[m+1]){
						oldE = Energies[m+1];
						p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p;
						GetEnergies();
						newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.5){ // attempt to delete Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(InvPermutationIndex[Sq[m]]==Sq[m+1]){
						oldE = Energies[m]; oldE2 = Energies[m+1];
						for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
						GetEnergies(); Rfactor /= Nop;
						newWeight = UpdateWeightDel(oldE,oldE2);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.75){ // attempt to insert Sq[m] and Sq[m+1]
				if(q+2<qmax){
					m = int(val(rng)*(q+1)); // m is between 0 and q
					p = diceNop(rng);
					for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = p; Sq[m+1] = InvPermutationIndex[p];
					GetEnergies(); Rfactor *= Nop;
					newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
					if(REALW(currD) == 0) cont = 1;
				} else qmax_achieved = 1;
			} else{ // attempting a fundamental cycle completion
				int j = 0, s, t, inv, inv_pr; double wfactor;
				u = geometric_int(rng); // a random integer u is picked according to geometric distribution
				if(q >= u+rmin){
					inv = dice2(rng); inv_pr = min(rmax,q-u)-(rmin)+1;
					r = int(val(rng)*inv_pr) + (rmin);  // r is random integer between rmin and min(rmax,q-u)
					PickSubsequence(r+u); // indexes of the subsequence are min_index, min_index+1,..., max_index=min_index+r+u-1
					std::shuffle(Sq_subseq,Sq_subseq+r+u,rng);
					for(i=0;i<u;i++) Sq_gaps[i] = Sq_subseq[i+r];
					m = FindCycles(r,inv);
					if(found_cycles > 0){ // cycles[m] is one of the found cycles, containing all the operators of Sq_subseq
						P = cycles[m]; for(i=0;i<Nop;i++) if(P[i]!=0) P[i] -= operator_count[i]; // P is remaining part of the cycle
						if(inv==0){ // if inv = 0, then P is inverse of the remaining part of the cycle
							for(i=0;i<Nop;i++) if(MainPermutationIndex[i]){
								t = P[i];
								P[i] = P[InvPermutationIndex[i]];
								P[InvPermutationIndex[i]] = t;
							}
						}
						p = vector_sum(P); // here, p is length of the sequence S' to which S will be replaced
						inv = 1-inv;
						if(q+p-r >= qmax) qmax_achieved = 1; else if(p >= rmin && p <= rmax){
							if(r<p)	     for(i=q-1;i>max_index;i--) Sq[i+p-r] = Sq[i]; // shift the values to the right
							else if(r>p) for(i=max_index+1;i<q;i++) Sq[i+p-r] = Sq[i]; // shift the values to the left
							i=0; while(i<p){ 
								while(P[j]==0) j++; t = i + P[j];
								for(s=i; s<t; s++) Sq_subseq[s] = Sq[min_index+s] = j; j++; i = t;
							}
							for(i=0;i<u;i++) Sq[min_index+p+i] = Sq_gaps[i]; // S' contains the remaining operators
							std::shuffle(Sq+min_index,Sq+min_index+p+u,rng);
							q += p-r; // the length q may have changed
							newWeight = UpdateWeight();
							wfactor = found_cycles/permutation_count; FindCycles(p,inv); wfactor /= found_cycles/permutation_count;
							wfactor *= inv_pr; inv_pr = min(rmax,q-u)-(rmin)+1; wfactor /= inv_pr;
							Rfactor *= wfactor;
							cycles_used[m] = 1;
							if(REALW(currD) == 0) cont = 1;
						}
					}
				}
			}
		} while(cont && !save_data_flag
#ifdef COMPOSITE_UPDATE_BREAK_PROBABILITY
			&& val(rng) > COMPOSITE_UPDATE_BREAK_PROBABILITY
#endif
			);
		if(REALW(currD) != 0 && val(rng) < Metropolis(newWeight*Rfactor)){
			currWeight = newWeight;
		} else{
			q = oldq;
			memcpy(Sq,Sq_backup,q*sizeof(int));
			cycles_used = cycles_used_backup;
			currWeight = UpdateWeight();
		}
	} else if(v < 0.9 && q>=2){ // attempting a block swap
		m = q==2 ? 0 : int(val(rng)*(q-1)); // m is between 0 and (q-2)
		oldE = currEnergy; for(i=0;i<=m;i++) ApplyOperator(Sq[i]);
		for(i=0;i<=m;i++)  Sq_backup[q-1-m+i] = Sq[i];
		for(i=m+1;i<q;i++) Sq_backup[i-m-1] = Sq[i];
		for(i=0;i<q;i++) { p = Sq[i]; Sq[i] = Sq_backup[i]; Sq_backup[i] = p;}
		memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
		GetEnergies(); newWeight = UpdateWeightReplace(oldE,currEnergy);
		if(val(rng) < Metropolis(newWeight)){
			z = lattice; currWeight = newWeight;
		} else{ UpdateWeightReplace(currEnergy,oldE); currEnergy = oldE;
			lattice = z; memcpy(Sq,Sq_backup,q*sizeof(int));
			memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
		}
	} else{ // flip of a random spin
		p = diceN(rng); u = lattice[p]; lattice_flip(p); newWeight = UpdateWeight();
		if(val(rng) < Metropolis(newWeight)) { z = lattice; currWeight = newWeight;}
			else { lattice[p] = u; currWeight = UpdateWeight();}
	}
}

#endif

#ifdef MEASURE_CUSTOM_OBSERVABLES

std::complex<double> calc_MD0(int n){ // calculate <z | MD_0 | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD0_size[n];i++) sum -= double(2*(int((MD0_product[n][i] & (~lattice)).count())%2)-1) * MD0_coeff[n][i];
	return sum;
}

std::complex<double> calc_MD(int n, int k){ // calculate d_k = <z | MD_k | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD_size[n][k];i++) sum -= double(2*(int((MD_product[n][k][i] & (~lattice)).count())%2)-1) * MD_coeff[n][k][i];
	return sum;
}

#endif

double measure_H(){
	double R = d->z[q]/(-beta);
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_H2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta));
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Hdiag(){
	return currEnergy;
}

double measure_Hdiag2(){
	return currEnergy*currEnergy;
}

double measure_Hoffdiag(){
	double R = 0;
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_Hoffdiag2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta)) + currEnergy*(currEnergy - 2*measure_H());
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Z_magnetization(){
	int sum = 0; double m;
	for(int i = 0; i < lattice.size(); i++) sum += lattice[i];
	m = double(sum)/(Svalue*(2*N)); // this value is between 0 and 1
	return m;
}

double measure_Z_magnetization2(){
	double m = measure_Z_magnetization();
	return m*m;
}

std::string name_of_observable(int n){
	std::string s;
	if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
		s = Mnames[n];
#endif
	} else switch(n-Nobservables){
			case 0: s = "H";             break;
			case 1: s = "H^2";           break;
			case 2: s = "H_{diag}";      break;
			case 3: s = "H_{diag}^2";    break;
			case 4: s = "H_{offdiag}";   break;
			case 5: s = "H_{offdiag}^2"; break;
			case 6: s = "Z_magnetization"; break;
			case 7: s = "Z_magnetization^2"; break;
	}
	return s;
}

double measure_observable(int n){
	double R = 0;
	if(valid_observable[n]){
		if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
			int i,k,len,cont;
			std::complex<double> T = calc_MD0(n);
			for(k=0;k<MNop[n];k++){
				P = MP[n][k]; len = P.count(); if(len>q) continue;
				if(!NoRepetitionCheck(Sq+(q-len),len)) continue;
				cont = 0; for(i=0;i<len;i++) if(!P.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
				T +=	(d->divdiffs[q-len]/d->divdiffs[q]).get_double() *
				        (beta_pow_factorial[q-len]/beta_pow_factorial[q]).get_double()/factorial[len] *
					(currD_partial[q-len]/currD) * calc_MD(n,k);
			}
			if(abs_weights)
				R = std::abs(T)*currWeight.sgn()*cos(std::arg(T)+std::arg(currD));
			else 	R = std::real(currD*T)/std::real(currD); // we importance-sample Re(W_C A_C)/Re(W_C)
#endif
		} else  switch(n-Nobservables){
				case 0:	R = measure_H(); break;
				case 1:	R = measure_H2(); break;
				case 2:	R = measure_Hdiag(); break;
				case 3:	R = measure_Hdiag2(); break;
				case 4:	R = measure_Hoffdiag(); break;
				case 5:	R = measure_Hoffdiag2(); break;
				case 6: R = measure_Z_magnetization(); break;
				case 7: R = measure_Z_magnetization2(); break;
		}
	}
	return R;
}

void measure(){
	double R, sgn; int i;
	currWeight = GetWeight();
	if(abs_weights)	sgn = currWeight.sgn() * cos(std::arg(currD)); // arg(W) = arg(currD) + arg(currWeight), where arg(currWeight) = either 0 or Pi
		else 	sgn = currWeight.sgn();
	meanq += q; if(maxq < q) maxq = q; in_bin_sum_sgn += sgn;
	if((measurement_step+1) % bin_length == 0){
		in_bin_sum_sgn /= bin_length; bin_mean_sgn[measurement_step/bin_length] = in_bin_sum_sgn; in_bin_sum_sgn = 0;
	}
	for(i=0;i<N_all_observables;i++){
		R = measure_observable(i); in_bin_sum[i] += R*sgn;
		if((measurement_step+1) % bin_length == 0){
			in_bin_sum[i] /= bin_length; bin_mean[i][measurement_step/bin_length] = in_bin_sum[i]; in_bin_sum[i] = 0;
		}
	}
}

double mean_O[N_all_observables], stdev_O[N_all_observables], mean_O_backup[N_all_observables];

const int N_derived_observables = 2;  // we define number of derived observables

std::string name_of_derived_observable(int n){ // we define names of derived observables
	std::string s;
	switch(n){
		case 0 : s = "specific heat"; break;
		case 1 : s = "magnetic susceptibility"; break;
	}
	return s;
}

int valid_derived_observable(int n){ // we define which observables are needed for each derived observable
	int r = 0;
	switch(n){
		case 0: r = valid_observable[Nobservables+0] && valid_observable[Nobservables+1]; break;
		case 1: r = valid_observable[Nobservables+6] && valid_observable[Nobservables+7]; break;
	}
	return r;
}

double compute_derived_observable(int n){ // we compute the derived observables
	double R = 0;
	switch(n){
		case 0 : R = beta*beta*(mean_O[Nobservables+1] - mean_O[Nobservables+0]*mean_O[Nobservables+0]); break;
		case 1 : R = beta*(mean_O[Nobservables+7] - mean_O[Nobservables+6]*mean_O[Nobservables+6]); break;
	}
	return R;
}
