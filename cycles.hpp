vector<int> good_cycle;

void PrintList(ExExFloat* list, int len){
	int i;
	printf("{");
	for(i=0;i<len;i++){
		list[i].print();
		if(i<len-1) printf(",");
	}
	printf("}");
}

void PrintList(double* list, int len){
	int i;
	printf("{");
	for(i=0;i<len;i++){
		printf("%.17g",list[i]);
		if(i<len-1) printf(",");
	}
	printf("}");
}

void PrintList(int* list, int len){
	int i;
	printf("{");
	for(i=0;i<len;i++){
		printf("%d",list[i]);
		if(i<len-1) printf(",");
	}
	printf("}");
}

void PrintList(std::vector<int> list){
	int i; int len = list.size();
	printf("{");
	for(i=0;i<len;i++){
		printf("%d",list[i]);
		if(i<len-1) printf(",");
	}
	printf("}");
}

vector<int> modulo_add(vector<int> vec1 , vector<int> vec2){
    vector<int> result;
    if (vec1.size() != vec2.size()){ cout << "The size of the added vectors are different!" << endl; exit(1);}
    for(int i = 0; i < vec1.size(); i++) result.push_back((vec1[i] + vec2[i])%twoSplus1);
    return result;
}

int Int_sum(vector<int> vec){
    int sum = 0;
    // for(int i = 0; i < vec.size(); i++) sum += vec[i];
    for(int i = 0; i < vec.size(); i++) if(vec[i] != 0) sum++;
    return sum;
}

bool compare_null(const vector<int> & a, const vector<int> & b){
	return Int_sum(a) < Int_sum(b);
}

int test_cycle(vector<int> cycle){
	int i,j,k,l,m, passed = 0; q = 0;
	for(i=0;i<Nop;i++) for(m=0;m<cycle[i];m++) Sq[q++] = i;
	for(k=0;k<N_CYCLE_TESTS;k++){
		shuffle(Sq,Sq+q,rng);
		for(i=0;i<N;i++) lattice[i] = diceS(rng);
		GetEnergies();
		if(currD.real() != 0){
			passed = 1;
			break;
		}
	}
	return passed;
}

int cycle_minimize(vector<vector<int>>& null_eigs, int test_cycles){
	int nullsize, k, m, null_k, changes_made = 0; vector<int> curr;
	nullsize = null_eigs.size();
	sort(null_eigs.begin(), null_eigs.end(), compare_null);
	for(k=nullsize-1; k>0; k--){
		null_k = Int_sum(null_eigs[k]);
		for(m=0; m<k; m++){
			curr = modulo_add( null_eigs[k] , null_eigs[m]);
			if(Int_sum(curr) < null_k) if(!test_cycles || test_cycle(curr)){
				null_eigs[k] = curr; changes_made = 1;
				break;
			}
		}
	}
	return changes_made;
}

int fix_cycle(int k){
	vector<int> curr; int i, j, success = 0;
	for(i=0;i<Ncycles;i++){
		if(i==k) continue;
		curr = cycles[k];
		for(j=0;j<twoSplus1;j++){
			if(test_cycle(curr)){
				success = 1; cycles[k] = curr; good_cycle[k] = 1; break;
			} else curr = modulo_add(curr, cycles[i]);
		}
		if(success) break;
	}
	return success;
}

int fix_cycle(int k, int m){
	vector<int> curr; int i, j, p, success = 0;
	vector<int> involved(Ncycles,0);
	for(i=0;i<N_CYCLE_FIX_ATTEMPTS;i++){
		fill(involved.begin(),involved.end(),0);
		for(j=0;j<m;j++){
			p = int(val(rng)*(Ncycles-1)); if(p>=k) p++; // p is between 0 and Ncycles-1, p is not equal to k.
			involved[p]++;
		}
		curr = cycles[k];
		for(p=0;p<Ncycles;p++) for(j=0;j<involved[p];j++) curr = modulo_add(curr, cycles[p]);
		if(test_cycle(curr)){
			success = 1; cycles[k] = curr; good_cycle[k] = 1; break;
		}
	}
	return success;
}

int fix_cycles(){
	int i,k,status; lattice.resize(N); good_cycle.resize(Ncycles);
	while(cycle_minimize(cycles,false));
	for(i=0;i<Ncycles;i++) good_cycle[i] = test_cycle(cycles[i]);
	for(i=0;i<Ncycles;i++) if(!good_cycle[i]) fix_cycle(i); // attempting to correct i-th cycle
	for(k=2;k<5;k++) for(i=0;i<Ncycles;i++) if(!good_cycle[i])  fix_cycle(i,k); // attempting to correct i-th cycle
	while(cycle_minimize(cycles,true));
	for(i=0;i<Ncycles;i++) good_cycle[i] = test_cycle(cycles[i]);
	status = 1; for(i=0;i<Ncycles;i++) status &= good_cycle[i];
	return status;
}
