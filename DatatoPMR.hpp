#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bitset>

using namespace std;

typedef vector<vector<complex<double>>> Coeffs;
typedef vector<vector<int>> ParticlePerm;
typedef vector<vector<pair<int,int>>> TotalPerms;
typedef pair<int,int> ParticleDiag;
struct TotalDiag {
    int ztype , k , particle;
};
typedef vector<vector<vector<ParticleDiag>>> ParticleDVecs;
typedef vector<vector<vector<TotalDiag>>> TotalDVecs;
typedef pair<vector<int> , pair<vector<ParticleDiag> , vector<complex<double>>>> PauliCDPs; // This typedef is to summarize the Paulis as sum of DPs (with corresponding coefficients);
struct PDdata {
    int twoSplusOne , NumOfParticles;
    TotalPerms Permutations;
    TotalDVecs Diagonals;
    Coeffs Coefficients;
    vector<vector<TotalDiag>> D0;
    vector<complex<double>> D0Coeffs;
};

template<typename T>
void printMatrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    int n = matrix[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void Print_diagonals(vector<vector<TotalDiag>> Diags , vector<complex<double>> Cs){
    for(int j = 0; j < Diags.size(); j++){
        cout <<  Cs[j];
        for(int k = 0; k < Diags[j].size() ; k++){
            if(Diags[j][k].ztype == 1){
                cout << " D(z," << Diags[j][k].k << "," << "spin " << Diags[j][k].particle << ")";
            }
            else{
                cout << " D(" << Diags[j][k].k << "," << "spin " << Diags[j][k].particle << ")";
            }
        }
        if(j < Diags.size()-1){
            cout << " + ";
        }
        else{
            cout << endl;
        }
    }
    cout << endl;
}

void Print_data(PDdata CDPdata){
    Coeffs C = CDPdata.Coefficients;
    TotalDVecs D = CDPdata.Diagonals;
    TotalPerms P = CDPdata.Permutations;
    vector<vector<TotalDiag>> D0 = CDPdata.D0;
    vector<complex<double>> D0Coeffs = CDPdata.D0Coeffs;

    int D0size = D0.size();
    if(D0size > 0){
        cout << "The diagonal D0 is " << endl;
        for(int i = 0; i < D0size; i++){
            cout << D0Coeffs[i];
            for(int j = 0; j < D0[i].size(); j++){
                if(D0[i][j].ztype == 1){
                    cout <<" D(z," << D0[i][j].k << "," << "spin # " << D0[i][j].particle <<")";
                }
                else{
                    cout << " D(" << D0[i][j].k << "," << "spin # " << D0[i][j].particle << ")";
                }
            }
            if(i < D0size - 1){
                cout << " + ";
            }
        }
        cout << endl;
    }
    else{
        cout << "There are no purely diagonal terms.." << endl;
    }
    cout << endl;

    for(int i = 0; i < P.size() ; i++){
        cout << "The permutation is ";
        for(int l = 0; l < P[i].size(); l++){
            cout << "< spin # " << P[i][l].first << " , P^" << P[i][l].second << " > ";
        }
        cout << endl;
        for(int j = 0; j < D[i].size(); j++){
            cout <<  C[i][j];
            for(int l = 0; l < D[i][j].size() ; l++){
                if(D[i][j][l].ztype == 1){
                    cout << " D(z," << D[i][j][l].k << "," << "spin # " << D[i][j][l].particle <<")";
                }
                else{
                    cout << " D(" << D[i][j][l].k << "," << "spin # " << D[i][j][l].particle << ")";
                }
            }
            if(j < D[i].size()-1){
                cout << " + ";
            }
            else{
                cout << endl;
            }
        }
        cout << endl;
    }
}

vector<pair<complex<double>, vector<int>>> data_extract(const string& fileName){
    vector<pair<complex<double>, vector<int>>> data;

    ifstream inputFile(fileName);
    if (!inputFile) {
        cout << "Failed to open the input file!" << endl; exit(1);
    }

    string line;
    while (getline(inputFile, line)) {
        // The first non-empty element is the coefficient
        istringstream iss(line);
        pair<complex<double>,vector<int>> linedata;

        // Extracting the complex coefficient:
        double realpart, imagpart=0;
        char sign;
        string complexPart;
        iss >> complexPart; 
	if(complexPart.length() == 0) continue; // empty lines, as well as lines containing nothing but spaces, are ignored

        istringstream complexIss(complexPart);
        complexIss >> realpart >> imagpart;
        linedata.first = complex<double> (realpart , imagpart);

	if(linedata.first == 0.0) continue; // if the coefficient is equal to zero, the line is ignored

        //Extracting the integer vectors of qubits and paulis:
        string token;
        vector<int> integers;
        while (iss >> token){
            // int num = std::stoi(token);
            int num = token == "X" || token == "x" ? 1 : 
                      token == "Y" || token == "y" ? 2 : 
                      token == "Z" || token == "z" ? 3 : std::stoi(token);
            integers.push_back(num);
        }
        if(integers.size()%4 != 0){
              cout << endl << "Error: the number of values in an input file line is incorrect" << endl; exit(1);
        }
        linedata.second = integers;

        data.push_back(linedata);
    }
    inputFile.close();
    if(data.size() == 0){
         cout << endl << "Error: input file is empty" << endl; exit(1);
    }
    return data;
}

// This function finds an instance of an integer in a vector of integers!
pair<bool , int> Find_number(int Num , vector<int> NumVec){
    for(int i = 0; i < NumVec.size() ; i++){
        if(NumVec[i] == Num){
            return {true , i};
        }
    }
    return {false , 0};
}

bool PartDiag_equal(ParticleDiag Da , ParticleDiag Db){
    if(Da == Db)
        return true;
    else
        return false;
}

bool PartDiag_compare(vector<ParticleDiag> DiagA , vector<ParticleDiag> DiagB){
    vector<ParticleDiag> DiagB2;
    bool DiagAiMatched;
    if(DiagA.size() != DiagB.size()){
        return false;
    }
    for(int i = 0; i < DiagA.size(); i++){
        DiagAiMatched = false;
        //DiagB2 = DiagB;
        for(int j = 0; j < DiagB2.size(); j++){
            if(PartDiag_equal(DiagA[i] , DiagB2[j])){
                DiagB2.erase(DiagB2.begin() + j);
                DiagAiMatched = true;
                break;
            }
        }
        if(!DiagAiMatched){
            return false;
        }
    }
    return true;
}

// Test out these Diag compares!

void PartDiag_append(vector<vector<ParticleDiag>>& AllDiagsi , vector<complex<double>>& AllCoesi , vector<vector<ParticleDiag>> NewDiags , vector<complex<double>> NewCoes){
    for(int i = 0; i < NewDiags.size() ; i++){
        bool NewDiagiFound = false;
        for(int j = 0; j < AllDiagsi.size(); j++){
            if(PartDiag_compare(NewDiags[i] , AllDiagsi[j])){
                AllCoesi[j] += NewCoes[i];
                NewDiagiFound = true;
                if(abs(AllCoesi[j]) < 1E-8){
                    AllDiagsi.erase(AllDiagsi.begin() + j);
                    AllCoesi.erase(AllCoesi.begin() + j);
                }
                break;
            }
        }
        if(!NewDiagiFound){
            AllDiagsi.push_back(NewDiags[i]);
            AllCoesi.push_back(NewCoes[i]);
        }
    }
}

void Permutation_append(vector<int>& AllPermsOnParticle , ParticleDVecs& AllDiagsOnParticle , Coeffs& AllCoeffsOnParticle , PauliCDPs CurrentPnDnCs , int twosp1){
    vector<int> CurrentPs = CurrentPnDnCs.first;
    vector<ParticleDiag> CurrentDs = CurrentPnDnCs.second.first;
    vector<complex<double>> CurrentCs = CurrentPnDnCs.second.second;

    // We have to go through AllPermsOnParticle and generate a new term by adding the c
    vector<int> AllPInit = AllPermsOnParticle;
    AllPermsOnParticle.clear();
    ParticleDVecs AllDInit = AllDiagsOnParticle;
    AllDiagsOnParticle.clear();
    Coeffs AllCoes = AllCoeffsOnParticle;
    AllCoeffsOnParticle.clear();

    // Initializing the diagonal element to be moved to the left of the existing permutations
    ParticleDiag CurrDj;
    CurrDj = {CurrentDs[0].first , 0};

    for(int j=0; j < CurrentPs.size(); j++){
        vector<int> AllP2 = AllPInit;
        ParticleDVecs AllDiags2 = AllDInit;
        Coeffs AllCoes2 = AllCoes;
        vector<int> RemovePermsIndices={};
        for(int i=0; i < AllP2.size(); i++){
            AllP2[i] = (AllP2[i] + CurrentPs[j]) % twosp1;
            CurrDj.second = (CurrentDs[j].second + AllPInit[i]) % twosp1; // Moving the new D matrix to the right of the existing permutation matrices!
            // Adding the new D matrix to the list of the existing D matrices!
            for(int k = 0; k < AllDiags2[i].size() ; k++){
                AllDiags2[i][k].push_back(CurrDj);
            }
            // If the new permutation already exists, only add the new diagonal to the list of existing ones
            pair<bool , int> PermFound = Find_number(AllP2[i] , AllPermsOnParticle); 
            if(PermFound.first){
                for(int k=0; k < AllCoes2[i].size(); k++){
                    AllCoes2[i][k] = AllCoes2[i][k]*CurrentCs[j];
                }
                PartDiag_append(AllDiagsOnParticle[PermFound.second] , AllCoeffsOnParticle[PermFound.second] , AllDiags2[i], AllCoes2[i]);
                // Remove the permutation, and the added diagonals with corresponding coefficients
                RemovePermsIndices.push_back(i);
            }
        }
        for(int k=0; k < RemovePermsIndices.size(); k++){
            AllP2.erase(AllP2.begin() + RemovePermsIndices[k] - k);
            AllDiags2.erase(AllDiags2.begin() + RemovePermsIndices[k] - k);
            AllCoes2.erase(AllCoes2.begin() + RemovePermsIndices[k] - k);
        }
        AllPermsOnParticle.insert(AllPermsOnParticle.end() , AllP2.begin() , AllP2.end());
        AllDiagsOnParticle.insert(AllDiagsOnParticle.end() , AllDiags2.begin() , AllDiags2.end());

        // Multiplying all the existing coefficients with the new coefficients!
        for(int i = 0; i < AllCoes2.size(); i++){
            for(int k = 0; k < AllCoes2[i].size() ; k++){
                AllCoes2[i][k] = AllCoes2[i][k] * CurrentCs[j];
            }
        }
        AllCoeffsOnParticle.insert(AllCoeffsOnParticle.end() , AllCoes2.begin() , AllCoes2.end());
    }
}

template<typename T>
vector<T> Concat_one(vector<T> A , vector<T> B){
    vector<T> output = A;
    for(int i = 0; i < B.size(); i++){
        output.push_back(B[i]);
    }
    return output;
}

template<typename T>
vector<vector<T>> Concat_two(vector<vector<T>> A , vector<vector<T>> B){
    vector<vector<T>> output = A;
    for(int i = 0; i < B.size(); i++){
        output.push_back(B[i]);
    }
    return output;
}

void Diag_Multiply(vector<vector<TotalDiag>>& ExistingDiags , vector<complex<double>>& ExistingCoeffs , vector<vector<TotalDiag>> NewDiags , vector<complex<double>> NewCoeffs){
    vector<vector<TotalDiag>> ExistingDiags0 = ExistingDiags;
    vector<complex<double>> ExistingCoeffs0 = ExistingCoeffs;
    int NumExisting = ExistingDiags.size() , NumNew = NewDiags.size(); 
    ExistingDiags.clear();
    ExistingCoeffs.clear();
    for(int i = 0; i < NumExisting; i++){
        for(int j = 0; j < NumNew; j++){
            ExistingDiags.push_back(Concat_one(ExistingDiags0[i] , NewDiags[j]));
            ExistingCoeffs.push_back(ExistingCoeffs0[i]*NewCoeffs[j]);
        }
    }
}

// This function converts the string of diagonals on a single particle into a string of total diagonals with the particle number appended
vector<vector<TotalDiag>> Diag_convert(vector<vector<ParticleDiag>> D , int ParticleNo){
    vector<vector<TotalDiag>> TotalDs;
    for(int i = 0; i < D.size(); i++){
        vector<TotalDiag> TotalDsi;
        for(int j = 0; j < D[i].size(); j++){
            TotalDiag TotalDsij;
            TotalDsij.ztype = D[i][j].first;
            TotalDsij.k = D[i][j].second;
            TotalDsij.particle = ParticleNo;
            if(TotalDsij.ztype != -1)
                TotalDsi.push_back(TotalDsij);
        }
        TotalDs.push_back(TotalDsi);
    }
    return TotalDs;
}

void PMR_otimes(TotalPerms& PermutationSet , TotalDVecs& DiagonalSet , Coeffs& CoeffSet , vector<int> NewPermutations , ParticleDVecs NewDiagonals, Coeffs NewCoeffs, int NewParticleNo , int twoSplusOne){
    TotalPerms PermutationSet0 = PermutationSet;
    TotalDVecs DiagonalSet0 = DiagonalSet;
    Coeffs CoeffSet0 = CoeffSet;
    int len = PermutationSet.size();
    PermutationSet.clear();
    DiagonalSet.clear();
    CoeffSet.clear();

    for(int j=0; j < NewPermutations.size(); j++){
        TotalPerms PermutationSetj = PermutationSet0;
        TotalDVecs DiagonalSetj = DiagonalSet0;
        Coeffs CoeffSetj = CoeffSet0;
        for(int k=0; k < len; k++){
            PermutationSetj[k].push_back({NewParticleNo , NewPermutations[j]}); // Here, we are multiplying the permutation matrices
            // Convert individual diagonal to total diagonal:
            vector<vector<TotalDiag>> NewDj = Diag_convert(NewDiagonals[j] , NewParticleNo); 
            Diag_Multiply(DiagonalSetj[k] , CoeffSetj[k] , NewDj , NewCoeffs[j]); // Here, we are multiplying the diagonal matrices with corresponding coefficients!
        }
        PermutationSet.insert(PermutationSet.end() , PermutationSetj.begin() , PermutationSetj.end());
        DiagonalSet.insert(DiagonalSet.end() , DiagonalSetj.begin() , DiagonalSetj.end());
        CoeffSet.insert(CoeffSet.end() , CoeffSetj.begin() , CoeffSetj.end());
    }
}

vector<pair<int,int>> Strip_zeros(vector<pair<int,int>> A){
    vector<pair<int,int>> A_stripped;
    for(int i=0; i < A.size(); i++){    
        if(A[i].second != 0){
            A_stripped.push_back(A[i]);
        }
    }
    return A_stripped;
}

// This function checks whether the permutation vector A is identity
bool Is_identity(vector<pair<int,int>> A){
    if(A.size() == 0){
        return true;
    }
    for(int i=0; i < A.size(); i++){
        if(A[i].second != 0){
            return false;
        }
    }
    return true;
}

void Perm_sort(vector<pair<int, int>>& vec) {
    // Sort the vector using a custom comparator function
    sort(vec.begin(), vec.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
}

bool Perm_compare(vector<pair<int,int>> A , vector<pair<int,int>> B){
    if(A.size() != B.size())
        return false;
    else{
        Perm_sort(A);
        Perm_sort(B);
        for(int i = 0; i < A.size(); i++){
            if(A[i]!= B[i]){
                return false;
            }
        }
    }
    return true;
}

pair<bool , int> Find_permutation(vector<pair<int,int>> LinePerms , vector<vector<pair<int,int>>> AllPerms){
    for(int i=0; i<AllPerms.size(); i++){
        if(Perm_compare(LinePerms , AllPerms[i])){
            return {true, i};
        }
    }
    return {0,0};
}


bool Diag_equal(TotalDiag Da , TotalDiag Db){
    if(Da.k != Db.k){
        return false;
    }
    else if(Da.ztype != Db.ztype){
        return false;
    }
    else if(Da.particle != Db.particle){
        return false;
    }
    return true;
}

bool Diag_compare(vector<TotalDiag> DiagA , vector<TotalDiag> DiagB){
    vector<TotalDiag> DiagB2 = DiagB;
    bool DiagAiMatched;
    if(DiagA.size() != DiagB.size()){
        return false;
    }
    for(int i = 0; i < DiagA.size(); i++){
        DiagAiMatched = false;
        //DiagB2 = DiagB;
        for(int j = 0; j < DiagB2.size(); j++){
            if(Diag_equal(DiagA[i] , DiagB2[j])){
                DiagB2.erase(DiagB2.begin() + j);
                DiagAiMatched = true;
                break;
            }
        }
        if(!DiagAiMatched){
            return false;
        }
    }
    return true;
}

// Test out these Diag compares!

void Diag_append(vector<vector<TotalDiag>>& AllDiagsi , vector<complex<double>>& AllCoesi , vector<vector<TotalDiag>> NewDiags , vector<complex<double>> NewCoes){
    for(int i = 0; i < NewDiags.size() ; i++){
        bool NewDiagiFound = false;
        for(int j = 0; j < AllDiagsi.size(); j++){
            if(Diag_compare(NewDiags[i] , AllDiagsi[j])){
                AllCoesi[j] += NewCoes[i];
                NewDiagiFound = true;
                if(abs(AllCoesi[j]) < 1E-8){
                    AllDiagsi.erase(AllDiagsi.begin() + j);
                    AllCoesi.erase(AllCoesi.begin() + j);
                }
                break;
            }
        }
        if(!NewDiagiFound){
            AllDiagsi.push_back(NewDiags[i]);
            AllCoesi.push_back(NewCoes[i]);
        }
    }
}

void PMR_append(PDdata& pdData, PDdata pdDataLine){
    TotalPerms AllPerms = pdData.Permutations , LinePerms = pdDataLine.Permutations;
    TotalDVecs AllDiags = pdData.Diagonals , LineDiags = pdDataLine.Diagonals;
    Coeffs AllCoes = pdData.Coefficients , LineCoes = pdDataLine.Coefficients;
    vector<vector<TotalDiag>> D0 = pdData.D0 , D0Line = pdDataLine.D0;
    vector<complex<double>> D0Coeffs = pdData.D0Coeffs , D0CoeffsLine = pdDataLine.D0Coeffs;

    for(int i=0; i < LinePerms.size(); i++){
        pair<bool , int> PermFound = Find_permutation(Strip_zeros(LinePerms[i]) , AllPerms);
        if(PermFound.first){
            // Do diag found and add coefficient
            Diag_append(AllDiags[PermFound.second] , AllCoes[PermFound.second] , LineDiags[i], LineCoes[i]);
        }
        else{
            AllPerms.push_back(Strip_zeros(LinePerms[i]));
            AllDiags.push_back(LineDiags[i]);
            AllCoes.push_back(LineCoes[i]);
        }
    }

    for(int i=0; i < D0Line.size();i++){
        D0.push_back(D0Line[i]);
        D0Coeffs.push_back(D0CoeffsLine[i]);
    }

    // Remove permutations which have no coefficients (coefficients have been canceled out!)
    vector<int> RemovedIndices = {};
    for(int i = 0; i < AllCoes.size(); i++){
        if(AllCoes[i].size() < 1E-8){
            RemovedIndices.push_back(i);
        }
    }

    for(int i=0; i < RemovedIndices.size(); i++){
        AllPerms.erase(AllPerms.begin() + RemovedIndices[i] - i);
        AllCoes.erase(AllCoes.begin() + RemovedIndices[i] - i);
        AllDiags.erase(AllDiags.begin() + RemovedIndices[i] - i);
    }

    pdData.Permutations = AllPerms;
    pdData.Coefficients = AllCoes;
    pdData.Diagonals = AllDiags;
    pdData.D0 = D0;
    pdData.D0Coeffs = D0Coeffs;
}

PDdata CDPconvert(const vector<pair<complex<double>,vector<int>>> data){
    int NumLines = data.size() , twoSplusOne;
    int NumOfParticles = 0;
    PDdata pdData;

    // Definining the diagonal matrices and their complex coefficients
    complex<double> one(1,0) , plusi(0,1) , minusi(0,-1);
    ParticleDiag Dplus = {0 , 0} , Dminus = {0 , -1} , Dz = {1 , 0} , Identity = {-1 , -1};

    // Defining the diagonal matrices with the pair convection. If the first integer is zero, we have D^{(k)}, and if the
    //      first integer is one, we have D^{(z,k)}. The second integer specifies k!
    for (int l = 0; l < NumLines; l++){
        complex<double> coeffl = data[l].first;
        vector<int> datal = data[l].second; // Extracts the array of spins and paulis for every line of input!
        vector<int> ParticlNos;
        vector<vector<int>> Permsl; // This is a vector (for each particle number) of (vector<int> The set of permutations acting on that particle)
        vector<ParticleDVecs> Diagsl; // The first vector is for each particle, the second vector is for the set of permutations, and the third vector is for the set of diagonal vectors for the specified permutation
        vector<Coeffs> Coeffsl; // Since we can get multiple permutation matrices per line of data, we need to keep track of the coefficients for each 

        TotalPerms PMatricesLine;  // The first in the pair is the set of particle numbers and the second vector is the set of powers of permutations
        TotalDVecs DMatricesLine; //This vector maps the Zs to Ps it is a many to one mapping!
        Coeffs CsLine;
        vector<vector<TotalDiag>> D0Line;
        vector<complex<double>> D0CoeffsLine;
        PDdata pdDataLine;

        int twoSpinplus1l; 
        for (size_t i = 0; i < datal.size()/4; i++) {
            int Particlei = datal[4 * i] , Pauli = datal[4 * i + 1] , Poweri = datal[4 * i + 2] , twoSpinPlusOnei = datal[4 * i + 3];
            if(Particlei <= 0){
                cout << endl << "Error: " << Particlei << " is incorrect particle index. Particle indices must be positive integers." << endl; exit(1);
            }
            twoSpinplus1l = twoSpinPlusOnei; // Assuming only a single spin species for now. In future, if there are multiple spin species, additional matrices should be created to separate these permutation operators!
            if (Particlei > NumOfParticles)
                NumOfParticles = Particlei;
            PauliCDPs Operator;
            if (Pauli == 1){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {one , one}}};  // Operator = S_x!
            } else if (Pauli == 2){
                Operator = {{1 , twoSpinPlusOnei - 1} , {{Dplus , Dminus} , {minusi , plusi}}}; // Operator = S_y!
            } else if (Pauli == 3){
                Operator = {{0} , {{Dz} , {one}}}; // Operator = S_z!
            } else if (Pauli == 0){
                Operator = {{0}, {{Identity} , {one}}};
            } else{
                cout << endl << "Error: incorrectly specified spin matrix" << endl; exit(1);
            }
            pair<bool, int> PartFound = Find_number(Particlei , ParticlNos);
            if(PartFound.first){
                while(Poweri > 0){
                    Permutation_append(Permsl[PartFound.second] , Diagsl[PartFound.second], Coeffsl[PartFound.second] , Operator , twoSpinPlusOnei);
                    Poweri--;
                }
            }
            else{
                ParticlNos.push_back(Particlei);
                Permsl.push_back(Operator.first);
                ParticleDVecs NewD;
                Coeffs NewC;
                for(int p = 0; p < Operator.second.first.size() ; p++){
                    ParticleDiag NewDp = Operator.second.first[p];
                    NewDp.second = (NewDp.second + twoSpinplus1l) % twoSpinplus1l; // Replacing k = -1 with k = 2S 
                    NewD.push_back({{NewDp}});
                    NewC.push_back({{Operator.second.second[p]}}); // This takes care of the front coefficient multiplication!
                }
                Diagsl.push_back(NewD);
                Coeffsl.push_back(NewC);
                int lastcoeff = Permsl.size()-1;
                while(Poweri > 1){
                    Permutation_append(Permsl[lastcoeff] , Diagsl[lastcoeff], Coeffsl[lastcoeff] , Operator , twoSpinPlusOnei);
                    Poweri--;
                }
            }
        }
        twoSplusOne = twoSpinplus1l;
        // Combining all the P and D matrices from a single line into corresponding PRM forms (there could be multiple PMR terms from each line)
        for(int k = 0; k < Permsl.size(); k++){
            if(k == 0){
                for(int kk = 0; kk < Permsl[k].size(); kk++){
                    PMatricesLine.push_back({{ParticlNos[k] , Permsl[k][kk]}});
                    DMatricesLine.push_back(Diag_convert(Diagsl[k][kk] , ParticlNos[k]));
                    CsLine.push_back(Coeffsl[k][kk]);
                }
            }
            else{
                PMR_otimes(PMatricesLine , DMatricesLine , CsLine , Permsl[k] , Diagsl[k] , Coeffsl[k] , ParticlNos[k] , twoSpinplus1l);
            }
        }
        // Multiply all coefficients by the coefficient at front:
        for(int c = 0; c < CsLine.size(); c++){
            for(int cc = 0; cc < CsLine[c].size(); cc++){
                CsLine[c][cc] = CsLine[c][cc]*coeffl;
            }
        }

        // Removing the identity permutation and extracting D0 and D0Coeffs
        vector<int> IdenPermsIndices;
        for(int i = 0; i < PMatricesLine.size(); i++){
            if(Is_identity(PMatricesLine[i])){
                D0Line.insert(D0Line.end() , DMatricesLine[i].begin() , DMatricesLine[i].end());
                D0CoeffsLine.insert(D0CoeffsLine.begin(), CsLine[i].begin() , CsLine[i].end());
                IdenPermsIndices.push_back(i);
            }
        }
        for(int i = 0; i < IdenPermsIndices.size(); i++){
            PMatricesLine.erase(PMatricesLine.begin() + IdenPermsIndices[i] - i);
            DMatricesLine.erase(DMatricesLine.begin() + IdenPermsIndices[i] - i);
            CsLine.erase(CsLine.begin() + IdenPermsIndices[i] - i);
        }

        // Going over the permutations and simplifying
        // Adding the new CDP from the line to the entire set of total CDPs
        pdDataLine.Coefficients = CsLine;
        pdDataLine.Permutations = PMatricesLine;
        pdDataLine.Diagonals = DMatricesLine;
        pdDataLine.D0 = D0Line;
        pdDataLine.D0Coeffs = D0CoeffsLine;

        PMR_append(pdData, pdDataLine);
    }

    pdData.twoSplusOne = twoSplusOne;
    pdData.NumOfParticles = NumOfParticles;

    return pdData;
}

vector<vector<int>> Convert_perms(vector<vector<pair<int,int>>> PMatrices , int NumOfParticles){
    vector<vector<int>> ColPermMatrix(NumOfParticles , vector<int>(PMatrices.size() , 0));
    for(int i=0; i< PMatrices.size(); i++){
        for(int j = 0; j < PMatrices[i].size(); j++){
            int particle = PMatrices[i][j].first , power = PMatrices[i][j].second;
            ColPermMatrix[particle-1][i] = power;
        }
    }
    return ColPermMatrix;
}

/*int main(int argc , char* argv[]){
    string fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> data = data_extract(fileName);
    cout << "Now processing ... " << endl;
    cout << endl;
    PDdata CDPdata = CDPconvert(data);
    Coeffs Cs = CDPdata.Coefficients;
    TotalDVecs DMatrices = CDPdata.Diagonals;
    vector< vector<pair<int,int>>> PMatrices = CDPdata.Permutations;
    int N = CDPdata.NumOfParticles;

    cout << "The following is the breakdown of the data " << endl;
    cout << endl;
    Print_data(CDPdata);

    // Converting the PMatrices into vector<vector<int>> to make matrix of column permutations

    vector<vector<int>> PermMatrixColumn = Convert_perms(PMatrices, N);
    cout << endl;
    printMatrix(PermMatrixColumn);


    return 0;
}*/
