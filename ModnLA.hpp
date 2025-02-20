#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;


template<typename T>
void Print_matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<int>> Modp_matrixadd(vector<vector<int>> A , vector<vector<int>> B , int p){
    int lA1 = A.size() , lB1 = B.size() , lA2 = A[0].size() , lB2 = B[0].size();
    if(lA1 != lB1 || lA2 != lB2){
        cout << "Matrix addition cannot be performed correctly, as the matrices provided do not have identical dimensions." << endl;
        return {{}};
    }

    vector<vector<int>> ApB = A;
    for(int i = 0; i < lA1; ++i){
        for(int j = 0; j < lA2; ++j){
            ApB[i][j] = (ApB[i][j] + B[i][j])%p;
        }
    }
    return ApB;
}

// This function mod p multiplies the matrix A by a scalar c!
vector<vector<int>> Modp_scalmult(const int c , const vector<vector<int>> A , int p){
    int rowsA = A.size();
    if(rowsA == 0){
        return A;
    }
    int colsA = A[0].size();
    vector<vector<int>> cA(rowsA , vector<int> (colsA , 0));
    for(int i = 0; i < rowsA ; i++){
        for(int j = 0; j < colsA ; j++){
            if(c >= 0){
                cA[i][j] = c*A[i][j] % p;
            }
            else{
                cA[i][j] = (p + (c*A[i][j] % p)) % p;
            }
        }
    }
    return cA;
}

// This function mod p multiplies the matrix A by a scalar c!
vector<vector<int>> Modp_scalmult_double(const double c , const vector<vector<int>> A , int p){
    int rowsA = A.size();
    if(rowsA == 0){
        return A;
    }
    int colsA = A[0].size();
    if(colsA == 0){
        return A;
    }
    vector<vector<int>> cA(rowsA , vector<int> (colsA , 0));
    for(int i = 0; i < rowsA ; i++){
        for(int j = 0; j < colsA ; j++){
            int num = c*A[i][j];
            if(c >= 0){
                cA[i][j] = num % p;
            }
            else{
                cA[i][j] = (p + (num % p)) % p;
            }
        }
    }
    return cA;
}

// This function computes matrix multiplication of A and b mod p.
vector<vector<int>> Modp_mult(const vector<vector<int>> A , const vector<vector<int>> B , int p){
    int rowsA = A.size(), colsA = A[0].size() , rowsB = B.size() , colsB;

    if (colsA != rowsB) {
        if(rowsB == 0){
            return {};
        }
        else{
            colsB = B[0].size();
            if(colsB == 0){
                return {{}};
            }
            else{
                cerr << "Matrices are not compatible for multiplication." << std::endl;
                return {};
            }
        }
    }
    colsB = B[0].size();

    vector<vector<int>> AB(rowsA, vector<int>(colsB, 0));
    
    // Perform element-wise multiplication and summation
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                AB[i][j] += A[i][k] * B[k][j];
            }
            AB[i][j] = AB[i][j] % p;
        }
    }

    return AB;
}

// This function divides num by divider (mod p): this only works if p is prime, since multiplication cannot produce division when p is not a prime!
int Modp_divide(int num , int divider , int p){
    int answer;
    for (int i = 0; i < p; i++){
        if( (divider * i) % p == num ){
            answer = i;
        }
    }
    return answer;
}

// This function transposes matrix A
vector<vector<int>> Transpose(vector<vector<int>> A){
    vector<vector<int>> A_transpose;
    int rowsA = A.size();
    if(rowsA == 0){
        return A_transpose;
    }
    else{
        int colsA = A[0].size();
        for(int j=0; j < colsA ; j++){
            vector<int> A_transpose_j(rowsA , 0);
            for(int i = 0; i < rowsA ; i++){
                A_transpose_j[i] = A[i][j];
            }
            A_transpose.push_back(A_transpose_j);
        }
        return A_transpose;
    }
}

// Divide_vector divides the given matrix in 2 vertically!
vector<vector<vector<int>>> Divide_vector(const vector<vector<int>> Vec , int vecDim){
    int rowsVec = Vec.size();
    vector<vector<int>> Vec1, Vec2;
    vector<vector<vector<int>>> VecDiv;

    for (int i = 0; i < rowsVec ; i++){
        if(i < vecDim){
            Vec1.push_back(Vec[i]);
        }
        else{
            Vec2.push_back(Vec[i]);
        }
    }

    VecDiv.push_back(Vec1);
    VecDiv.push_back(Vec2);
    
    return VecDiv;
}

vector<vector<vector<int>>> Divide_blocks(const vector<vector<int>> A , int blockDim){
    int rowsA = A.size(), colsA = A[0].size();

    vector<vector<vector<int>>> BlockMatrix;
    vector<vector<int>> TopLeft(blockDim , vector<int> (blockDim , 0)), BottomRight(rowsA - blockDim , vector<int> (rowsA - blockDim , 0));
    vector<vector<int>> TopRight(blockDim , vector<int> (rowsA - blockDim , 0));

    for (int i = 0; i < blockDim; i++) {
        for (int j = i; j < blockDim; j++) {
            TopLeft[i][j] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = i; j < rowsA; j++) {
            BottomRight[i - blockDim][j - blockDim] = A[i][j];
        }
    }
    for (int i = 0; i < blockDim; i++) {
        for (int j = blockDim; j < rowsA; j++) {
            TopRight[i][j - blockDim] = A[i][j];
        }
    }

    BlockMatrix.push_back(TopLeft);
    BlockMatrix.push_back(TopRight);
    BlockMatrix.push_back(BottomRight);

    return BlockMatrix;
}

vector<vector<vector<int>>> Divide_blocks_lower(const vector<vector<int>> A , int blockDim){
    int rowsA = A.size(), colsA = A[0].size();

    vector<vector<vector<int>>> BlockMatrix;
    vector<vector<int>> TopLeft(blockDim , vector<int> (blockDim , 0)), BottomRight(rowsA - blockDim , vector<int> (rowsA - blockDim , 0));
    vector<vector<int>> BottomLeft(rowsA-blockDim , vector<int> (blockDim , 0));

    for (int i = 0; i < blockDim; i++) {
        for (int j = 0; j <= i; j++) {
            TopLeft[i][j] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = 0; j <= i; j++) {
            BottomRight[i - blockDim][j - blockDim] = A[i][j];
        }
    }
    for (int i = blockDim; i < rowsA; i++) {
        for (int j = 0; j < blockDim; j++) {
            BottomLeft[i-blockDim][j] = A[i][j];
        }
    }

    BlockMatrix.push_back(TopLeft);
    BlockMatrix.push_back(BottomLeft);
    BlockMatrix.push_back(BottomRight);

    return BlockMatrix;
}

// This function vertically concatenates to matrices X1 and X2
vector<vector<int>> Vert_conc(vector<vector<int>> X1 , vector<vector<int>> X2){
    vector<vector<int>> X = X1;
    for(int i = 0; i < X2.size(); i++){
        X.push_back(X2[i]);
    }
    return X;
}

// This function horizontally concatenates matrices X1 and X2
vector<vector<int>> Horz_conc(vector<vector<int>> X1 , vector<vector<int>> X2){
    vector<vector<int>> X = X1;
    if(X1.size() == 0){
        return X2;
    }
    else if(X2.size() == 0){
        return X1;
    }
    else{
        for(int i = 0; i < X.size(); i++){
            for(int j = 0; j < X2[0].size(); j++){
                X[i].push_back(X2[i][j]);
            }
        }
        return X;
    }
}

vector<vector<int>> Vec_sub(vector<vector<int>> Vec1, vector<vector<int>> Vec2 , int p){
    vector<vector<int>> Vecdiff;
    int diff;
    for(int i =0; i < Vec1.size() ; i++){
        vector<int> Vecdiff_i = {};
        for(int j = 0; j < Vec1[0].size(); j++){
            diff = (Vec1[i][j] - Vec2[i][j])%p;
            if(diff < 0){
                diff += p;
            }
            Vecdiff_i.push_back(diff);
        }
        Vecdiff.push_back(Vecdiff_i);
    }
    return Vecdiff;
}

vector<vector<int>> Vec_add(vector<vector<int>> Vec1, vector<vector<int>> Vec2 , int p){
    vector<vector<int>> Vecsum;
    int sum;
    for(int i =0; i < Vec1.size() ; i++){
        vector<int> Vecsum_i = {};
        for(int j = 0; j < Vec1[0].size(); j++){
            sum = (Vec1[i][j] + Vec2[i][j])%p;
            Vecsum_i.push_back(sum);
        }
        Vecsum.push_back(Vecsum_i);
    }
    return Vecsum;
}

bool All_zeros(vector<vector<int>> A){
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[0].size(); j++){
            if(A[i][j] != 0){
                return false;
            }
        }
    }
    return true;
}

bool Num_found_inVec(int num , vector<int> Vec){
    bool found = false;
    for (int i=0; i < Vec.size(); i++){
        if(num == Vec[i]){
            found = true;
            break;
        }
    }
    return found;
}

typedef vector<vector<int>> Matrix;

// Mod p Gaussian Elimination:
void Modp_GE(Matrix &A , int p) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int lead = 0;

    for (int r = 0; r < rowsA; r++) {
        if (lead >= colsA) {
            break;
        }
        int i = r;
        while (A[i][lead] == 0) {
            i++;
            if (i == rowsA) {
                i = r;
                lead++;
                if (lead == colsA) {
                    return;
                }
            }
        }
        
        swap(A[i], A[r]);

        int lv = A[r][lead] , multiplier = Modp_divide(1, A[r][lead], p); // Multiplier is the number that multiplies the leading term to 1.
        for (int j = 0; j < colsA; j++) {
            A[r][j] = (A[r][j] * multiplier) % p;
        }
        for (int i = 0; i < rowsA; i++) {
            if (i != r) {
                int lv2 = A[i][lead];
                for (int j = 0; j < colsA; j++) {
                    A[i][j] = (A[i][j] - (A[r][j] * lv2) % p + p) % p;
                }
            }
        }
        lead++;
    }
}

int Sum(vector<int> A, int p){
    int sum = 0;
    for(int i=0; i < A.size() ; i++){
        sum =+ A[i];
    }
    return sum%p;
}

// The following function decomposes an integer into its prime factors
vector<pair<int , int>> Prime_decomp(int n){
    vector<pair<int, int>> primes;

    int power2 = 0;
    if(n%2 == 0){
        while(n % 2 == 0){
            n /= 2;
            power2++;
        }
        primes.push_back({2 , power2});
    }

    for(int i = 3; i*i < n ; i += 2){
        int poweri = 0;
        if(n%i == 0){
            while(n % i == 0){
                n /= i;
                poweri++;
            }
            primes.push_back({i , poweri});
        }
    }

    // In case n is itself a prime:
    if(n > 1){
        primes.push_back({n , 1});
    }

    return primes;
}

typedef vector<vector<int>> Matrix;

int modInverse(int a, int mod) {
    // Assumes mod is a prime number
    int m = mod - 2;  // Fermat's little theorem
    int result = 1;
    while (m) {
        if (m & 1)
            result = (result * a) % mod;
        a = (a * a) % mod;
        m >>= 1;
    }
    return result;
}

bool Not_null(const vector<int> v) {
    return any_of(v.begin(), v.end(), [](int i) { return i != 0; });
}

bool Check_null(Matrix Perms , vector<int> null , int modn){
    Matrix PermTrans = Transpose(Perms);
    vector<int> eig = Modp_scalmult(null[0] , {PermTrans[0]} , modn)[0];
    for(int i = 1; i < null.size() ; ++i){
        eig = Vec_add({eig} , {Modp_scalmult(null[i] , {PermTrans[i]} , modn)} , modn)[0];
    }
    return !Not_null(eig);
}

Matrix Modp_nullspace(const Matrix& A, int p) {
    int n = A.size();        // number of rows
    int m = A[0].size();     // number of columns
    Matrix B = A;            // Make a copy of A to transform into RREF
    Matrix nullspace;
    
    for(int i = 0; i < n ; ++i){
        for(int j = 0; j < m; ++j){
            B[i][j] = B[i][j]%p;
        }
    }

    
    if(All_zeros(B)){
        return {{0}};
    }


    // Step 1: Transform B into RREF form mod p
    int row = 0;
    vector<int> lead_col(m, -1);  // Track leading columns
    for (int col = 0; col < m && row < n; col++) {
        int sel = row;
        // Find a non-zero entry in the current column
        for (int i = row; i < n; i++) {
            if (B[i][col] % p != 0) {
                sel = i;
                break;
            }
        }

        if (B[sel][col] % p == 0)
            continue;  // Skip this column if no non-zero entries

        // Swap rows if necessary
        if (sel != row) {
            swap(B[sel], B[row]);
        }

        // Scale row to make the leading coefficient 1
        int inv = modInverse(B[row][col], p);
        for (int j = 0; j < m; j++) {
            B[row][j] = (B[row][j] * inv) % p;
            if (B[row][j] < 0) B[row][j] += p;
        }

        // Zero out all other entries in this column
        for (int i = 0; i < n; i++) {
            if (i != row) {
                int c = B[i][col];
                for (int j = 0; j < m; j++) {
                    B[i][j] = (B[i][j] - c * B[row][j] % p + p) % p;
                }
            }
        }

        lead_col[row] = col;  // Mark the leading column
        row++;
    }
    // Step 2: Identify free variables and form nullspace basis vectors
    nullspace.clear();
    vector<int> free_var_index;
    for (int j = 0; j < m; j++) {
        if (find(lead_col.begin(), lead_col.end(), j) == lead_col.end()) {
            free_var_index.push_back(j);  // j is a free variable
        }
    }

    for (int idx : free_var_index) {
        vector<int> basis_vec(m, 0);
        basis_vec[idx] = 1;  // Set the free variable
        for (int i = 0; i < row; i++) {
            if (lead_col[i] != -1) {
                basis_vec[lead_col[i]] = -B[i][idx];
                if (basis_vec[lead_col[i]] < 0)
                    basis_vec[lead_col[i]] += p;
            }
        }
        nullspace.push_back(basis_vec);
    }

    return Transpose(nullspace);
}


// The following function solves the system of equations Ax = b (mod p)
// Even though b is a vector<vector<int>>, it must be only a single column vector!
Matrix Modp_solver(Matrix A , Matrix b , int p){
    
    vector<vector<int>> Ab = Horz_conc(A , Modp_scalmult(-1 , b , p));
    Modp_GE(Ab , p);

    Matrix x = Modp_nullspace(Ab , p) , x_final;
    Matrix xt = Transpose(x);

    int n = xt[0].size() , m = xt.size();
    for(int i = 0; i < m; i++){ 
        if(xt[i][n-1] != 0){
            int mult = Modp_divide(xt[i][n-1] , 1 , p);
            vector<int> x_i = Modp_scalmult(mult , {xt[i]} , p)[0];

            x_i.erase(x_i.begin() + n - 1);
            x_final.push_back(x_i);
        }
    }
    return Transpose(x_final);
}

// This function computes nullspace basis of A mod p^r! when r = 1, it simply makes a single call to Modp_Nullspace(A , p)
Matrix Modp_nullspace_r(Matrix A , int p , int r){
    if(r == 1){
        return Modp_nullspace(A , p);
    }
    else{
        Matrix x = Modp_nullspace_r(A , p , r-1) , xt , b;
        xt = Transpose(x);
        double factor = -1.0/pow(p , r-1);
        // In case all vectors are in nullspace when computing in lower mod
        bool IdentityNull = false;
        if(x.size() > 0){
            if(x.size() == 1 & x[0].size() == 1){
                IdentityNull = true;
            }
        }

        if(IdentityNull){
            // This means that no nullspace was found, and one must divide by p!
            int m = A[0].size();
            Matrix nulliden(m , vector<int> (m , 0));
            for(int i = 0; i < m; ++i){
                nulliden[i][i] = pow(p , r-1);
            }
            x = Horz_conc(nulliden , Modp_nullspace_r(Modp_scalmult_double( 1.0/p , A , pow(p,r)) , p , r-1));
            //xt = Transpose(x);
            //factor = factor*p;
            //r++;
            return x;
        }

        b = Modp_mult(A , x , pow(p , r));
        b = Modp_scalmult_double( factor , b , pow(p , r));
        Matrix btrans = Transpose(b) , xrs;

        if(!All_zeros(btrans) & btrans.size() > 0){
            for(int i = 0; i < btrans.size(); i++){
                Matrix xr = Modp_solver(A , Transpose({btrans[i]}) , p);
                if (Transpose(xr).size() == 1){
                    // add p*xr with xi
                    xt[i] = Vec_add({xt[i]} , Modp_scalmult(pow(p , r-1) , Transpose(xr) , pow(p , r)) , pow(p,r))[0];
                }
            }
        }
        x = Horz_conc(Modp_scalmult(p , x , pow(p , r)) ,  Transpose(xt));
        return x;
    }
}

// This function determines if the input vector is divisible by m (i.e. all entries are multiples of m)
bool Vec_divisible(vector<int> vec , int m){
    for(int i=0; i < vec.size();++i){
        if(vec[i]%m != 0){
            return false; 
        }
    }
    return true;
}

pair<Matrix , vector<int>> Find_divisibles(const Matrix A , int p , int r){
    // Assuming we are given column vectors (A), in order to find if column is a multiple of p^k,
    //     we will transform into rows and determine if rows are multiples of p^k

    // This function returns the column matrix of all divisible columns and a vector<int> of their indices
    //      corresponding to the original columns of A

    // This function would return a non-empty matrix only if there are subsets of A that have lower dimensions than A,
    //    i.e. the function would return empty even if the entire A is itself a multiple!
    int m = pow(p , r);
    vector<int> Indices;
    Matrix At = Transpose(A) , Divisibles;
    for(int i=0; i < At.size();++i){
        if(Vec_divisible(At[i] , m)){
            Divisibles.push_back(At[i]);
            Indices.push_back(i);
        }
    }

    if(Divisibles.size() == At.size() & r < 2){
        return { Matrix {} , vector<int> {} };
    }

    return {Transpose(Divisibles) , Indices};
}

Matrix Modp_rnull_full(const Matrix A , int p , int r){
    Matrix FullNulls = Modp_nullspace_r(A , p , r);
    for(int i=1; i < r; i++){
        int n = pow(p , i);
        pair<Matrix , vector<int>> ApandInd = Find_divisibles(A , p , i);
        Matrix Ap = ApandInd.first , ApNull;
        vector<int> Indices = ApandInd.second;
        if(Ap.size() > 0){
            ApNull = Modp_nullspace_r( Modp_scalmult_double(1.0/n , Ap , pow(p,r)) , p , r - i);
        }
        // Building the actual nullspace eigenvector using the conversion indices of A and Ap
        Matrix ApNullt = Transpose(ApNull) , FullNullp;
        for(int j=0; j < ApNullt.size(); ++j){
            vector<int> eigvec(A[0].size() , 0);
            for(int k = 0; k < ApNullt[j].size(); ++k){
                int num = ApNullt[j][k];
                if(num != 0){
                    eigvec[Indices[k]] = num;
                }
            }
            FullNullp.push_back(eigvec);
        }
        FullNulls = Horz_conc(FullNulls, Transpose(FullNullp)); 
    }
    return FullNulls;
}


int GCD_vec(vector<int> vec) {
    int current_gcd = 0;  // Start with 0, which is neutral for GCD computation.
    for (int num : vec) {
        if (num != 0) {  // Only consider non-zero elements
            if (current_gcd == 0) {
                current_gcd = num;  // Initialize current_gcd with the first non-zero element
            } else {
                current_gcd = gcd(current_gcd, num);
            }
        }
    }
    return abs(current_gcd);  // Return the absolute value of the GCD to ensure it's non-negative
}

vector<int> Divide_vec( vector<int> vec, int num){
    vector<int> vec_div = vec;
    for(int i =0 ; i < vec.size() ; ++i){
        vec_div[i] = vec[i]/num;
    }
    return vec_div;
}

void Simplify_nullspace(Matrix PermMat , Matrix& Nulls , int modn){
    // Assumption here is that the rows of the Nulls are eigenvectors of the nullspace of the column matrix PermMat
    int Nrows = Nulls.size();
    for(int i=0; i < Nrows ; ++i){
        int gcdi = GCD_vec(Nulls[i]);
        if(gcdi > 1){
            vector<int> nullsimp_i = Divide_vec(Nulls[i] , gcdi);
            if(Check_null(PermMat , nullsimp_i , modn)){
                Nulls[i] = nullsimp_i;
            }
        }
    }
}

void Remove_redundancy(Matrix& A , int n){
    // Assumption here is that matrix A is a column matrix of nullspace basis

    Matrix At = Transpose(A);
    Matrix Atcop = At;
    vector<int> ToBeRemoved;
    // while keep adding columns
    for(int i = 0; i < At.size(); ++i){
        for(int j = i + 1; j < At.size() ; ++j){
            if(All_zeros( Vec_sub({Atcop[j]} , {At[i]} , n) ) ){
                ToBeRemoved.push_back(j);
            }
        }
    }

    for(int i = 0; i < ToBeRemoved.size(); i++){
        At.erase(At.begin() + ToBeRemoved[i]);
    }
    A = Transpose(At);
}

// This function computes the mod n nullspace basis of A, where n is any integer!
Matrix Nullspace_n(const Matrix A , int n){
    vector<pair<int, int>> ps = Prime_decomp(n);
    Matrix Null;

    if(A.size() > 0){
        for(int i=0 ; i < ps.size() ; i++){
            int prime = ps[i].first , power = ps[i].second;
            //Matrix Nullspacei = Transpose(Modp_rnull_full(A , prime , power)) , NullsiValid;
            Matrix Nullspacei = Transpose(Modp_rnull_full(A , prime , power)) , NullsiValid;

            // Checking if eigenvectors are valid , as for non-primes, there could be cases of invalid eigenvectors produced
            if(power > 1){
                for(auto& vec : Nullspacei){
                    if(Check_null(A , vec , pow(prime , power))){
                        NullsiValid.push_back(vec);
                    }
                }
            }
            else{
                NullsiValid = Nullspacei;
            }
            NullsiValid = Modp_scalmult(n/pow(prime , power) , NullsiValid, n);
            Null = Vert_conc(Null , NullsiValid);
        }
    }

    Simplify_nullspace(A , Null , n);
    Null = Transpose(Null);
    Remove_redundancy(Null , n);
    return Null;
}



/*
vector<int> Matvec_mult(const Matrix& mat, const vector<int>& vec, int mod) {
    int n = mat.size();     // number of rows
    int m = vec.size();     // should match the number of columns in mat
    vector<int> result(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[i] = (result[i] + mat[i][j] * vec[j]) % mod;
        }
    }
    return result;
}

void RefineBasis(const Matrix& A, Matrix& basis, int p, int pi, int pii) {
    Matrix newBasis;
    for (auto& vec : basis) {
        std::vector<int> w = Matvec_mult(A, vec, pii);  // A * v mod p^(i+1)
        for (int j = 0; j < w.size(); ++j) {
            if (w[j] % pi != 0) { // Check if congruent to 0 mod p^i
                for (int k = 0; k < vec.size(); ++k) {
                    int adjustment = (pii + w[j]) / pi;  // Adjust to make it zero mod p^(i+1)
                    vec[k] = (vec[k] - adjustment * vec[k]) % pii;
                }
            }
        }
        newBasis.push_back(vec);
    }
    basis = newBasis;
}

Matrix Modp_nullspace_r_new(const Matrix& A, int p, int r) {
    Matrix result;
    Matrix basis = Modp_nullspace_new(A, p);  // Start with the nullspace mod p

    int pi = p;
    for (int i = 1; i < r; ++i) {
        int pii = pi * p;
        RefineBasis(A, basis , p, pi, pii);
        pi = pii;
    }

    result = basis;

    return result;
}
*/

/*int main(){
    // Testing Nullspace finder!
    vector<vector<int>> A = Transpose({{1 , 2 , 1} , {2 , 1 , 0} , {0 , 2  , 1} , {1 , 0 , 1}}), Arr;
    vector<vector<int>> Iden = Transpose({{1 , 0 , 0 , 0} , {0 , 1 , 0 , 0} , {0 , 0  , 1 , 0} , {0 , 0 , 0 , 1}});
    vector<vector<int>> Empty = Transpose({{}});

    Arr = A;
    int p = 7;
    Modp_GE(Arr , p); // Gaussian Elimination performed on A!
    cout << "A is " << endl;
    Print_matrix(A);
    cout << endl;

    cout << "A in row echelon form " << endl;
    Print_matrix(Arr);

    cout << endl;

    cout << "The nullspace calculation is starting:" << endl;
    vector<vector<int>> x = Nullspace_n(A , p);
    cout << endl;
    cout << "The nullspace basis of A " << endl;
    Print_matrix(x);

    return 0;
}*/
