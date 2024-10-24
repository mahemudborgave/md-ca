#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>  // Include the chrono library

using namespace std;
using namespace std::chrono; // Use chrono namespace for easier access

const int N = 8; // Matrix size (8x8)

// Function to print a matrix
void printMatrix(double M[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

// Function to print vector
void printVector(double v[N]) {
    for (int i = 0; i < N; i++) {
        cout << v[i] << endl;
    }
}

// Function to print the permutation matrix P
void printPermutationMatrix(double P[N][N]) {
    cout << "Permutation matrix P:" << endl;
    printMatrix(P);
}

// Function to perform LUP decomposition
void LUPDecomposition(double A[N][N], double L[N][N], double U[N][N], double P[N][N]) {
    int i, j, k, maxIndex;
    double maxValue;

    // Initialize permutation matrix P as an identity matrix and L
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            P[i][j] = (i == j) ? 1 : 0; // Identity matrix for P
            L[i][j] = (i == j) ? 1 : 0; // Identity matrix for L
            U[i][j] = A[i][j]; // Start U as a copy of A
        }
    }

    printPermutationMatrix(P);

    for (k = 0; k < N; k++) {
        maxValue = 0.0;
        maxIndex = k;

        // Find the pivot element
        for (i = k; i < N; i++) {
            if (fabs(U[i][k]) > maxValue) {
                maxValue = fabs(U[i][k]);
                maxIndex = i;
            }
        }



        // Pivot if necessary
        if (maxIndex != k) {
            // Swap rows in permutation matrix P
            for (j = 0; j < N; j++) {
                swap(P[k][j], P[maxIndex][j]);
            }

            // Swap rows in U
            for (j = 0; j < N; j++) {
                swap(U[k][j], U[maxIndex][j]);
            }

            // Swap rows in L (only the first k columns)
            for (j = 0; j < k; j++) {
                swap(L[k][j], L[maxIndex][j]);
            }
        }

        // Perform Gaussian elimination to form L and U
        for (i = k + 1; i < N; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (j = k; j < N; j++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    printPermutationMatrix(P);
}

// Function to solve the system using LUP decomposition
vector<double> LUPSolve(double U[N][N], double L[N][N], double P[N][N], double b[N]) {
    double y[N], x[N];
    int i, j;

    // Apply permutation matrix P to b (Pb)
    double Pb[N] = {0};
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Pb[i] += P[i][j] * b[j]; // Multiply P by b
        }
    }

    // Forward substitution for Ly = Pb
    for (i = 0; i < N; i++) {
        y[i] = Pb[i];
        for (j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    // Back substitution for Ux = y
    for (i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return vector<double>(x, x + N);
}

int main() {
    // Input 8x8 matrix A
    double A[N][N] = {
        {10, -1, 2, 0, 0, 0, 0, 0},
        {-1, 11, -1, 3, 0, 0, 0, 0},
        {2, -1, 10, -1, 0, 0, 0, 0},
        {0, 3, -1, 8, 2, 0, 0, 0},
        {0, 0, 0, 2, 6, -1, -1, 0},
        {0, 0, 0, 0, -1, 9, 0, -1},
        {0, 0, 0, 0, -1, 0, 8, -1},
        {0, 0, 0, 0, 0, -1, -1, 9}
    };

    // Input 8-dimensional vector b
    double b[N] = {6, 25, -11, 15, -20, 12, 17, 10};

    // Permutation matrix P
    double P[N][N];
    
    // Matrices for L and U
    double L[N][N], U[N][N];

    // Measure time for LUP decomposition
    auto start = high_resolution_clock::now();  // Start timing

    // Perform LUP decomposition on matrix A
    LUPDecomposition(A, L, U, P);

    // Measure time for solving the system
    auto mid = high_resolution_clock::now();  // Time after decomposition
    vector<double> x = LUPSolve(U, L, P, b);
    auto end = high_resolution_clock::now();  // End timing

    // Calculate elapsed time
    auto durationDecomp = duration_cast<microseconds>(mid - start);
    auto durationSolve = duration_cast<microseconds>(end - mid);
    auto totalDuration = duration_cast<microseconds>(end - start);

    // Output the solution
    cout << "Solution vector x:" << endl;
    printVector(x.data());

    // Output L and U matrices
    cout << "Matrix L:" << endl;
    printMatrix(L);

    cout << "Matrix U:" << endl;
    printMatrix(U);

    // Output elapsed time
    cout << "Elapsed time for LUP decomposition: " << durationDecomp.count() << " microseconds" << endl;
    cout << "Elapsed time for solving the system: " << durationSolve.count() << " microseconds" << endl;
    cout << "Total elapsed time: " << totalDuration.count() << " microseconds" << endl;

    return 0;
}
