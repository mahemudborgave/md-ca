#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>  // Include the chrono library

using namespace std;
using namespace std::chrono; // Use chrono namespace for easier access

const int N = 8; // Matrix size (8x8)

void LUPDecomposition(double A[N][N], int P[N]) {
    int i, j, k, maxIndex;
    double maxValue, temp;

    for (i = 0; i < N; i++) {
        P[i] = i;  // Initialize permutation matrix
    }

    for (k = 0; k < N; k++) {
        maxValue = 0.0;
        maxIndex = k;

        for (i = k; i < N; i++) {
            if (fabs(A[i][k]) > maxValue) {
                maxValue = fabs(A[i][k]);
                maxIndex = i;
            }
        }

        // Pivot if necessary
        if (maxIndex != k) {
            for (j = 0; j < N; j++) {
                temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }

            temp = P[k];
            P[k] = P[maxIndex];
            P[maxIndex] = temp;
        }

        for (i = k + 1; i < N; i++) {
            A[i][k] /= A[k][k];
            for (j = k + 1; j < N; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

vector<double> LUPSolve(double A[N][N], int P[N], double b[N]) {
    double y[N], x[N];
    int i, j;

    // Forward substitution for Ly = Pb
    for (i = 0; i < N; i++) {
        y[i] = b[P[i]];
        for (j = 0; j < i; j++) {
            y[i] -= A[i][j] * y[j];
        }
    }

    // Back substitution for Ux = y
    for (i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (j = i + 1; j < N; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
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

    // Input 8-dimensional vector b (inequalities)
    double b[N] = {6, 25, -11, 15, -20, 12, 17, 10};

    // Permutation matrix P
    int P[N];

    // Measure time for LUP decomposition
    auto start = high_resolution_clock::now();  // Start timing

    // Perform LUP decomposition on matrix A
    LUPDecomposition(A, P);

    // Measure time for solving the system
    auto mid = high_resolution_clock::now();  // Time after decomposition
    vector<double> x = LUPSolve(A, P, b);
    auto end = high_resolution_clock::now();  // End timing

    // Calculate elapsed time
    auto durationDecomp = duration_cast<microseconds>(mid - start);
    auto durationSolve = duration_cast<microseconds>(end - mid);
    auto totalDuration = duration_cast<microseconds>(end - start);

    // Output the solution
    cout << "Solution vector x:" << endl;
    for (int i = 0; i < N; i++) {
        cout << x[i] << endl;
    }

    // Output elapsed time
    cout << "Elapsed time for LUP decomposition: " << durationDecomp.count() << " microseconds" << endl;
    cout << "Elapsed time for solving the system: " << durationSolve.count() << " microseconds" << endl;
    cout << "Total elapsed time: " << totalDuration.count() << " microseconds" << endl;

    return 0;
}
