#include "mex.h"

void approximate_solution(int it, int nr, double dt, double* G_array, int nt, double* I, double* C) {
    int i, j;

    for (i = 0; i < nr; i++) { // distance
        C[i] = 0.0;
        for (j = 0; j < it; j++) { // time
            /*
            i : index for the first dimension (distance r).
            it : index for the second dimension (current time step).
            j : index for the third dimension (previous time step).
            
            Given the dimensions:
            nr : number of distances
            nt : number of time steps
            */
            C[i] += I[j] * G_array[i + (it * nr) + (j * nr * nt)] * dt;
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check for proper number of arguments
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:approximate_solution:invalidNumInputs", "Six input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:approximate_solution:maxlhs", "Too many output arguments.");
    }

    // Input variables
    int it = mxGetScalar(prhs[0]);
    int nr = mxGetScalar(prhs[1]);
    double* I = mxGetPr(prhs[2]);
    double dt = mxGetScalar(prhs[3]);
    double* G_array = mxGetPr(prhs[4]);
    int nt = mxGetScalar(prhs[5]);

    // Output variable
    plhs[0] = mxCreateDoubleMatrix(1, nr, mxREAL);
    double* C = mxGetPr(plhs[0]);

    // Call the computation function
    approximate_solution(it, nr, dt, G_array, nt, I, C);
}
