#include "mex.h"
#include <math.h>

// For some compilers this is defined in math.h, for others its not
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void precompute_greens(double *r, double *it, double D, int nr, int nt, double *G_array) {
    double time_diff;
    for (int i = 0; i < nr; i++) {
        for (int t = 0; t < nt; t++) {
            for (int t_prime = 0; t_prime < t; t_prime++) {
                time_diff = it[t] - it[t_prime];
                G_array[i + nr * (t + nt * t_prime)] = pow(4 * M_PI * D * time_diff, -1.5) * exp(-r[i] * r[i] / (4 * D * time_diff));
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Check for proper number of arguments. */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:precompute_greens:invalidNumInputs", "Three input arguments required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:precompute_greens:invalidNumOutputs", "One output argument required.");
    }
    
    /* Input parameters */
    double *r = mxGetPr(prhs[0]);
    double *it = mxGetPr(prhs[1]);
    double D = mxGetScalar(prhs[2]);
    int nr = mxGetN(prhs[0]); /* Number of distances */
    int nt = mxGetN(prhs[1]); /* Number of time points */
    
    /* Create output matrix */
    mwSize dims[3] = {nr, nt, nt};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *G_array = mxGetPr(plhs[0]);
    
    /* Call the computational routine */
    precompute_greens(r, it, D, nr, nt, G_array);
}
