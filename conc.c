#include "mex.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void compute_concentration(double *C, const double *I, const double *r, const double *t, double D, double dt, int nr, int nt, int it) {
    for (int j = 0; j < nr; j++) {
        for (int k = 0; k <= it; k++) {
            double tau = t[k];
            double delta_t = t[it] - tau;
            if (delta_t > 0) {
                double term1 = I[k] / pow(4 * M_PI * D * delta_t, 1.5);
                double term2 = erfc(-pow(r[j], 2) / (4 * D * delta_t));
                C[it + j * nt] += term1 * term2 * dt;
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs
    double *C = mxGetPr(prhs[0]);
    const double *I = mxGetPr(prhs[1]);
    const double *r = mxGetPr(prhs[2]);
    const double *t = mxGetPr(prhs[3]);
    double D = mxGetScalar(prhs[4]);
    double dt = mxGetScalar(prhs[5]);
    int nr = mxGetN(prhs[2]);
    int nt = mxGetN(prhs[3]);
    int it = (int)mxGetScalar(prhs[6]);

    // Call the compute_concentration function
    compute_concentration(C, I, r, t, D, dt, nr, nt, it);
}