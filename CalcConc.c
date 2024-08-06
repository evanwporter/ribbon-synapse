#include "mex.h"
#include <math.h>

/* 
Compile with the command
mex CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CalcConc.c
*/

#define USE_OPENMP 1  // Set to 1 to use OpenMP, 0 to disable

#if USE_OPENMP
#include <omp.h>
#endif

// For some compilers this is defined in math.h, for others its not
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:e_iterate:invalidNumInputs",
                          "Six input arguments required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:e_iterate:maxlhs",
                          "One output argument required.");
    }

    double *r = mxGetPr(prhs[0]);
    size_t nr = (size_t) mxGetScalar(prhs[1]);
    double D = mxGetScalar(prhs[2]);
    double dt = mxGetScalar(prhs[3]);
    double *current = mxGetPr(prhs[4]);
    size_t it = (size_t) mxGetScalar(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(1, nr, mxREAL);
    double *C = mxGetPr(plhs[0]);

    for (size_t i = 0; i < nr; i++) {
        C[i] = 0;
    }

    double four_pi_D = 4 * M_PI * D;
    double t = it * dt;

    double r_i, r_i_squared, t_prime, diff_time, G;

    #if USE_OPENMP
    #pragma omp parallel for private(r_i, r_i_squared, t_prime, diff_time, G) shared(C, r, current)
    #endif
    for (size_t i = 0; i < nr; i++) {
        r_i = r[i];
        r_i_squared = r_i * r_i;
        for (size_t j = 0; j < it; j++) {
            t_prime = j * dt;
            diff_time = t - t_prime;
            G = pow(four_pi_D * diff_time, -1.5) * exp(-r_i_squared / (four_pi_D * diff_time));
            
            #if USE_OPENMP
            #pragma omp atomic
            #endif
            C[i] += current[j] * G * dt;
        }
    }
}