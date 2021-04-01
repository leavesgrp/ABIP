#include "mex.h"
#include "abip.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    if (nrhs != 0) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    plhs[0] = mxCreateString(abip_version());
    return;
}
