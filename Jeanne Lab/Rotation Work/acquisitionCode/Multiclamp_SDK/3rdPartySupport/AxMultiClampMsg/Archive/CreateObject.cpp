#include <windows.h>
#include <stdint.h>
#include "mex.h"
#include "AxMultiClampMsg.h"

void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )

{
    // create DLL handle
    int nError;
    nError = 0;
    
    HMCCMSG hMCCmsg;
    hMCCmsg = MCCMSG_CreateObject(&nError);
    
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *((uint32_t *)mxGetData(out)) = (uint32_t)hMCCmsg;
    
    plhs[0] = out;
}