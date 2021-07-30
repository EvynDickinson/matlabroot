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

    char szError[256] = "";
    UINT uModel = 0; // Identifies MultiClamp 700A or 700B model
    char szSerialNum[16] = ""; // Serial number of MultiClamp 700B
    UINT uCOMPortID = 0; // COM port ID of MultiClamp 700A (1-16)
    UINT uDeviceID = 0; // Device ID of MultiClamp 700A (1-8)
    UINT uChannelID = 0; // Headstage channel ID
    
    // find the first MultiClamp
    MCCMSG_FindFirstMultiClamp(hMCCmsg, &uModel, szSerialNum,
            sizeof(szSerialNum), &uCOMPortID,
            &uDeviceID, &uChannelID, &nError); 

    // find the next MultiClamp
    MCCMSG_FindNextMultiClamp(hMCCmsg, &uModel, szSerialNum,
            sizeof(szSerialNum), &uCOMPortID,
            &uDeviceID, &uChannelID, &nError); 

    // select the first MultiClamp
    MCCMSG_SelectMultiClamp(hMCCmsg, uModel, szSerialNum,
            uCOMPortID, uDeviceID, uChannelID, &nError); 
        
    // Primary Gain
    double d1Gain = 0;
    MCCMSG_GetPrimarySignalGain(hMCCmsg, &d1Gain, &nError);
    
    mxArray *out0 = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    *((double *)mxGetData(out0)) = (double)d1Gain;
    plhs[0] = out0;

    // Primary Signal
    UINT u1SignalID = 0;
    MCCMSG_GetPrimarySignal(hMCCmsg, &u1SignalID, &nError);

    mxArray *out1 = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    *((uint32_t *)mxGetData(out1)) = (uint32_t)u1SignalID;
    plhs[1] = out1;

    // Secondary Gain
    double d2Gain = 0;
    MCCMSG_GetSecondarySignalGain(hMCCmsg, &d2Gain, &nError);

    mxArray *out2 = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    *((double *)mxGetData(out2)) = (double)d2Gain;
    plhs[2] = out2;
    
    // Secondary Signal
    UINT u2SignalID = 0;
    MCCMSG_GetSecondarySignal(hMCCmsg, &u2SignalID, &nError);

    mxArray *out3 = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    *((uint32_t *)mxGetData(out3)) = (uint32_t)u2SignalID;
    plhs[3] = out3;

    MCCMSG_DestroyObject(hMCCmsg);
    hMCCmsg = NULL;

}