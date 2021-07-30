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
    
    UINT uMode = 0;
    MCCMSG_GetMode(hMCCmsg, &uMode, &nError);

    MCCMSG_DestroyObject(hMCCmsg);
    hMCCmsg = NULL;

    mxArray *out = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    *((uint32_t *)mxGetData(out)) = (uint32_t)uMode;
    //mxArray *outstr = mxCreateString(uMode);
        
    plhs[0] = out;

}
