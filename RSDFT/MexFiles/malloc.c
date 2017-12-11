#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    mxArray *allocatedMatrix=mxCreateDoubleMatrix(0,0,mxREAL);
    size_t rows,columns;
    double *allocatedMemory;
    
    
    if (nrhs!=2){
        mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    rows=*mxGetPr(prhs[0]);
    columns=*mxGetPr(prhs[1]);
    
    allocatedMemory=mxMalloc(rows*columns*sizeof(double));
    
    mxSetPr(allocatedMatrix,allocatedMemory);
    mxSetM(allocatedMatrix,rows);
    mxSetN(allocatedMatrix,columns);
     
    plhs[0]=allocatedMatrix;
    
    return;
}
