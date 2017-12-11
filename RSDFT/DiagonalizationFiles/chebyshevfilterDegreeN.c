#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* A,x,y,c,t1,t2 */
    
    
    double *xArray,*yArray;
    mxArray *solution=mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);    
    double *solutionArray,*onePastEnd;
    mxArray *ATimesy,*inputVariables[2];
    double *ATimesyArray;
    
    size_t numberOfElements;
    double c,t1,t2;
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=6){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    xArray=mxGetPr(prhs[1]);
    yArray=mxGetPr(prhs[2]);
    
    c=*mxGetPr(prhs[3]);
    t1=*mxGetPr(prhs[4]);
    t2=*mxGetPr(prhs[5]);
    
    inputVariables[0]=prhs[0];
    inputVariables[1]=prhs[2];
    mexCallMATLAB(1,&ATimesy,2,inputVariables,"*");
    ATimesyArray=mxGetPr(ATimesy);
    
    solutionArray=mxGetPr(solution);
    
    numberOfElements=mxGetM(prhs[1])*mxGetN(prhs[1]);
    onePastEnd=solutionArray+numberOfElements;
    
    for (;solutionArray<onePastEnd;solutionArray++){
        *solutionArray=(*ATimesyArray-c*(*yArray))*t1-t2*(*xArray);
        
        ATimesyArray++;
        yArray++;
        xArray++;
    }
    
    plhs[0]=solution;

    return;
}
