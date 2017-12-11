#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* A,x,c,sigma1,e */
    /* evaluates the equation (A*x - c*x) .* (sigma1/e);*/
    
    double *xArray;
    mxArray *solution=mxCreateDoubleMatrix(mxGetM(prhs[1]),mxGetN(prhs[1]),mxREAL);    
    double *solutionArray,*onePastEndOfArray;
    mxArray *ATimesx,*inputVariables[2];
    double *ATimesxArray;
    
    size_t numberOfElements;
    double c,sigma1,e,sigma1DividedBye;
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=5){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    /* initialize variables */
    xArray=mxGetPr(prhs[1]);
    
    c=*mxGetPr(prhs[2]);
    sigma1=*mxGetPr(prhs[3]);
    e=*mxGetPr(prhs[4]);
    
    inputVariables[0]=prhs[0];
    inputVariables[1]=prhs[1];
    mexCallMATLAB(1,&ATimesx,2,inputVariables,"*");
    ATimesxArray=mxGetPr(ATimesx);
    
    solutionArray=mxGetPr(solution);
    
    /* precalc values */
    numberOfElements=mxGetM(prhs[1])*mxGetN(prhs[1]);
    sigma1DividedBye=sigma1/e;
    
    onePastEndOfArray=solutionArray+numberOfElements;

       
    for (;solutionArray<onePastEndOfArray;solutionArray++){
        *solutionArray=(*ATimesxArray-c*(*xArray))*sigma1DividedBye;
        ATimesxArray++;
        xArray++;
    }
    
    mxDestroyArray(ATimesx);
   
    plhs[0]=solution;

    return;
}
