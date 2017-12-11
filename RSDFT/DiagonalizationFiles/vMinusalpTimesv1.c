#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* takes in 3 args and returns 1 */
    /* evaluates the equation v-alp*v1 */
    
    mxArray *v,*v1;
    double alp;
    
    double *vArray,*v1Array,*onePastLastElementInv;
    double *solution,*tempPointer;
    size_t numberOfElements;
    
     /*checks on inputs and number of outputs*/
    if (nrhs!=3){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    /* initializae variables*/
    v=prhs[0];
    v1=prhs[1];
    alp=*mxGetPr(prhs[2]);
    
    vArray=mxGetPr(v);
    v1Array=mxGetPr(v1);
    numberOfElements=mxGetM(v)*mxGetN(v);
    
    solution=mxCalloc(numberOfElements,sizeof(double));
    tempPointer=solution;/* save a pointer to the start of solution*/
    
    /* find the stoping condition */
    onePastLastElementInv=vArray+numberOfElements;
    
    /* iterate through each element in the arrays */
    for (;vArray<onePastLastElementInv;vArray++){
        *solution=*vArray-alp*(*v1Array);
        v1Array++;
        solution++;
    }
        
    solution=tempPointer;
    
    plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(plhs[0],solution);
    mxSetM(plhs[0],mxGetM(v));
    mxSetN(plhs[0],mxGetN(v));
    
    return;
}
