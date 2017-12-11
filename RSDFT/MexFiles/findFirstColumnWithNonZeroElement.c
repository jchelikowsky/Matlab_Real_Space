#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
   /* takes in a sparse matrix and finds which colum has the first
    * non-zero element
    */
    
   mwIndex *jc;
   size_t columns;
   size_t counter;
   size_t index=0;
   
   /*checks on inputs and number of outputs*/
   if (nrhs!=1){
      mexErrMsgTxt("Incorrect number of input arguments"); 
   }
   if (nlhs!=1){
       mexErrMsgTxt("Incorrect number of output arguments"); 
   }
   
   jc=mxGetJc(prhs[0]);
   columns=mxGetN(prhs[0]); 
    
   /* find first column that has an element*/
   for (counter=0;counter<columns;counter++){
        if (jc[counter+1]-jc[counter]>0){
            index=counter;
            break;
        }
    }
   
   plhs[0]=mxCreateDoubleScalar(index);
   
   return;
}
