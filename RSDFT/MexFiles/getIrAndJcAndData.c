#include "mex.h"
#include "string.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* takes in a sparse matrix and returns its ir and jc array and the data array */
    mxArray *sparseMatrix,*solution1,*solution2,*solution3;
    int *ir,*jc,*tempPointer;
    int *solutionArray1,*solutionArray2;
    int counter;
    double *data,*solutionArray3;   

    sparseMatrix=prhs[0];
    ir=mxGetIr(sparseMatrix);
    jc=mxGetJc(sparseMatrix);
    data=mxGetPr(sparseMatrix);

    solution1=mxCreateNumericMatrix(1,jc[mxGetN(sparseMatrix)],mxINT32_CLASS,mxREAL);
    solution2=mxCreateNumericMatrix(1,mxGetN(sparseMatrix)+1,mxINT32_CLASS,mxREAL);
    solution3=mxCreateDoubleMatrix(1,jc[mxGetN(sparseMatrix)],mxREAL); 

    solutionArray1=mxGetData(solution1);    
    solutionArray2=mxGetData(solution2);
    solutionArray3=mxGetData(solution3);

    for (counter=0;counter<jc[mxGetN(sparseMatrix)];counter++){
        solutionArray1[counter]=ir[counter];
        solutionArray3[counter]=data[counter];
    }
    for (counter=0;counter<mxGetN(sparseMatrix)+1;counter++){
        solutionArray2[counter]=jc[counter];
    }
  
    
    plhs[0]=solution1;
    plhs[1]=solution2;
    plhs[2]=solution3;    

    return;
}

