#include "mex.h"
#include "matrixAndVectorOps.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* takes in 5 arguments and returns 1 mxArray 
     * function evaluates the equation B*v1-bet*v0
     */
    /* input arguments */
    /* B is a sparse matrix. v1,v0 are vectors */
    mxArray *B,*v1,*v0;
    double bet;
    size_t column;/* column should always be >= to 0*/
    
    /* output argument, solution is a vector */
    mxArray *solution;
    
    /* local variables */
    double *solutionArray,*solutionArrayPermanent,*v0Array,*v1Array,*onePastEndOfv0Array;
    double *matrixArray,*matrixArrayPermanent,*firstItemInNextColumn;
    mwIndex *jc,*ir;
    
     /*checks on inputs and number of outputs*/
    if (nrhs!=5){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    /* take arguments and transfer them to local variables */
    B=prhs[0];
    matrixArray=mxGetPr(B);
    matrixArrayPermanent=matrixArray;
    ir=mxGetIr(B);
    jc=mxGetJc(B);
    
    
    v0=prhs[1];
    v0Array=mxGetPr(v0);
    
    v1=prhs[2];
    v1Array=mxGetPr(v1);
    
    bet=*mxGetPr(prhs[3]);
    column=*mxGetPr(prhs[4]);
    
    /* allocate solution */
    solution=mxCreateDoubleMatrix(mxGetM(v0),1,mxREAL);
    solutionArray=mxGetPr(solution);
    solutionArrayPermanent=solutionArray;
    
    /* find ending criteria and increment jc so it points to the 
     * second element in the array.  See MATLAB help to understand 
     * how sparse matrices are stored
     */
    onePastEndOfv0Array=v0Array+mxGetM(v0);
    jc+=column+1;
    
    /* outer for loop iterates through each column of the matrix
     * the inner for loop iterates through all the elements in the
     * current column of the matrix
     * loops have been converted as much as possible to pointer code
     */
    for (;v0Array<onePastEndOfv0Array;v0Array++){
        
        firstItemInNextColumn=matrixArrayPermanent+(*jc);
        for (;matrixArray<firstItemInNextColumn;matrixArray++){
            /* solutionArrayPermanent[*ir] finds the appropriate row
             * in the solution array based on the row that the element in
             * the matrixArray is in
             */
            solutionArrayPermanent[*ir]+=(*matrixArray)*(*v1Array);
        
            ir++;
        }
        jc++;
        v1Array++;
        
        *solutionArray-=bet*(*v0Array);
        solutionArray++;
    }
    
    
    plhs[0]=solution;
    
    return;
}


