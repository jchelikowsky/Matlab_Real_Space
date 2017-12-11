#include "mex.h"
/*#include "pthread.h
 
/* below are two versions of a threaded matrix vector multiple and
 * at the bottom is the non-threaded version
 * the threaded versions have NOT been tested
/* input arguments 
mxArray *sparseMatrix,*vector;
int *columnIndices,numberOfThreads;
/* output arguments
mxArray *solution;
pthread_mutex_t *mutexsOnSolution;

void *threadedMultiply(void *arguments){
    int *argPointer=(int *)arguments;
    int start=argPointer[0];
    int end=argPointer[1];
    int counter,*ir;
    double *localSolution=mxCalloc(mxGetM(sparseMatrix),sizeof(double));
    double *matrixArray,*vectorArray;
    
    ir=mxGetIr(sparseMatrix);
    matrixArray=mxGetPr(sparseMatrix);
    vectorArray=mxGetPr(vector);
    
    /* do multiplication 
    for (counter=start;counter<end;counter++){
        localSolution[ir[counter]]+=matrixArray[counter]*vectorArray[columnIndices[counter]];
    }
    
    
    /* save localSolution 
    for (counter=0;counter<mxGetM(sparseMatrix);counter++){
        pthread_mutex_lock(mutexsOnSolution+counter);
        mxGetPr(solution)[counter]+=localSolution[counter];
        pthread_mutex_unlock(mutexsOnSolution+counter);
    }
    
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){

    /* local variables 
    pthread_t *threads;
    int startingIndex=0;
    int counter,elementsPerThread;
    int *jc,errorCode;
    int *tempArray;
    
    sparseMatrix=prhs[0];
    vector=prhs[1];
    columnIndices=mxGetData(prhs[2]);
    numberOfThreads=*mxGetPr(prhs[3]);
    jc=mxGetJc(sparseMatrix);
    
    
    solution=mxCreateDoubleMatrix(mxGetM(sparseMatrix),1,mxREAL);
   
    mutexsOnSolution=mxCalloc(mxGetM(sparseMatrix),sizeof(pthread_mutex_t));
    for (counter=0;counter<mxGetM(sparseMatrix);counter++){
        pthread_mutex_init(mutexsOnSolution+counter,NULL);
    }
    
    threads=mxCalloc(numberOfThreads,sizeof(pthread_t));
    elementsPerThread=jc[mxGetN(sparseMatrix)]/numberOfThreads;
    for (counter=0;counter<numberOfThreads-1;counter++){
        tempArray=mxCalloc(2,sizeof(int));
        tempArray[0]=startingIndex;
        startingIndex+=elementsPerThread;
        tempArray[1]=startingIndex;
        errorCode=pthread_create(threads+counter,NULL,threadedMultiply,(void *)tempArray);
        if (errorCode!=0){
            mexErrMsgTxt("FAILURE TO CREATE THREAD!");
        }
    }
    tempArray=mxCalloc(2,sizeof(int));
    tempArray[0]=startingIndex;
    tempArray[1]=jc[mxGetN(sparseMatrix)];
    errorCode=pthread_create(threads+counter,NULL,threadedMultiply,(void *)tempArray);
    /*
    

    /* wait for all threads to return 
    for (counter=0;counter<numberOfThreads;counter++){
        pthread_join(threads[counter],NULL);
    }
    
    /* clean up mutexs 
    for (counter=0;counter<mxGetM(sparseMatrix);counter++){
        pthread_mutex_destroy(mutexsOnSolution+counter);
    }
    
    plhs[0]=solution;

    return;
}*/





void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* single threaded version */
    /* only slightly faster than MATLAB's * operator */
    
    /* Input arguments */
    mxArray *sparseMatrix,*vector;
    /* output arguments */
    mxArray *solution;
    /* local variables */
    mwIndex *ir,*jc;
    double *matrixArray,*vectorArray,*solutionArray;
    double *matrixArrayPermanent,*onePastLastColumnInVector,*firstItemInNextColumn;
    size_t currentColumn;
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=3){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }

    /* initialize variables */
    sparseMatrix=prhs[0];
    vector=prhs[1];
    currentColumn=*mxGetPr(prhs[2]);

    matrixArray=mxGetPr(sparseMatrix);
    matrixArrayPermanent=matrixArray;
    vectorArray=mxGetPr(vector);
    ir=mxGetIr(sparseMatrix);
    jc=mxGetJc(sparseMatrix);
    
    solution=mxCreateDoubleMatrix(mxGetM(sparseMatrix),1,mxREAL);
    solutionArray=mxGetPr(solution); 
    
    /* setup needed for forloop */
    onePastLastColumnInVector=vectorArray+mxGetN(sparseMatrix);
    jc+=currentColumn+1;

    /* outter for loop goes through each column of the matrix and
     * the inner for loop goes through each element in the current column
     */
    for (vectorArray+=currentColumn;vectorArray<onePastLastColumnInVector;vectorArray++){
        firstItemInNextColumn=matrixArrayPermanent+(*jc);
        for (;matrixArray<firstItemInNextColumn;matrixArray++){
            solutionArray[*ir]+=(*matrixArray)*(*vectorArray);
        
            ir++;
        }
        jc++;
    }
    
    plhs[0]=solution;
    
    return;
}
