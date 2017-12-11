#include "mex.h"
#include "pthread.h"

/* another attempt at threaded matrix vector multiply
 * not used in current version of RSDFT
/*pthread_mutex_t mutexOnSolution;
mxArray *solution;

void *threadedMultiply(void *args){
    mxArray **internalArgs=(mxArray **)args;
    mxArray *input[2];
    mxArray *output[2];
    
    input[0]=internalArgs[0];
    input[1]=internalArgs[1];
    /*mexCallMATLAB(1,output,2,input,"*");*/
    
    
    /*save localSolution *
    pthread_mutex_lock(&mutexOnSolution);
        input[0]=output[0];
        input[1]=solution;
        mexCallMATLAB(1,output,2,input,"+");
        solution=output[0];
    pthread_mutex_unlock(&mutexOnSolution);*
    return;
}    

int *shiftJcArray(int *jc,int numberOfColumns){
    int *newJc=mxMalloc((numberOfColumns+1)*sizeof(int));
    int counter;
    
    for (counter=0;counter<numberOfColumns+1;counter++){
        newJc[counter]=jc[counter]-jc[0];
    }
    return newJc;
}
    
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    
    /* Input arguments *
    mxArray *sparseMatrix,*vector;
    int currentColumn,numberOfThreads;
    /* output arguments *
    /* local variables *
    int counter,columnsPerThread;
    int numberOfElements,errorFlag;
    int numberOfRows;
    int *ir,*jc;
    mxArray *subMatrix,*subVector;
    double *matrixArray,*vectorArray,*solutionArray;
    double *matrixArrayPermanent,*onePastLastColumnInVector,*firstItemInNextColumn;;
    mxArray **argPasser;
    
    pthread_attr_t attr;
    pthread_t *threads;
   
    
    /*checks on inputs and number of outputs*
    if (nrhs!=4){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    /* initialize variables *
    sparseMatrix=prhs[0];
    vector=prhs[1];
    currentColumn=*mxGetPr(prhs[2]);
    numberOfThreads=*mxGetPr(prhs[3]);
    
    matrixArray=mxGetPr(sparseMatrix);
    matrixArrayPermanent=matrixArray;
    vectorArray=mxGetPr(vector);
    ir=mxGetIr(sparseMatrix);
    jc=mxGetJc(sparseMatrix);
    numberOfRows=mxGetM(sparseMatrix);
    
    solution=mxCreateDoubleMatrix(numberOfRows,1,mxREAL);
    solutionArray=mxGetPr(solution); 
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    threads=malloc(numberOfThreads*sizeof(pthread_t));
    pthread_mutex_init(&mutexOnSolution, NULL);
    
    columnsPerThread=mxGetN(sparseMatrix)/numberOfThreads;
    
    
    for (counter=0;counter<numberOfThreads-1;counter++){
       subMatrix=mxCreateDoubleMatrix(0,0,mxREAL);
       mxSetPr(subMatrix,matrixArray+jc[currentColumn]);
       mxSetM(subMatrix,numberOfRows);
       mxSetN(subMatrix,columnsPerThread);
       mxSetIr(subMatrix,ir+jc[currentColumn]);
       mxSetJc(subMatrix,shiftJcArray(jc+currentColumn,columnsPerThread));
       
       subVector=mxCreateDoubleMatrix(0,0,mxREAL);
       mxSetPr(subVector,vectorArray);
       mxSetM(subVector,columnsPerThread);
       mxSetN(subVector,1);
       vectorArray+=columnsPerThread;
       
       argPasser=malloc(3*sizeof(mxArray*));
       argPasser[0]=subMatrix;
       argPasser[1]=subVector;
       argPasser[2]=solution;
       
       pthread_create(&threads[counter],&attr,&threadedMultiply,(void *) argPasser);
       
       currentColumn+=columnsPerThread;
    }
    
    subMatrix=mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(subMatrix,matrixArray+jc[currentColumn]);
    mxSetM(subMatrix,numberOfRows);
    mxSetN(subMatrix,mxGetN(sparseMatrix)-currentColumn);
    mxSetIr(subMatrix,ir+jc[currentColumn]);
    mxSetJc(subMatrix,shiftJcArray(jc+currentColumn,mxGetN(sparseMatrix)-currentColumn));
    
    subVector=mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(subVector,vectorArray);
    mxSetM(subVector,mxGetN(sparseMatrix)-currentColumn);
    mxSetN(subVector,1);
       
    argPasser=malloc(3*sizeof(mxArray*));
    argPasser[0]=subMatrix;
    argPasser[1]=subVector;
    argPasser[2]=solution;
       
    pthread_create(&threads[numberOfThreads-1],&attr,&threadedMultiply,(void *) argPasser);
    
    /*for (counter=0;counter<numberOfThreads;counter++){
        pthread_join(threads[counter], NULL);
    }*
    
    pthread_mutex_destroy(&mutexOnSolution);
    pthread_attr_destroy(&attr);
    
    plhs[0]=solution;
    
    return;
}*/


/*
/* input arguments */
mxArray *sparseMatrix,*vector;
int *columnIndices,numberOfThreads;
/* output arguments */
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
    
    /* do multiplication */
    for (counter=start;counter<end;counter++){
        localSolution[ir[counter]]+=matrixArray[counter]*vectorArray[columnIndices[counter]];
    }
    
    
    /* save localSolution */
    for (counter=0;counter<mxGetM(sparseMatrix);counter++){
        pthread_mutex_lock(mutexsOnSolution+counter);
        mxGetPr(solution)[counter]+=localSolution[counter];
        pthread_mutex_unlock(mutexsOnSolution+counter);
    }
    
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){

    /* local variables */
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
    

    /* wait for all threads to return */
    for (counter=0;counter<numberOfThreads;counter++){
        pthread_join(threads[counter],NULL);
    }
    
    /* clean up mutexs */
    for (counter=0;counter<mxGetM(sparseMatrix);counter++){
        pthread_mutex_destroy(mutexsOnSolution+counter);
    }
    
    plhs[0]=solution;

    return;
}
