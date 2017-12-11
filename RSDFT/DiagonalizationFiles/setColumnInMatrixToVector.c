#include "mex.h"

/* number of rows in matrix should equal number of elements in vector */
void setColumnInMatrixToVector(mxArray *matrix,mxArray *vector, size_t column){
	size_t numberOfRows=mxGetM(matrix);
    size_t columnTimesNumberOfRows=column*numberOfRows;
    double *matrixPointer=mxGetPr(matrix);
    double *vectorPointer=mxGetPr(vector);
    double *onePastEndOfvector=vectorPointer+numberOfRows;
    
    matrixPointer+=columnTimesNumberOfRows;
	for (;vectorPointer<onePastEndOfvector;){
		*matrixPointer=*vectorPointer;
        vectorPointer++;
        matrixPointer++;
	}
}


