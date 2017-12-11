#include "mex.h"
#include "math.h"

/* TODO
// this function can be inlined*/
int convertFromIndicesToIndex(int i, int j, int rows, int columns){
	/* returns column major order index*/
	return j*rows+i;
}

bool compare(double *v1,double *v2, int numberOfElements){
    int counter;   

    for (counter=0;counter<numberOfElements;counter++){
        if (fabs(v1[counter]-v2[counter])>0.0000000001){
            return false;
        }
    }
    return true;
}




/*double norm(mxArray *vector){
    int length=mxGetM(vector)*mxGetN(vector);
    int counter;
    double sum=0;
    double *pointer=mxGetPr(vector);
    
    for (counter=0;counter<length;counter++){
        sum+=pointer[counter]*pointer[counter];
    }
    return sqrt(sum);
}*/
double norm(double *vector, int length){
    int counter;
    double sum=0;
    
    for (counter=0;counter<length;counter++){
        sum+=vector[counter]*vector[counter];
    }
    return sqrt(sum);
}


double absMax(double* vector, int numberOfElements){
	double max=vector[0];
	int counter;
	
	if (max<0.0){
		max=max*-1.0;
	}
	
	for (counter=1;counter<numberOfElements;counter++){
		if (vector[counter]>max){
			max=vector[counter];
		}
		else if (-1.0*vector[counter]>max){
			max=-1.0*vector[counter];
		}
	}
	return max;
}

mxArray* getDiagonal(double *matrix,int rows){
	double *diag=mxCalloc(rows,sizeof(double));
	int counter;
    mxArray *solution=mxCreateDoubleMatrix(0,0,mxREAL);
    
	for (counter=0;counter<rows;counter++){
        diag[counter]=matrix[counter*rows+counter];
	}

    mxSetPr(solution,diag);
    mxSetM(solution,rows);
    mxSetN(solution,1);
    
	return solution;
}

/* end is inclusive */
double sumVector(double* vector, int start , int end){
	int counter;
	double sum=0;
	
	for (counter=start;counter<=end;counter++){
		sum+=vector[counter];
	}	
	return sum;
}


double* subtractMatrices(double* matrix1, double* matrix2, int rows, int columns){

	int i;
	int numberOfElements=rows*columns;
	double* solution=mxCalloc(numberOfElements,sizeof(double));
	double temp=0;

	/*// new method takes the most advantage of cache behavior
	// TODO
	// loop can be unrolled*/
	for (i=0;i<numberOfElements;i++){
        temp=matrix1[i]-matrix2[i];
        solution[i]=temp;
	}

	return solution;
}

double* multiplyMatrixByScaler(double* matrix, int rows, int columns, double scaler){
	int i;
	int numberOfElements=rows*columns;
	double* solution=mxCalloc(numberOfElements,sizeof(double));
    
	/*// new method takes the most advantage of cache behavior
	// TODO
	// loop can be unrolled*/
	for (i=0;i<numberOfElements;i++){
        solution[i]=scaler*matrix[i];
	}

	return solution;
}



/* starts at 0 and goes to row/column inclusive*/
mxArray* getSubMatrix(double *matrix, size_t row, size_t column, size_t numberOfRows){
	double *submatrix=mxCalloc((row+1)*(column+1),sizeof(double));
	size_t i,j;
	mxArray *solution=mxCreateDoubleMatrix(0,0,mxREAL);
    size_t rowPlusOne=row+1,jTimesRowPlusOne,jTimesNumberOfRows;
    
	for (j=0;j<=column;j++){
		jTimesRowPlusOne=j*rowPlusOne;
        jTimesNumberOfRows=j*numberOfRows;
        for (i=0;i<=row;i++){
			submatrix[i+jTimesRowPlusOne]=matrix[i+jTimesNumberOfRows];
		}
	}
    
    mxSetPr(solution,submatrix);
    mxSetM(solution,rowPlusOne);
    mxSetN(solution,column+1);
    
	return solution;
} 
