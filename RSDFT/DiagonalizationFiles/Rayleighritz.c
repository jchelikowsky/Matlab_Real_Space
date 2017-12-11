#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
	/* The name of this function may not be correct.  used the name from a nearby comment 
	   in the .m file
	*/
	/*  input:
		Vin,W,n
		output:
		G
	*/
	/* input variables*/
	mxArray *Vin,*W;
	size_t n;
	/* ouput variables*/
	mxArray *G;
	
	/* local variables*/
	size_t i,j;
	double *GArray;
	double *VinArray,*WArray;
	size_t numberOfRowsInVin;
	size_t jTimesn;
	
    double sum=0;
    size_t counter;
    size_t iTimesRows;
    size_t jTimesRows;
    
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=3){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=1){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    
	Vin=prhs[0];
	W=prhs[1];
	n=*mxGetPr(prhs[2]);
	
	G=mxCreateDoubleMatrix(n,n,mxREAL);
	GArray=mxGetPr(G);
	
	VinArray=mxGetPr(Vin);
	WArray=mxGetPr(W);
	
	
	numberOfRowsInVin=mxGetM(Vin);
	
    /* MATLAB version
     *for j=1:mev
        	for i=1:j
            	G(i,j) = Vin(:,i)'*W(:,j);
            	G(j,i) = G(i,j);
        	end
    	end
     */
    
	for (j=0;j<n;j++){
		jTimesn=j*n;
        jTimesRows=j*numberOfRowsInVin;
		for (i=0;i<=j;i++){
			iTimesRows=i*numberOfRowsInVin;
            sum=0;
            for (counter=0;counter<numberOfRowsInVin;counter++){
                sum+=VinArray[iTimesRows+counter]*WArray[jTimesRows+counter];
            }
            
            GArray[jTimesn+i]=sum;
			GArray[i*n+j]=sum;
		}
	}
	
	plhs[0]=G;
	
	return;
}
