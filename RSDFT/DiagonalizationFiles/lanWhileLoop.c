#include "mex.h"
#include "matrixAndVectorOps.c"
#include "setColumnInMatrixToVector.c"
#include "math.h"
#include "time.h"

double min(double a, double b){
	if (a<b){
		return a;
	}
	else {
		return b;
	}
}
double max(double a, double b){
	if (a>b){
		return a;
	}
	else {
		return b;
	}
}


long *cycles;
long *cyclesStart;
int *numberOfExecutions;
int numberOfItems=0;

void setupProfiling(int numberOfItemsToProfile){
    cycles=mxCalloc(numberOfItemsToProfile,sizeof(long));
    numberOfExecutions=mxCalloc(numberOfItemsToProfile,sizeof(int));
    
    cyclesStart=mxCalloc(numberOfItemsToProfile,sizeof(long));
    numberOfItems=numberOfItemsToProfile;
}

void ticTimer(int item){
    cyclesStart[item]=clock();
}
void tocTimer(int item){
    cycles[item]+=clock()-cyclesStart[item];
    numberOfExecutions[item]++;
}

void printStats(){
    int counter;
    
    for (counter=0;counter<numberOfItems;counter++){
        mexPrintf("item %d took %f seconds\n",counter,((double)cycles[counter])/CLOCKS_PER_SEC);
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
	/*input
		m,v0,v1,B,reorth,nev,VV,Tmat,tr,tol,enableConvergenceTest,percentageComplete,handle
	*/
	/* output
		rr,X1,indx,k
     *
     *VV is directly changed, does not need to be returned
	*/
	/* input variables */
	size_t m,reorth,enableConvergenceTest;
	double nev;
	double tolValue,percentageComplete;
    mxArray *v0,*v1,*handle,*barHandle;
	mxArray *B;
	mxArray *VV,*Tmat,*tr,*tol;
	/* output variables */
	mxArray *rr,*X1,*indx;
	/* local variables */
	size_t k,NTest,ll, numberOfRows,numberOfColumns;
	double alp,bet,tr1,inverseOfbet;
	double *trPointer;
	mxArray *v=mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *inputVariables[5],*outputVariables[5];
	double *TmatPointer;
    mxArray *columnArray;
    bool enableProgressBar=true;
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=14 && nrhs!=11){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nrhs==11){
        enableProgressBar=false;
    }    
    if (nlhs!=4){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
	/* initialize variables */
    k=0;
	bet=0;
	ll=0;
	
	m=*mxGetPr(prhs[0]);
	v0=prhs[1];
	v1=prhs[2];
	B=prhs[3];
	reorth=*mxGetPr(prhs[4]);
	nev=*mxGetPr(prhs[5]);
	VV=prhs[6];
	Tmat=prhs[7];
	tr=prhs[8];
	tol=prhs[9];
    enableConvergenceTest=*mxGetPr(prhs[10]);
    if (enableProgressBar){
        percentageComplete=*mxGetPr(prhs[11]);
        handle=prhs[12];
        barHandle=prhs[13];
    }
    
	TmatPointer=mxGetPr(Tmat);
	trPointer=mxGetPr(tr);
	
    tolValue=*mxGetPr(tol);
    
    mexCallMATLAB(1,&columnArray,1,&B,"findFirstColumnWithNonZeroElement");

    
    
	while (k < m){
    	k++;
    	
        /* VV(:,k) = v1;*/ 
        /* directly changes VV */
        setColumnInMatrixToVector(VV,v1,k-1);
        /* VV(:,k) = v1; */
        
        
        /*v   =  B*v1 - bet*v0;*/
        inputVariables[0]=B;
        inputVariables[1]=v0;
        inputVariables[2]=v1;
        inputVariables[3]=mxCreateDoubleScalar(bet);
        inputVariables[4]=columnArray;
        mexCallMATLAB(1,outputVariables,5,inputVariables,"BTimesv1MinusbetTimesv0");
        mxDestroyArray(v);
        v=outputVariables[0];
        /*v   =  B*v1 - bet*v0; */
        
        
        
        
        /*  alp = v1'*v;*/
        /* transposes v1 */
        numberOfRows=mxGetM(v1);
        numberOfColumns=mxGetN(v1);
        mxSetM(v1,numberOfColumns);
        mxSetN(v1,numberOfRows);
        
        inputVariables[0]=v1;
        inputVariables[1]=v;
        mexCallMATLAB(1,outputVariables,2,inputVariables,"*");
        alp=*mxGetPr(outputVariables[0]);
        
        /* undoes the transpose  */
        mxSetM(v1,numberOfRows);
        mxSetN(v1,numberOfColumns);
        /*  alp = v1'*v; */
        
        
        
        /* v   = v-alp*v1 ;*/
        inputVariables[0]=v;
        inputVariables[1]=v1;
        inputVariables[2]=mxCreateDoubleScalar(alp);
        mexCallMATLAB(1,outputVariables,3,inputVariables,"vMinusalpTimesv1");
        mxDestroyArray(v);
        v=outputVariables[0];
        mxDestroyArray(inputVariables[2]);
        /* v   = v-alp*v1 ; */
        
        
        /* %%-----------------------------
        %% reorth  -- test for reorth. needed!
        %%-----------------------------*/
        /* if reorth is going to be used on a regular basis, this code
         * must be sped up.  This code is only here for compatability
         */
        if (reorth==1) {
            mexPutVariable("caller","VVCopy",VV);
            mexPutVariable("caller","vCopy",v);
            mexPutVariable("caller","kCopy",mxCreateDoubleScalar(k));
            mexEvalString("newv=vCopy-VV(:,1:kCopy)*(VV(:,1:kCopy)'*vCopy);");
            v=mexGetVariable("caller","newv");
        }
        /*%%-------------------- normalize and store v1*/
        
        /* bet = norm(v);*/
        inputVariables[0]=v;
        mexCallMATLAB(1,outputVariables,1,inputVariables,"norm");
        bet=*mxGetPr(outputVariables[0]);
        /* bet = norm(v);   */ 
        
        if (k>2){
            /* destroying v0 on the first or second iteration of the loop
             * will deallocate a variable in the MATLAB workspace which can
             * cause MATLAB to crash.
             */
            mxDestroyArray(v0);
        }
        v0=v1;
               
        /* v1 = v/bet;*/
        inverseOfbet=1.0/bet;
        inputVariables[0]=v;
        inputVariables[1]=mxCreateDoubleScalar(inverseOfbet);
        mexCallMATLAB(1,outputVariables,2,inputVariables,"*");
        v1=outputVariables[0];
        
        mxDestroyArray(inputVariables[1]);
        /* v1 = v/bet; */
        
        
        /* %%-------------------- update tmat */
        TmatPointer[(k-1)+(m+1)*(k-1)]   = alp;
        TmatPointer[k+(m+1)*(k-1)] = bet;
        TmatPointer[(k-1)+(m+1)*k] = bet;
        NTest  = min(8*nev,m) ;     /*%% when to   start testing
        %%-------------------- tr, ll,  == for plotting*/
        
        if (enableConvergenceTest && (((k >= NTest) && (k%10 == 0 )) || k == m)){
            /* Tmat(1:k,1:k) */
            inputVariables[0]=getSubMatrix(mxGetPr(Tmat),k-1,k-1,mxGetM(Tmat));
            /* Tmat(1:k,1:k) */
             
             /* must have gone into this if statement atleast once for indx
              * to not be null
              */
            if (ll>0){
                mxDestroyArray(indx);
            }
            
            if (k!=m){
                /* rr  = eig(Tmat(1:k,1:k));*/
                mexCallMATLAB(1,outputVariables,1,inputVariables,"eig");
                rr=outputVariables[0];
                mxDestroyArray(inputVariables[0]);
                /* rr  = eig(Tmat(1:k,1:k)); */
                
                
                /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
                inputVariables[0]=rr;
                mexCallMATLAB(2,outputVariables,1,inputVariables,"sort");
                mxDestroyArray(rr);
                rr=outputVariables[0];
                indx=outputVariables[1];
                /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
            }
            else {
                /* [X1,rr]  = eig(Tmat(1:k,1:k)); */
                mexCallMATLAB(2,outputVariables,1,inputVariables,"eig");
                X1=outputVariables[0];
                rr=outputVariables[1];
                mxDestroyArray(inputVariables[0]);
                /* [X1,rr]  = eig(Tmat(1:k,1:k));*/ 
                
                /* rr=diag(rr,0); */
                inputVariables[0]=rr;
                inputVariables[1]=mxCreateDoubleScalar(0);
                mexCallMATLAB(1,outputVariables,2,inputVariables,"diag");
                mxDestroyArray(rr);
                mxDestroyArray(inputVariables[1]);
                rr=outputVariables[0];
                /* rr=diag(rr,0); */
                
                /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
                inputVariables[0]=rr;
                mexCallMATLAB(2,outputVariables,1,inputVariables,"sort");
                mxDestroyArray(rr);
                rr=outputVariables[0];
                indx=outputVariables[1];
                /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
            }
            /*tr1 = sum(rr(1:nev));*/
            tr1=sumVector(mxGetPr(rr),0,nev-1);
            /*tr1 = sum(rr(1:nev));*/
            
            ll++;
            trPointer[ll-1] = tr1;
        }
        
        if (enableConvergenceTest && (ll>1 && (fabs(trPointer[ll-1]-trPointer[ll-2]) < tolValue*trPointer[ll-2]))){
            mxDestroyArray(indx);
            
            /* [X1,rr]  = eig(Tmat(1:k,1:k)); */
            inputVariables[0]=getSubMatrix(mxGetPr(Tmat),k-1,k-1,mxGetM(Tmat));
            mexCallMATLAB(2,outputVariables,1,inputVariables,"eig");
            X1=outputVariables[0];
            rr=outputVariables[1];
            mxDestroyArray(inputVariables[0]);
            /* [X1,rr]  = eig(Tmat(1:k,1:k));*/
            
            /* rr=diag(rr,0); */
            inputVariables[0]=rr;
            inputVariables[1]=mxCreateDoubleScalar(0);
            mexCallMATLAB(1,outputVariables,2,inputVariables,"diag");
            mxDestroyArray(rr);
            mxDestroyArray(inputVariables[1]);
            rr=outputVariables[0];
            /* rr=diag(rr,0); */
            
            /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
            inputVariables[0]=rr;
            mexCallMATLAB(2,outputVariables,1,inputVariables,"sort");
            mxDestroyArray(rr);
            rr=outputVariables[0];
            indx=outputVariables[1];
            /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
            break;
        }
        
        if (enableProgressBar && k%20==0){
            inputVariables[0]=mxCreateDoubleScalar(percentageComplete+.1*k/m);
            inputVariables[1]=barHandle;
            inputVariables[2]=handle;
            mexCallMATLAB(0,NULL,3,inputVariables,"updatePercentageComplete");
        }   
        
        
    }
	
    
    if (enableConvergenceTest==0){
     /* [X1,rr]  = eig(Tmat(1:k,1:k)); */
            inputVariables[0]=getSubMatrix(mxGetPr(Tmat),k-1,k-1,mxGetM(Tmat));
            mexCallMATLAB(2,outputVariables,1,inputVariables,"eig");
            X1=outputVariables[0];
            rr=outputVariables[1];
            mxDestroyArray(inputVariables[0]);
            /* [X1,rr]  = eig(Tmat(1:k,1:k));*/
            
            /* rr=diag(rr,0); */
            inputVariables[0]=rr;
            inputVariables[1]=mxCreateDoubleScalar(0);
            mexCallMATLAB(1,outputVariables,2,inputVariables,"diag");
            mxDestroyArray(rr);
            mxDestroyArray(inputVariables[1]);
            rr=outputVariables[0];
            /* rr=diag(rr,0); */
            
            /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
            inputVariables[0]=rr;
            mexCallMATLAB(2,outputVariables,1,inputVariables,"sort");
            mxDestroyArray(rr);
            rr=outputVariables[0];
            indx=outputVariables[1];
            /*[rr, indx]  = sort(rr) ;      %% sort increasingly*/
    }   
    
	plhs[0]=rr;
	plhs[1]=X1;
	plhs[2]=indx;
	plhs[3]=mxCreateDoubleScalar(k);
	
	
	return;
}
