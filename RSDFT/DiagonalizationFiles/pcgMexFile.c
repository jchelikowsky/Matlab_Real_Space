#include "mex.h"
#include "string.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* function [x, its] = pcg (A, rhs, x0, m, tol) */
    mxArray *A,*rhs,*x0,*x,*z,*r,*p;
    double tol;
    
    mxArray *ap;
    double ro1,tol1,ro,alp,bet,temp,*tempPointer;
    size_t its,m;
    
    double *AArray,*AArrayPermanent,*rhsArray,*xArray,*rArray,*rArrayPermanent,*pArray,*apArray,*zArray;
    double *onePastEndOfrhsArray,*firstItemInNextColumn,*onePastEndOfArray;
    double apTransposedTimesp;
    mwIndex *ir,*jc;
    size_t rows,columns;
    size_t firstColumnWithNonzeroElement;
    
    mxArray *inputVariables[4];
    mxArray *outputVariables[2];
    
    /*checks on inputs and number of outputs*/
    if (nrhs!=5){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=2){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    A=prhs[0];
    AArray=mxGetPr(A);
    AArrayPermanent=AArray;
    ir=mxGetIr(A);
    jc=mxGetJc(A);
    
    rhs=prhs[1];
    rhsArray=mxGetPr(rhs);
    
    
    x0=prhs[2];
    m=*mxGetPr(prhs[3]);
    tol=*mxGetPr(prhs[4]);
    
    
    x = x0;
    xArray=mxGetPr(x);
    
    mexCallMATLAB(1,outputVariables,1,&A,"findFirstColumnWithNonZeroElement");
    firstColumnWithNonzeroElement=*mxGetPr(outputVariables[0]);
    
    
    r=mxCreateDoubleMatrix(mxGetM(x),mxGetN(x),mxREAL);
    rArray=mxGetPr(r);
    rArrayPermanent=rArray;
    /*r = rhs - A * x;*/
    onePastEndOfrhsArray=rhsArray+mxGetM(rhs);
    jc+=firstColumnWithNonzeroElement+1;
    
    for (;rhsArray<onePastEndOfrhsArray;rhsArray++){
        firstItemInNextColumn=AArrayPermanent+(*jc);
        for (;AArray<firstItemInNextColumn;AArray++){
            rArrayPermanent[*ir]-=(*AArray)*(*xArray);
        
            ir++;
        }
        jc++;
        xArray++;
        
        *rArray+=*rhsArray;
        rArray++;
    }
    /*r = rhs - A * x;*/
    
    /* z=r */
    z=mxCreateDoubleMatrix(0,0,mxREAL);
    tempPointer=mxMalloc(mxGetM(r)*sizeof(double));
    mxSetPr(z,memcpy(tempPointer,mxGetPr(r),mxGetM(r)*sizeof(double)));
    mxSetM(z,mxGetM(r));
    mxSetN(z,mxGetN(r));
    /*z=r;*/
    
    /*p = z ;*/
    p=mxCreateDoubleMatrix(0,0,mxREAL);
    tempPointer=mxMalloc(mxGetM(z)*sizeof(double));
    mxSetPr(p,memcpy(tempPointer,mxGetPr(z),mxGetM(z)*sizeof(double)));
    mxSetM(p,mxGetM(z));
    mxSetN(p,mxGetN(z));
    /*p = z ;*/
    
    /*ro1 = z' * r;*/
    /* transpose z */
    rows=mxGetM(z);
    columns=mxGetN(z);
    mxSetM(z,columns);
    mxSetN(z,rows);
    
    inputVariables[0]=z;
    inputVariables[1]=r;
    mexCallMATLAB(1,outputVariables,2,inputVariables,"*");
    ro1=*mxGetPr(outputVariables[0]);
    
    /* reverse the transpose */
    mxSetM(z,rows);
    mxSetN(z,columns);
    /*ro1 = z' * r;*/
    
    tol1 = tol*tol*ro1;
    
    
    its = 0 ;
    while (its < m && ro1 > tol1) {
        its++;
        
        ro = ro1;
        
        /*ap = A * p;*/
        jc=mxGetJc(A);
        ir=mxGetIr(A);
        pArray=mxGetPr(p);
        AArray=mxGetPr(A);
        AArrayPermanent=AArray;
        
        /* cannot destroy ap on first iteration because it has not been allocated yet*/
        if (its>1){
            mxDestroyArray(ap);
        }
        
        ap=mxCreateDoubleMatrix(mxGetM(A),1,mxREAL);
        apArray=mxGetPr(ap);
        
        onePastEndOfArray=pArray+mxGetN(A);
    
        jc+=firstColumnWithNonzeroElement+1;

        for (pArray+=firstColumnWithNonzeroElement;pArray<onePastEndOfArray;pArray++){
            firstItemInNextColumn=AArrayPermanent+(*jc);
            for (;AArray<firstItemInNextColumn;AArray++){
                apArray[*ir]+=(*AArray)*(*pArray);
        
                ir++;
            }
            jc++;
        }
        /*ap = A * p;*/
        
        
        /*alp = ro / ( ap'* p ) ;*/
        /* transpose ap */
        rows=mxGetM(ap);
        columns=mxGetN(ap);
        mxSetM(ap,columns);
        mxSetN(ap,rows);

        inputVariables[0]=ap;
        inputVariables[1]=p;
        mexCallMATLAB(1,outputVariables,2,inputVariables,"*");
        apTransposedTimesp=*mxGetPr(outputVariables[0]);
        alp=ro/apTransposedTimesp;
        
        /* reverse the transpose */
        mxSetM(ap,rows);
        mxSetN(ap,columns);        
        /*alp = ro / ( ap'* p ) ;*/
        
        
        
        /*x = x + alp * p ;*/
        pArray=mxGetPr(p);
        xArray=mxGetPr(x);
        onePastEndOfArray=xArray+mxGetM(x);
        for (;xArray<onePastEndOfArray;xArray++){
            *xArray+=alp*(*pArray);
            pArray++;
        }
        /*x = x + alp * p ;*/
        
        /*r = r - alp * ap;*/
        apArray=mxGetPr(ap);
        rArray=mxGetPr(r);
        onePastEndOfArray=rArray+mxGetM(r);
        for (;rArray<onePastEndOfArray;rArray++){
            *rArray-=alp*(*apArray);
            apArray++;
        }
        /*r = r - alp * ap;*/
        
        /* z=r */
        mxDestroyArray(z);
        z=mxCreateDoubleMatrix(0,0,mxREAL);
        tempPointer=mxMalloc(mxGetM(r)*sizeof(double));
        mxSetPr(z,memcpy(tempPointer,mxGetPr(r),mxGetM(r)*sizeof(double)));
        mxSetM(z,mxGetM(r));
        mxSetN(z,mxGetN(r));
        /*z=r;*/
        
        /*ro1 = z' * r;*/
        /* transpose z */
        rows=mxGetM(z);
        columns=mxGetN(z);
        mxSetM(z,columns);
        mxSetN(z,rows);

        inputVariables[0]=z;
        inputVariables[1]=r;
        mexCallMATLAB(1,outputVariables,2,inputVariables,"*");
        ro1=*mxGetPr(outputVariables[0]);

        /* reverse the transpose */
        mxSetM(z,rows);
        mxSetN(z,columns);
        /*ro1 = z' * r;*/
        
        bet = ro1 / ro ;
        
        /*p = z + bet * p;*/
        pArray=mxGetPr(p);
        zArray=mxGetPr(z);
        onePastEndOfArray=pArray+mxGetM(p);
        for (;pArray<onePastEndOfArray;pArray++){
            temp=(*zArray)+bet*(*pArray);
            *pArray=temp;
            zArray++;
        }
        /*p = z + bet * p;*/
    }

    plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    tempPointer=mxMalloc(mxGetM(x)*sizeof(double));
    mxSetPr(plhs[0],memcpy(tempPointer,mxGetPr(x),mxGetM(x)*sizeof(double)));
    mxSetM(plhs[0],mxGetM(x));
    mxSetN(plhs[0],mxGetN(x));
    
    plhs[1]=mxCreateDoubleScalar(its);
    
    return;
}

