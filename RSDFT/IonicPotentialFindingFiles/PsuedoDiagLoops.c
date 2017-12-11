#include "mex.h"

void fsplevalIO(mxArray *zArray,mxArray *cArray, mxArray *dArray, mxArray *xiArray, double x, double j_temp, double *returnArray);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    /* takes in 17 args and outputs 3 arrays*/
	/* Input Variables */
	mxArray *z_p_s,*c_p_s,*d_p_s,*x_pot_s;
	mxArray *z_chg,*c_chg,*d_chg,*x_charg;
	mxArray *z_vht,*c_vht,*d_vht,*x_vhart;
	mxArray *r1;
	mxArray *j_p_s,*j_ch,*j_vht;
	size_t nx;

	/* output variables*/
    mxArray *ppot,*rrho,*hpot00;

	/* local variables */
	size_t i,iPlusOne;
	double *returnArray=mxCalloc(2,sizeof(double));
	double *r1Pointer;
	double *ppotArray;
	double *rrhoArray;
	double *hpot00Array;
	double *j_p_sArray,*j_chArray,*j_vhtArray;/* named ...Array, but actually only have one elment in them*/

    
    /*checks on inputs and number of outputs*/
    if (nrhs!=17){
       mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=3){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
	/* initialize variables */
    z_p_s=prhs[0];
	c_p_s=prhs[1];
	d_p_s=prhs[2];
	
	z_chg=prhs[3];
	c_chg=prhs[4];
	d_chg=prhs[5];
	
	z_vht=prhs[6];
	c_vht=prhs[7];
	d_vht=prhs[8];
	
	x_pot_s=prhs[9];
	x_charg=prhs[10];
	x_vhart=prhs[11];
	
	r1=prhs[12];
	r1Pointer=mxGetPr(r1);
	
	j_p_s=prhs[13];
	j_ch=prhs[14];
	j_vht=prhs[15];
	
	j_p_sArray=mxGetPr(j_p_s);
	j_chArray=mxGetPr(j_ch);
	j_vhtArray=mxGetPr(j_vht);
	
	nx=*mxGetPr(prhs[16]);

	ppot=mxCreateDoubleMatrix(nx,1,mxREAL);
	rrho=mxCreateDoubleMatrix(nx,1,mxREAL);
	hpot00=mxCreateDoubleMatrix(nx,1,mxREAL);
	
	ppotArray=mxGetPr(ppot);
	rrhoArray=mxGetPr(rrho);
	hpot00Array=mxGetPr(hpot00);
	
	/* loop has been unrolled once because nx is ALWAYS even */
    for (i=0;i<nx;i+=2){
		fsplevalIO(z_p_s,c_p_s,d_p_s,x_pot_s,r1Pointer[i],*mxGetPr(j_p_s),returnArray);
		ppotArray[i]=returnArray[0];
		*j_p_sArray=returnArray[1];
		
		fsplevalIO(z_chg,c_chg,d_chg,x_charg,r1Pointer[i],*mxGetPr(j_ch),returnArray);
		rrhoArray[i]=returnArray[0];
		*j_chArray=returnArray[1];
		
		fsplevalIO(z_vht,c_vht,d_vht,x_vhart,r1Pointer[i],*mxGetPr(j_vht),returnArray);
		hpot00Array[i]=returnArray[0];
		*j_vhtArray=returnArray[1];
        
        iPlusOne=i+1;
        
        fsplevalIO(z_p_s,c_p_s,d_p_s,x_pot_s,r1Pointer[iPlusOne],*mxGetPr(j_p_s),returnArray);
		ppotArray[iPlusOne]=returnArray[0];
		*j_p_sArray=returnArray[1];
		
		fsplevalIO(z_chg,c_chg,d_chg,x_charg,r1Pointer[iPlusOne],*mxGetPr(j_ch),returnArray);
		rrhoArray[iPlusOne]=returnArray[0];
		*j_chArray=returnArray[1];
		
		fsplevalIO(z_vht,c_vht,d_vht,x_vhart,r1Pointer[iPlusOne],*mxGetPr(j_vht),returnArray);
		hpot00Array[iPlusOne]=returnArray[0];
		*j_vhtArray=returnArray[1];
	}

	plhs[0]=ppot;
	plhs[1]=rrho;
	plhs[2]=hpot00;


	return;
}

void fsplevalIO(mxArray *zArray,mxArray *cArray, mxArray *dArray, mxArray *xiArray, double x, double j_temp, double *returnArray){

    /*function [y j_out] = fsplevalIO(z,c,d,xi,x,j_in)
    %% [y j_out] = fsplevalB(z,c,d,xi,x,j_in)
    %% evaluates free spline function at x.
    %% x a scalar (only).
    %% z, c, d, as output from fspline..
    %% j_in  = input value for interval to try first
    %%         this interval is [xi(jin) xi(j_in+1)]
    %% j_out =  output interval for next call
    %%
    %% algorithm tests if x is in [x_[j_in] x_[j_in+1] ]
    %% if yes -- then computes the value -- if not
    %% then finds the correct interval.
    %%------------------------------------------
    */
    size_t n;
    double *xi,*c,*d,*z;
    long long ind_low,ind_high, ind_middle;
    long long j_in,j_out,j_outPlusOne;
    double t1,t2,h_j_out,y,val_middle;
    
    
    z=mxGetPr(zArray);
    c=mxGetPr(cArray);
    d=mxGetPr(dArray);
    xi=mxGetPr(xiArray);
    j_in=j_temp;

    n=mxGetM(xiArray);/* number of rows in xi*/
    
    /*------------------- input j_in is incorrect-- restart:*/
    if (j_in <0 || j_in > n-2){
         j_in = 0;
    }
    /*-------------------- x is outside interval -- bring it to
    //                     one of the boundary points*/
    if (x < xi[0]){
         x = xi[0];
         j_in = 0;
    }
    else if (x > xi[n-1]){
         x = xi[n-1];
         j_in = n-2;
    }
    j_out = j_in;
    /*-------------------- if point not in input interval
    //                     do Binary search
    // x is not between xi[j_in] and xi[j_in+1]*/
	if (!(xi[j_in]<=x && x<=xi[j_in+1])){
        ind_low = 0;
        ind_high = n-1;


        if (x>=xi[ind_low] && x<=xi[ind_high]){
            while (ind_high-ind_low>1){
                ind_middle = (ind_high+ind_low)/2;

                val_middle = xi[ind_middle];
                if (x<val_middle){
                    ind_high = ind_middle;
                }
                else {
                    ind_low = ind_middle;
                }
            }
            j_out = ind_low;
        }
        else {
            j_out=0;
            mexErrMsgTxt(" SPLINE ERROR [ in binary search ] ");
        } 
    }

    j_outPlusOne=j_out+1;
    
    t1  = xi[j_outPlusOne] - x;
    t2  = x - xi[j_out];
    h_j_out = xi[j_outPlusOne] - xi[j_out];
    y   = t1*(z[j_out]*t1*t1/h_j_out+c[j_out])+t2*(z[j_outPlusOne]*t2*t2/h_j_out+d[j_out]);

    
    returnArray[0]=y;
    returnArray[1]=j_outPlusOne;
    

    return;
}
