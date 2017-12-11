#include "mex.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    
    /* input variables */
    mxArray *vxc;
    double ndim;
    double twovpia0;
    /* output variables*/
    
    
    /* local variables */
    double *vxcArray;
    double *newvxcArray;
    double p75vpi = 0.75/3.14159265358979323846;
    double third=1.0/3.0;
    double rs,sqrs,ec,alpha,rho;
    double exc=0.0;
    
    double g  =-0.2846, b1 = 1.0529;
    double b2 = 0.3334, c1 = 0.0622;
    double c2 = 0.096,  c3 = 0.004;
    double c4 = 0.0232, c5 = 0.0192;
    
    size_t i;
    
    if (nrhs!=3){
        mexErrMsgTxt("Incorrect number of input arguments"); 
    }
    if (nlhs!=2){
        mexErrMsgTxt("Incorrect number of output arguments"); 
    }
    
    vxc=prhs[0];
    ndim=*mxGetPr(prhs[1]);
    twovpia0=*mxGetPr(prhs[2]);
    
    vxcArray=mxGetPr(vxc);
    newvxcArray=mxMalloc(mxGetM(vxc)*mxGetN(vxc)*sizeof(double));
    
    for (i=0;i<ndim;i++){
         rho = vxcArray[i];
         newvxcArray[i] = 0;
         if (rho > 0){ 
            rs = pow((p75vpi/rho),third);/* to the 1/3 power */
            newvxcArray[i] = -twovpia0/rs;
            exc = exc + 0.75*rho*newvxcArray[i];
            if (rs >= 1){
               sqrs = sqrt(rs);
               ec = g/(1 + b1*sqrs + b2*rs);
               newvxcArray[i] +=ec*ec*(1+3.5*b1*sqrs*third+4*b2*rs*third)/g;
            }
            else{
               alpha = log(rs);
               ec = c1*alpha - c2 + (c3*alpha - c4)*rs;
               newvxcArray[i] += + ec - (c1 + (c3*alpha - c5)*rs)*third;
            }
            exc += rho*ec;
         }
    }

    plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(plhs[0],newvxcArray);
    mxSetM(plhs[0],mxGetM(vxc));
    mxSetN(plhs[0],mxGetN(vxc));    
    plhs[1]=mxCreateDoubleScalar(exc);
    
    return;
}    
    
    