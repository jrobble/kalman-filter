#include "mex.h"
#include <matrix.h>
#include <string.h>
#include <math.h>

#include "math_util.h"
#include "calc_particle_probabilities.h";

#define PI 3.1415926535897932384626433832795

// get color histogram profile for a pixel
// function [cprofile] = get_color_profile(numbins,p,rcenters,gcenters,bcenters)
int* getColorProfile(int numbins, double* p, double* rcenters, double* gcenters, double* bcenters) {
    
    
    double* rdiffs = new double[numbins];
    subtract(p[0],rcenters,numbins,rdiffs);
    abs(rdiffs,numbins,rdiffs);
    
    double* gdiffs = new double[numbins];
    subtract(p[1],gcenters,numbins,gdiffs);
    abs(gdiffs,numbins,gdiffs);
    
    double* bdiffs = new double[numbins];
    subtract(p[2],bcenters,numbins,bdiffs);
    abs(bdiffs,numbins,bdiffs);
    
    int* cprofile = new int[3];
    
    double rret[2];
    min(rdiffs,numbins,rret);
    cprofile[0] = (int)(rret[1]); // get min index
    
    double gret[2];
    min(gdiffs,numbins,gret);
    cprofile[1] = (int)(gret[1]); // get min index
    
    double bret[2];
    min(bdiffs,numbins,bret);
    cprofile[2] = (int)(bret[1]); // get min index
    
    delete rdiffs;
    delete gdiffs;
    delete bdiffs;
    
    // DEBUG
    /*
    double* cprofile = new double[3];
    cprofile[0] = 1;
    cprofile[1] = 1;
    cprofile[2] = 1;
    */
    
    return cprofile;
}

// calc the color distribution
// function [py] = get_color_distribution(imgwidth,imgheight,numbins,xcoord,ycoord,rad,var,img,rcenters,gcenters,bcenters)
double* getColorDistribution(
        int imgwidth, int imgheight, int numbins, int xcoord, int ycoord, int rad, double var, 
        uint8* img, double* rcenters, double* gcenters, double* bcenters) {
    
    // determine center of ROI
    int centerx = xcoord;
    int centery = ycoord;
    
    // determine boundaries of ROI to check
    int ymin = (int)(max(1,centery-rad));
    int ymax = (int)(min(imgheight,centery+rad));
    int xmin = (int)(max(1,centerx-rad));
    int xmax = (int)(min(imgwidth,centerx+rad));
    
    // mexPrintf("\n>> >> >> ymin: %d",ymin); // DEBUG
    // mexPrintf("\n>> >> >> ymax: %d",ymax); // DEBUG
    // mexPrintf("\n>> >> >> xmin: %d",xmin); // DEBUG
    // mexPrintf("\n>> >> >> xmax: %d",xmax); // DEBUG
    
    
    // calculate color distribution at location y w/in the ROI 
    // for all m = 8x8x8 = 512 bins (if numbins == 8)
    double* py = new double[(int)pow(numbins,3)];
    double normfactor = 0;
    
    
    int u;
    double e, k;
    for(int y = ymin; y <= ymax; y++) {
        for(int x = xmin; x <= xmax; x++) {
            // determine if pixel within circle ROI
            if( pow((x-centerx),2) + pow((y-centery),2) <= pow(rad,2) ) {
                // absolute Euclidean distance divided by radius
                e = distance(x,centerx,y,centery);
                e = abs(e)/rad;
                
                k = 0;
                if(e < 1) { // should always be true unless random coordinate chosen
                    k = 1 - pow(e,2);  
                } else {
                    // e DEBUG
                    // error('e >= 1, %d >= 1',e);
                }

                // get color histogram profile
                double p[3] = {img[0*imgheight*imgwidth + x*imgheight + y],  // red
                               img[1*imgheight*imgwidth + x*imgheight + y],  // green
                               img[2*imgheight*imgwidth + x*imgheight + y]}; // blue
                  
                
                int* hbins = getColorProfile(numbins,p,rcenters,gcenters,bcenters);

                // determine bin index (like applying Kronecker delta)
                u = (hbins[2] * pow(numbins,2)) + (hbins[1] * numbins) + hbins[0];

                // calculate the norm factor once
                normfactor = normfactor + k;

                // calculate the probability density
                py[u] = py[u] + k;
                
                
                delete hbins;
            }
        }
    }
    
    // apply norm factor to distribution
    normfactor = 1/normfactor;
    mult(normfactor,py,(int)pow(numbins,3),py);
    
    
    return py;
}


// The gateway function
// [obj] = calc_particle_probabilities(obj,img)
void mexFunction( int nlhs, mxArray *plhs[], // outputs
                  int nrhs, const mxArray *prhs[]){ // inputs
    // mexPrintf(">> C++ calc_particle_probabilities\n"); // DEBUG
    
    // double *outMatrix; // output matrix
    const mxArray *objArray;
    mxArray *newobjArray;
    
    double *sk, *qccenters, *q;
    uint8 *img;
    int numpts, numbins;
    int imgheight, imgwidth, rad;
    double var;
    
    const char numptsStr[] = "numpts";
    const char numbinsStr[] = "numbins";
    const char skStr[] = "sk";
    const char qccentersStr[] = "qccenters";
    const char qStr[] = "q";
    const char imgheightStr[] = "imgheight";
    const char imgwidthStr[] = "imgwidth";
    const char radStr[] = "rad";
    const char varStr[] = "var";
    
    // pull out parameters
    objArray = prhs[0];
    // mexPrintf("\n>> >> >> objArray class: %s",mxGetClassName(objArray)); // DEBUG
    // mexPrintf("\n>> >> >> objArray field 0: %s",mxGetFieldNameByNumber(objArray,0)); // DEBUG
    
    numpts = (int)((mxGetPr(mxGetField(objArray,0,numptsStr)))[0]);
    // mexPrintf("\n>> >> >> numpts: %d",numpts); // DEBUG
    
    numbins = (int)((mxGetPr(mxGetField(objArray,0,numbinsStr)))[0]);
    // mexPrintf("\n>> >> >> numbins: %d",numbins); // DEBUG
    
    qccenters = mxGetPr(mxGetField(objArray,0,qccentersStr));
    /*
    mexPrintf("\n>> >> >> qccenters: [ "); // DEBUG
    for(int c = 0; c < 3; c++) {
        mexPrintf("\n[ ");
        for(int b = 0; b < numbins; b++) {
            mexPrintf("%f, ",qccenters[b*3 + c]);
        }
        mexPrintf("]");
    }
    mexPrintf("]");
    */
    
    q = mxGetPr(mxGetField(objArray,0,qStr));
    /*
    mexPrintf("\n>> >> >> q: [ "); // DEBUG
    for(int i = 0; i < pow(numbins,3); i++) {
        mexPrintf("%f, ",q[i]);
    }
    mexPrintf("]");
    */
    
    sk = mxGetPr(mxGetField(objArray,0,skStr));
    /*
    mexPrintf("\n>> >> >> sk: [ "); // DEBUG
    for(int s = 0; s < numpts; s++) {
        mexPrintf("\n[ ");
        for(int i = 0; i < 8; i++) {
            mexPrintf("%f, ",sk[i*numpts + s]);
        }
        mexPrintf("]");
    }
    mexPrintf("]");
    */
    
    imgheight = (int)((mxGetPr(mxGetField(objArray,0,imgheightStr)))[0]);
    // mexPrintf("\n>> >> >> imgheight: %d",imgheight); // DEBUG
    
    imgwidth = (int)((mxGetPr(mxGetField(objArray,0,imgwidthStr)))[0]);
    // mexPrintf("\n>> >> >> imgwidth: %d",imgwidth); // DEBUG
    
    rad = (int)((mxGetPr(mxGetField(objArray,0,radStr)))[0]);
    // mexPrintf("\n>> >> >> rad: %d",rad); // DEBUG
    
    var = (mxGetPr(mxGetField(objArray,0,varStr)))[0];
    // mexPrintf("\n>> >> >> var: %f",var); // DEBUG
    
    img = (uint8*)(mxGetPr(prhs[1]));
    /*
    mexPrintf("\n>> >> >> prhs[1] class: %s",mxGetClassName(prhs[1])); // DEBUG
    mexPrintf("\n>> >> >> img: [ "); // DEBUG
    for(int y = 0; y < 10; y++) {
        for(int x = 0; x < 10; x++) {
            mexPrintf("\n[ ");
            for(int c = 0; c < 3; c++) {
                mexPrintf("%d, ",img[c*imgheight*imgwidth + x*imgheight + y]);
            }
        }
        mexPrintf("]");
    }
    mexPrintf("]");
    */
    
    double* rcenters = new double[numbins];
    for(int i = 0; i < numbins; i++) {
        rcenters[i] = qccenters[i*3 + 0];
    }

    double* gcenters = new double[numbins];
    for(int i = 0; i < numbins; i++) {
        gcenters[i] = qccenters[i*3 + 1];
    }
    
    double* bcenters = new double[numbins];
    for(int i = 0; i < numbins; i++) {
        bcenters[i] = qccenters[i*3 + 2];
    }
    
    // compute
    double* py;
    double rho, constant, power, b, c = 0;
    int xcoord, ycoord;
    int size = ((int)(pow(numbins,3)));
    for(int i = 0; i < numpts; i++) { // DEBUG - numpts
        // fprintf(1,'Calculating particle color distribution, [%d,%d] ...\n',sk(i,1),sk(i,2)); // DEBUG
        
        xcoord = (int)(sk[0*numpts + i]);
        ycoord = (int)(sk[1*numpts + i]);
        py = getColorDistribution(imgwidth,imgheight,numbins,xcoord,ycoord,rad,var,
                img,rcenters,gcenters,bcenters);
       
        // calculate Bhattacharyya coefficient (rho)
        mult(py,q,size,py);
        sqrt(py,size,py);
        rho = sum(py,size);
        // fprintf(1,'rho: %d\n',ro);
        
        
        // DEBUG
        /*
        mexPrintf("  q: [ ");
        for(int j = 0; j < size; j++) {
            mexPrintf("%1.2f ",q[j]);
        }
        mexPrintf("]\n");
        mexPrintf(" py: [ ");
        for(int j = 0; j < size; j++) {
            mexPrintf("%1.2f ",py[j]);
        }
        mexPrintf("]\n");  
        mexPrintf("rho: %1.2f\n",rho);
        */
        
        //////////////////////////////////////////////////////////////////
        // START DEBUG
        //////////////////////////////////////////////////////////////////
        /*
        mxArray* input0 = mxCreateDoubleMatrix(1,1,mxREAL);
        (mxGetPr(input0))[0] = xcoord;
        mxArray* input1 = mxCreateDoubleMatrix(1,1,mxREAL);
        (mxGetPr(input1))[0] = ycoord;
        mxArray* input2 = mxCreateString("w+");
        mxArray* tmpplhs0[] = {input0,input1,input2};
        mxArray* output0 = mxCreateScalarDouble(0);
        mxArray* tmpprhs0[] = {output0};
        mexCallMATLAB(1,tmpprhs0,3,tmpplhs0,"plot");
        // int h = (int)((mxGetPr(tmpprhs0[0]))[0]);
        // mexPrintf("h: %d\n",h);
        
        input0 = mxCreateDoubleMatrix(1,1,mxREAL);
        (mxGetPr(input0))[0] = rad;
        input1 = mxCreateDoubleMatrix(1,1,mxREAL);
        (mxGetPr(input1))[0] = xcoord;
        input2 = mxCreateDoubleMatrix(1,1,mxREAL);
        (mxGetPr(input2))[0] = ycoord;
        mxArray* tmpplhs1[] = {input0,input1,input2};
        output0 = mxCreateScalarDouble(0);
        mxArray* tmpprhs1[] = {output0};
        mexCallMATLAB(1,tmpprhs1,3,tmpplhs1,"plot_circle");
        // int g = (int)((mxGetPr(tmpprhs1[0]))[0]);
        // mexPrintf("g: %d\n",g);
         
        mexCallMATLAB(0,NULL,0,NULL,"pause");
        
        // input0 = mxCreateScalarDouble(h);
        input0 = tmpprhs0[0];
        mxArray* tmpplhs2[] = {input0};
        mexCallMATLAB(0,NULL,1,tmpplhs2,"delete");
        
        // input0 = mxCreateScalarDouble(g);
        input0 = tmpprhs1[0];
        mxArray* tmpplhs3[] = {input0};
        mexCallMATLAB(0,NULL,1,tmpplhs3,"delete");
        */
        //////////////////////////////////////////////////////////////////
        // END DEBUG
        //////////////////////////////////////////////////////////////////

        
        // calculate probability
        constant = 1/(sqrt(2*PI*var));
        power = -(1-rho)/(2*var);
        b = constant*exp(power);
        sk[2*numpts + i] = b;
        
        // calculate cumulative probability
        c += b; 
        sk[3*numpts + i] = c;
        
         
        delete py;
    }
    
    // create the output matrix
    // plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // get a pointer to the real data in the output matrix
    // outMatrix = mxGetPr(plhs[0]);
    
    // copy over data
    // newobjArray = mxDuplicateArray(objArray);
    // mxSetData(mxGetField(newobjArray,0,skStr),sk);
    // plhs[0] = newobjArray;
    
    // mexPrintf("<< C++ calc_particle_probabilities\n"); // DEBUG
    
    delete rcenters;
    delete gcenters;
    delete bcenters;        
    
    return;
}
