#ifndef CALC_PARTICLE_PROBABILITIES_H
#define CALC_PARTICLE_PROBABILITIES_H

typedef unsigned char uint8;

// get color histogram profile for a pixel
int* getColorProfile(int numbins, double* p, 
                        double* rcenters, double* gcenters, double* bcenters);

// calc the color distribution
double* getColorDistribution(
        int imgwidth, int imgheight, int numbins, int xcoord, int ycoord, int rad, double var, 
        uint8* img, double* rcenters, double* gcenters, double* bcenters);

// The gateway function
void mexFunction( int nlhs, mxArray *plhs[], // outputs
                  int nrhs, const mxArray *prhs[]); // inputs

#endif /* CALC_PARTICLE_PROBABILITIES_H */