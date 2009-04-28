#ifndef MATH_UTIL_H
#define MATH_UTIL_H
         
// limits library not fully implemented in Open Watcom, so use Visual C++ 2008 compiler
#include <limits>
#include <iostream>

// Visual C++ 2008 compiler has no pow(int,int) implementation
int pow(int a, int p) {
    return (int)(pow((double)a,p));
}

// element-by-element matrix square root
void sqrt(double* a, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = sqrt(a[i]);
    }
}

// matrix maximum
void max(double* a, int size, double* ret) {
    double max = -std::numeric_limits<double>::infinity();
    double maxindex = -1;
    for(int i = 0; i < size; i++) {
        if(a[i] > max) {
            max = a[i];
            maxindex = i;
        }
    }
    ret[0] = max;
    ret[1] = maxindex;
}
         
// matrix minimum
void min(double* a, int size, double* ret) {
    // double min = mxGetInf(); // corrupts obj data
    double min = std::numeric_limits<double>::infinity();
    double minindex = -1;
    for(int i = 0; i < size; i++) {
        if(a[i] < min) {
            min = a[i];
            minindex = i;
        }
    }
    ret[0] = min;
    ret[1] = minindex;
}

// element-by-element absolute matrix values
void abs(double* a, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = abs(a[i]);
    }
}
         
// element-by-element matrix multiplication
void mult(double a1, double* a2, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = a1 * a2[i];
    }
}

// element-by-element matrix multiplication
void mult(double* a1, double* a2, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = a1[i] * a2[i];
    }
}
         
// element-by-element matrix substraction
void subtract(double a1, double* a2, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = a1 - a2[i];
    }
}
     
// element-by-element matrix substraction
void subtract(double* a1, double* a2, int size, double* ret) {
    for(int i = 0; i < size; i++) {
        ret[i] = a1[i] - a2[i];
    }
}
         
// calc Euclidean distance
double distance(double x1, double x2, double y1, double y2) {
    return sqrt( pow(x2-x1,2) + pow(y2-y1,2) );
}
         
// calc maximum
double max(double first, double second) {
    if(first >= second) {
        return first;
    } else {
        return second;
    }
}
         
// calc minimum
double min(double first, double second) {
    if(first <= second) {
        return first;
    } else {
        return second;
    }
}

// calc product
double prod(double array[], int size) {
    double prod = 0;
    if(size > 0) {
        prod = 1;
        for(int i = 0; i < size; i++) {
            prod *= array[i];
        }
    }
    
    // DEBUG
    /*
    mexPrintf("prod size: %d\n",size); // DEBUG
    for(int i = 0; i < size; i++) {
        mexPrintf("array[%d]: %f\n",i,array[i]); // DEBUG
    }
    
    mexPrintf("prod: %f\n",prod); // DEBUG
    */
    return prod;
}

// calculate mean squared error
double calc_fitnessMSE(double absErr) {
    double fitnessMSE = 100000 * (1 / (1 + pow(absErr,2)));
    return fitnessMSE;
}

// round number
int round(double num) {
    return (int) floor(num+.5);
}

// calculate sum
double sum(double* array, int size) {
    double sum = 0;
    for(int i = 0; i < size; i++) {
        sum += array[i];
    }
    return sum;
}

// calculate mean
// note that there is no way to get the size of a dynamic array, so pass it in
double mean(double* array, int size) {
    return sum(array,size) / size;
}

// calculate median
double median(double* array, int size) {
    // quickSort(array, size); // actually, already sorted
    /*
    mexPrintf("pre-sorted:\n"); // DEBUG
    for(int i = 0; i < size; i++) {
        mexPrintf("%f, ", array[i]); // DEBUG
    }
    mexPrintf("\n"); // DEBUG
    
    mexPrintf("post-sorted:\n"); // DEBUG
    for(int i = 0; i < size; i++) {
        mexPrintf("%f, ", array[i]); // DEBUG
    }
    mexPrintf("\n"); // DEBUG
    */
    
    double median = -1;
    if (size % 2 == 1) { // odd
        median = array[(int)size/2];
    } else { // even
        median = (array[(int)size/2] + array[((int)(size/2))-1]) / 2;
    }
    return median;
}

#endif /* MATH_UTIL_H */
