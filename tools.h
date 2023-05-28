#ifndef __TOOLS_H__
#define __TOOLS_H__
#include <math.h>

#ifdef __cplusplus

extern "C" {

#endif

//#define PI 3.14159265358979
#define PI acos(-1)


typedef struct {
	double real;
	double imaginary;
} complex;

complex multiply_complex(const complex a, const complex b);
complex add_complex(const complex a, const complex b);
complex subtract_complex(const complex a, const complex b);
complex divide_complex_number(const complex a, const float b);
complex divide_complex(const complex a, const complex b);

void real_fft(complex* inout, int length);
void real_ifft(complex* inout, int length);
void hilbert(complex* input, int lenghth);
void rceps(double* input, int length, double* output);
void hanning(double* inout, int length);
void hanning_complex(complex* inout, int length);


#ifdef __cplusplus

}

#endif 
#endif