#include <malloc.h>

#include "tools.h"

complex multiply_complex(const complex a, const complex b)
{
	complex temp;
	temp.real = a.real * b.real - a.imaginary * b.imaginary;
	temp.imaginary = a.real * b.imaginary + a.imaginary * b.real;
	return temp;
}
complex add_complex(const complex a, const complex b)
{
	complex temp;
	temp.real = a.real + b.real;
	temp.imaginary = a.imaginary + b.imaginary;
	return temp;
}
complex subtract_complex(const complex a, const complex b)
{
	complex temp;
	temp.real = a.real - b.real;
	temp.imaginary = a.imaginary - b.imaginary;
	return temp;
}

complex divide_complex_number(const complex a, const float b)
{
	complex temp;
	temp.real = a.real / b;
	temp.imaginary = a.imaginary / b;
	return temp;
}

complex divide_complex(const complex a, const complex b)
{
	complex temp;
	temp.real = (a.real * b.real + a.imaginary * b.imaginary) / (b.real * b.real + b.imaginary * b.imaginary);
	temp.imaginary = (a.imaginary * b.real - a.real * b.imaginary) / (b.real * b.real + b.imaginary * b.imaginary);
	return temp;
}
static void reverse(complex* inout, int length)
{
	complex   temp;
	int   i = 0, j = 0, k = 0;
	double   t;
	for (i = 0; i < length; i++)
	{
		k = i; j = 0;
		t = log2(length);
		while ((t--) > 0)
		{
			j = j << 1;
			j |= (k & 1);
			k = k >> 1;
		}
		if (j > i)
		{
			temp = inout[i];
			inout[i] = inout[j];
			inout[j] = temp;
		}
	}
}

void real_fft(complex* inout, int length)
{
	int i = 0, j = 0, k = 0, l = 0, index;
	complex up, down, product, w;
	
	reverse(inout, length);
	//fft运算
	for (i = 0; i < log2(length); i++)
	{
		l = 1 << i;
		for (j = 0; j < length; j += 2 * l)
		{
			for (k = 0; k < l; k++)
			{
				index = length * k / 2.0 / l;
				w.real = cos(2 * PI / length * (index));
				w.imaginary = -1 * sin(2 * PI / length * (index));
				//product = multipy_complex(x[j + k + l] , W[length * k / 2 / l]);
				product = multiply_complex(inout[j + k + l], w);
				up = add_complex(inout[j + k] , product);
				down = subtract_complex(inout[j + k] , product);
				inout[j + k] = up;
				inout[j + k + l] = down;
			}
		}
	}
}


void real_ifft(complex* inout, int length)
{
	int   i = 0, j = 0, k = 0, l = length, index;
	complex   up, down, w;
	complex *x = inout;

	for (i = 0; i < log(length) / log(2); i++) 
	{
		l /= 2;
		for (j = 0; j < length; j += 2 * l)
		{
			for (k = 0; k < l; k++)
			{
				index = length * k / 2.0 / l;
				w.real = cos(2 * PI / length * (index));
				w.imaginary = -1 * sin(2 * PI / length * (index));

				up = add_complex(x[j + k] , x[j + k + l]);
				up = divide_complex_number(up, 2);
				down = subtract_complex(x[j + k], x[j + k + l]);
				down = divide_complex(divide_complex_number(down, 2), w);
				x[j + k] = up;
				x[j + k + l] = down;
			}
		}
	}
	reverse(x, length);
}

void hilbert(complex* inout, int length)
{
	complex* xn = inout;

	//fft
	real_fft(xn, length);

	for (int n = 1; n < length; n++)//得到Z(K)
	{
		if (n < length / 2)
		{
			xn[n].real = xn[n].real * 2.0;
			xn[n].imaginary = xn[n].imaginary * 2.0;
		} else {
			xn[n].real = 0;
			xn[n].imaginary = 0;
		}
	}

	real_ifft(xn, length);
}

void rceps(double* input, int length, double* output)
{
	complex* output_fft = (complex*)calloc(length, sizeof(complex));

	for (int i = 0; i < length; ++i) {
		output_fft[i].real = input[i];
		output_fft[i].imaginary = 0.0;
	}

	real_fft(output_fft, length);

	for (int i = 0; i < length; ++i) {
		output[i] = sqrt(output_fft[i].real * output_fft[i].real + output_fft[i].imaginary * output_fft[i].imaginary);
	}

	for (int i = 0; i < length; ++i) {
		if (output[i] != 0)
			output[i] = log(output[i]);
		output_fft[i].real = output[i];
		output_fft[i].imaginary = 0.0;
	}

	real_ifft(output_fft, length);
	
	for (int i = 0; i < length; ++i) {
		output[i] = output_fft[i].real;
	}

	free(output_fft);
}

void hanning_internal(double* window, int length) {
	int i;
	double pi, phase = 0, delta;

	pi = 4. * atan(1.0);
	delta = 2 * pi / ((double)length - 1);

	for (i = 0; i < length / 2; i++) {
		window[i] = (double)(0.5 * (1.0 - cos(phase)));
		phase += delta;
	}

	for (i = length / 2; i < length; i++) {
		window[i] = window[length - i - 1];
	}
}

void hanning(double* inout, int length) {
	int i;
	double* window = (double*)calloc(length, sizeof(double));

	hanning_internal(window, length);

	for (i = 0; i < length; ++i)
		inout[i] *= window[i];

	free(window);
}
static float* hanning_m(int N, short itype)
{
	int half, i, idx, n;
	float* w;

	w = (float*)calloc(N, sizeof(float));
	memset(w, 0, N * sizeof(float));

	if (itype == 1)	//periodic function
		n = N - 1;
	else
		n = N;

	if (n % 2 == 0)
	{
		half = n / 2;
		for (i = 0; i < half; i++) //CALC_HANNING   Calculates Hanning window samples.
			w[i] = 0.5 * (1 - cos(2 * PI * (i + 1) / (n + 1)));

		idx = half - 1;
		for (i = half; i < n; i++) {
			w[i] = w[idx];
			idx--;
		}
	}
	else
	{
		half = (n + 1) / 2;
		for (i = 0; i < half; i++) //CALC_HANNING   Calculates Hanning window samples.
			w[i] = 0.5 * (1 - cos(2 * PI * (i + 1) / (n + 1)));

		idx = half - 2;
		for (i = half; i < n; i++) {
			w[i] = w[idx];
			idx--;
		}
	}

	if (itype == 1)	//periodic function
	{
		for (i = N - 1; i >= 1; i--)
			w[i] = w[i - 1];
		w[0] = 0.0;
	}
	return(w);
}

void hanning_complex(complex* inout, int length)
{
	int i;

	float* window = hanning_m(length, 1);

	for (i = 0; i < length; ++i) {
		inout[i].real *= (window[i] * 2);//2 是修正系数 理论上应该计算完汉宁后，在外部乘。此处为节省计算量一并乘上
	}

	free(window);

}

void hamming_internal(double* window, int length) {
	int i;
	double phase = 0, delta;

	delta = 2 * PI / ((double)length - 1);

	for (i = 0; i < length / 2; i++) {
		window[i] = (double)(0.54 - 0.46 * cos(phase)); 
		phase += delta;
	}

	for (i = length / 2; i < length; i++) {
		window[i] = window[length - i - 1];
	}
	
}
