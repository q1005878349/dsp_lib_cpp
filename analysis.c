#include "analysis.h"

#include <malloc.h>
#include <math.h>

void data_integral_freq_one_mars(double* input, int length, double samp_freq, double lower_freq, double upper_freq, double *out)
{
	int num;
	double num2, num3;
	int i, num4, num5, min_num;;
	double* array3;
	complex* array2;
	complex* array4;
	complex* array5;
	complex zero, w_one, calc;
	zero.imaginary = 0;
	zero.real = 0;

	w_one.imaginary = 1.0;
	w_one.real = 0;

	num = (input != NULL) ? length : 0;

	array3 = (double*)calloc(num, sizeof(double));
	array2 = (complex*)calloc(num, sizeof(complex));
	array4 = (complex*)calloc(num, sizeof(complex));
	array5 = (complex*)calloc(num, sizeof(complex));

	for (i = 0; i < length; ++i) {
		array2[i].real = input[i];
		array2[i].imaginary = 0;
	}

	real_fft(array2, length);

	num2 = samp_freq / (float)num;
	num3 = PI * 2.0 * num2;
	for (i = 0; i < num; i++)
	{
		if (i <= num / 2 - 1)
			array3[i] = ((i == 0) ? 0.0 : (array3[i - 1] + num3));
		else
			array3[i] = ((i == num / 2) ? ((double)(-num) * num3 / 2.0) : (array3[i - 1] + num3));
	}

	array4[0].real = 0;
	array4[0].imaginary = 0;
	for (i = 1; i < num; i++)
	{
		calc.real = array3[i];
		calc.imaginary = 0;
		array4[i] = divide_complex(array2[i], multiply_complex(w_one, calc));
	}

	num4 = (int)floor(lower_freq / num2);
	num5 = (int)ceil(upper_freq / num2);
	min_num = (num - 1) > (num - num4) ? num - num4 : num - 1;

	for (int k = 0; k < num; k++)
		if (k >= num4 && k <= num5)
			array5[k] = array4[k];
		else if (k >= num - num5 && k <= min_num)
			array5[k] = array4[k];
		else
			array5[k] = zero;

	real_ifft(array5, num);

	for (i = 0; i < num; i++)
		out[i] = array5[i].real * 1000.0;

	free(array2);
	free(array3);
	free(array4);
	free(array5);
}

void data_integral_freq_one(double* input, int length, double samp_freq, double lower_freq, double upper_freq, double* out)
{
	int num;
	double num2, num3, last_w, cur_w;
	int i, num4, num5, min_num;;
	complex* array2;
	complex zero, w_one, calc;
	zero.imaginary = 0;
	zero.real = 0;

	w_one.imaginary = 1.0;
	w_one.real = 0;

	num = length;

	array2 = (complex*)calloc(num, sizeof(complex));

	for (i = 0; i < length; ++i)
		array2[i].real = input[i];

	real_fft(array2, length);

	num2 = samp_freq / (float)num;
	num3 = PI * 2.0 * num2;

	array2[0].real = 0;
	array2[0].imaginary = 0;
	last_w = 0.0;
	for (i = 1; i < num; i++)
	{
		if (i <= num / 2 - 1)
			cur_w = last_w + num3;
		else
			//cur_w = ((i == num / 2) ? ((double)(-num) * num3 / 2.0) : (last_w + num3));
			cur_w = ((i == num / 2) ? ((double)(-num) * num3 / 2.0 + num3) : (last_w + num3));

		last_w = cur_w;

		calc.real = cur_w;
		calc.imaginary = 0;
		array2[i] = divide_complex(array2[i], multiply_complex(w_one, calc));
	}

	num4 = (int)floor(lower_freq / num2);
	num5 = (int)ceil(upper_freq / num2);
	min_num = (num - 1) > (num - num4) ? num - num4 : num - 1;

	for (int k = 0; k < num; k++)
		if (!(k >= num4 && k <= num5) && !(k >= num - num5 && k <= min_num))
			array2[k] = zero;

	real_ifft(array2, num);

	for (i = 0; i < num; i++)
		out[i] = array2[i].real * 1000.0;

	free(array2);
}


void data_integral_freq_one_aes2(double* input, int length, double samp_freq, double lower_freq, double upper_freq, double* out)
{
	int num;
	double num2, num3, tmp;
	int i, num4, num5, min_num;;
	double* array3;
	complex* array2;
	complex* array4;
	complex* array5;
	complex zero, w_one, calc;
	zero.imaginary = 0;
	zero.real = 0;

	w_one.imaginary = 1.0;
	w_one.real = 0;

	num = (input != NULL) ? length : 0;

	array3 = (double*)calloc(num, sizeof(double));
	array2 = (complex*)calloc(num, sizeof(complex));
	array4 = (complex*)calloc(num, sizeof(complex));
	array5 = (complex*)calloc(num, sizeof(complex));

	for (i = 0; i < length; ++i) {
		array2[i].real = input[i];
		array2[i].imaginary = 0;
	}

	real_fft(array2, length);

	num2 = samp_freq / (float)num;
	num3 = PI * 2.0 * num2;
	for (i = 0; i < num; i++)
	{
		if (i <= num / 2 - 1)
			array3[i] = ((i == 0) ? 0.0 : (array3[i - 1] + num3));
		else
			array3[i] = ((i == num / 2) ? ((double)(num) * num3 / 2.0 - num3) : (array3[i - 1] - num3));
	}

	array4[0].real = 0;
	array4[0].imaginary = 0;
	for (i = 1; i < num; i++)
	{
		calc.real = array3[i];
		calc.imaginary = 0;
		array4[i] = divide_complex(array2[i], multiply_complex(w_one, calc));

		tmp = array4[i].real;
		array4[i].real = array4[i].imaginary;
		array4[i].imaginary = -tmp;
	}

	num4 = (int)floor(lower_freq / num2);
	num5 = (int)ceil(upper_freq / num2);
	min_num = (num - 1) > (num - num4) ? num - num4 : num - 1;

	for (int k = 0; k < num; k++)
		if (k >= num4 && k <= num5)
			array5[k] = array4[k];
		else if (k >= num - num5 && k <= min_num)
			array5[k] = array4[k];

	real_ifft(array5, num);

	for (i = 0; i < num; i++)
		out[i] = array5[i].real * 1000.0;

	free(array2);
	free(array3);
	free(array4);
	free(array5);
}

/* samplingRate 采样频率 input输入数据  length 输入数据长度  out 输出结果*/
void dtrend(double* input, int length,  double samplingRate, double* out)
{
	double* array = out;
	int num, i;
	double num2 = 1000.0 / samplingRate;
	double num3 = 0.0;
	double num4 = 0.0;
	double num5 = 0.0;
	double num6 = 0.0;
	double num7 = 0.0;
	double num8 = 0.0;

	num = length;

	for (i = 0; i < num; i++)
	{
		num3 = (double)((double)i * num2) + num3;
		num4 = (double)(((double)i * i) * num2 * num2) + num4;
		num5 += input[i];
		num6 = (double)((double)i * num2) * input[i] + num6;
	}
	num7 = (num3 * num5 - (double)num * num6) / (num3 * num3 - (double)num * num4);
	num8 = (num3 * num6 - num4 * num5) / (num3 * num3 - (double)num * num4);
	for (i = 0; i < num; i++)
	{
		array[i] = input[i] - num7 * (double)i * (double)num2 - num8;
	}

	//array[num - 1] = array[num - 2];
}

void cmpression_byte_array(double* input, int length, double coef, unsigned char *out)
{
	double num;
	unsigned short tvalue;

	for (int i = 0; i < length; i++)
	{
		num = input[i] / (double)coef;
		if (num < -32768.0)
		{
			num = -32768.0;
		}
		if (num > 32767.0)
		{
			num = 32767.0;
		}

		tvalue = (unsigned short)num;
		out[i * 2] = (unsigned char)tvalue;
		out[i * 2 + 1] = (unsigned char)(tvalue >> 8);
	}
}

void acc_to_vel(double* input, int length, double samp_freq, double lower_freq, double coef, double *out)
{
	data_integral_freq_one(input, length, samp_freq, lower_freq, samp_freq / 2, out);
	dtrend(out, length, samp_freq, out);
}


double rms(double* input, int length)
{
	double square_total = 0.0, ave;
	int i;

	ave = mean(input, length);

	for (i = 0; i < length; ++i)
	{
		square_total += (input[i] - ave) * (input[i] - ave);
	}
	return (double)sqrt(square_total / (double)length);
}

double mean(double* input, int length)
{
	int i;
	double total = 0;
	
	for (i = 0; i < length; ++i) {
		total += input[i];
	}

	return total / length;
}

void remove_mean(complex* inout, int length)
{
	int i;
	double total = 0;
	double mean = 0;

	for (i = 0; i < length; i++) {
		total += inout[i].real;
	}

	mean = total / length;

	for (i = 0; i < length; i++) {
		inout[i].real -= mean;
	}
}

/*input 去掉均值后*/
double kurtosis(double* input, int length)
{
	int i;
	double total = 0;

	if (length == 0)
		return 0;

	for (i = 0; i < length; ++i) {
		total += pow(fabs(input[i]), 4.0);
	}

	return total / (length * pow(rms(input, length), 4.0));
}

/*长度一定要2的n次方  output 长度如何确定*/
void fft_calc(complex* input, int length, double* output, int outlen)
{
	int i, num, num2, num3;
	double total, ave;
	
	num = length;
	num2 = ceil(log(num) / log(2.0));
	num3 = pow(2.0, num2);
	
	total = 0;
	for (i = 0; i < length; ++i)
		total += input[i].real;
	ave = total / length;

	for (i = 0; i < num3; ++i)
		input[i].real -= ave;

	real_fft(input, num3);
	if (num3 != 0)
	{
		for (i = 0; i < outlen; i++)
		{
			output[i] = sqrt(input[i].real * input[i].real + input[i].imaginary * input[i].imaginary) / (double)(num / 2);
		}
	}
}

int butterworth_bandpass_filter(double* input_data, int length, double sampling_rate, double low_cutoff_frequency, double high_cutoff_frequency, double *output_data)
{
	complex* array;
	int num, i;

	num = length;

	if ((low_cutoff_frequency > high_cutoff_frequency) || (num == 0))
		return -1;
	
	array = (complex*)calloc(num, sizeof(complex));

	for (i = 0; i < length; ++i) {
		array[i].real = input_data[i];
		array[i].imaginary = 0;
	}

	real_fft(array, length);

	for (int i = 0; i < num / 2; i++)
	{
		if (sampling_rate * (double)i / (double)num < low_cutoff_frequency)
		{
			array[i].real = 0;
			array[i].imaginary = 0;

			array[num - 1 - i].real = 0;
			array[num - 1 - i].imaginary = 0;
		}
		else if (sampling_rate * (double)i / (double)num > high_cutoff_frequency) {
			array[i].real = 0;
			array[i].imaginary = 0;

			array[num - 1 - i].real = 0;
			array[num - 1 - i].imaginary = 0;
		}
	}

	real_ifft(array, length);

	for (i = 0; i < length; ++i)
		output_data[i] = array[i].real;

	return 0;
}


int dm_calc(double* data, int length,  float sampling_rate, double max_freq, double center_freq, double filter_width, float width_rate, double** output, int* out_length)
{
	double num2, num3, *array3;
	int num, num4, i;
	complex* array4;

	num = length;
	array3 = (double *)calloc(num, sizeof(double));
	array4 = (complex*)calloc(num, sizeof(complex));

	num2 = center_freq + filter_width / 2.0;
	num3 = center_freq - filter_width / 2.0;
	if (num3 == 0.0)
		num3 = 0.001;

	if (num2 > max_freq || num3 < 0.0)
		return 1;
	else if (fabs(num2 - num3) / 2.0 <= (double)(sampling_rate / (float)num))
		return 2;
	else
	{
		butterworth_bandpass_filter(data, length, sampling_rate, num3, num2, array3);

		for (i = 0; i < num; ++i) {
			array4[i].real = array3[i];
			array4[i].imaginary = 0;
		}

		hilbert(array4, num);

		for (i = 0; i < num; i++)
		{
			array4[i].real = sqrt(array4[i].imaginary * array4[i].imaginary + array3[i] * array3[i]);
			array4[i].imaginary = 0;
		}

		num4 = (filter_width / 2.0 * (double)width_rate / (double)(sampling_rate / (float)num));
		if (num4 > num)
			num4 = num;

		if (num4 <= 0) {
			num4 = 0;
		}
		else {
			*output = (double*)malloc(sizeof(double) * num4);

			fft_calc(array4, num, *output, num4);
		}
		*out_length = num4;
	}

	free(array3);
	free(array4);

	return 0;
}

