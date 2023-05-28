#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "tools.h"

#ifdef __cplusplus
extern "C" {
#endif

	void data_integral_freq_one(double*, int, double, double, double, double*);
	void data_integral_freq_one_mars(double*, int, double, double, double, double*);
	void data_integral_freq_one_aes2(double*, int, double, double, double, double*);
	void dtrend(double*, int, double, double*);
	void ccmpression_byte_array(double*, int, double, unsigned char*);
	int butterworth_bandpass_filter(double* input_data, int length, double sampling_rate, double low_cutoff_frequency, double high_cutoff_frequency, double* output_data);

	void acc_to_vel(double*, int, double, double, double, double*);
	double rms(double*, int);
	double mean(double*, int);
	void remove_mean(complex*, int);
	double kurtosis(double*, int);
	int dm_calc(double*, int, float, double, double, double, float, double**, int*);
#ifdef __cplusplus
}
#endif 

#endif