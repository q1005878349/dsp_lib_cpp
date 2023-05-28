
/* Replace "dll.h" with the name of your header */
#include "dsp.h"
#include <iostream>
#include <malloc.h>

#include "analysis.h"

using namespace std;

JNIEXPORT void JNICALL Java_com_gobigg_cloud_common_util_iot_DspLib_fft(JNIEnv* env,
													jobject object,
													jdoubleArray inout)
{
	if (inout == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: inout NULL in fft");
		return;
	}
	int len = env->GetArrayLength(inout);

	if (!len || len % 1024 != 0) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: length error in fft");
		return;
	}

	complex* comp = (complex *)(env->GetDoubleArrayElements(inout, NULL));

	//¼ÓººÄþ´°
	hanning_complex(comp, len / 2);

	real_fft(comp, len / 2);

	env->SetDoubleArrayRegion(inout, 0, len, (const jdouble *)comp);
}

JNIEXPORT void JNICALL Java_com_gobigg_cloud_common_util_iot_DspLib_hanning_complex(JNIEnv* env, jobject object, jdoubleArray inout)
{
	if (inout == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: inout NULL in hanning complex");
		return;
	}
	int len = env->GetArrayLength(inout);

	if (!len || len % 1024 != 0) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: length error in hanning complex");
		return;
	}

	complex* comp = (complex*)(env->GetDoubleArrayElements(inout, NULL));

	//¼ÓººÄþ´°
	hanning_complex(comp, len / 2);

	env->SetDoubleArrayRegion(inout, 0, len, (const jdouble*)comp);
}

JNIEXPORT void JNICALL Java_com_gobigg_cloud_common_util_iot_DspLib_accToVel(JNIEnv* env, jobject object,
		jdoubleArray input, jdouble samp_freq, jdouble lower_freq, jdouble coef, jdoubleArray output)
{
	if (input == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: input NULL in accToVel");
		return;
	}
	int len = env->GetArrayLength(input);

	if (!len || len % 1024 != 0) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: length error in accToVel");
		return;
	}


	jdouble* data = env->GetDoubleArrayElements(input, NULL);
	jdouble* out = env->GetDoubleArrayElements(output, NULL);

	acc_to_vel(data, len, samp_freq, lower_freq, coef, out);

	env->SetDoubleArrayRegion(output, 0, len, out);
}


JNIEXPORT void JNICALL Java_com_gobigg_cloud_common_util_iot_DspLib_rceps(JNIEnv* env, jobject object, jdoubleArray inout)
{
	if (inout == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: inout NULL in rceps");
		return;
	}

	jdouble* data = env->GetDoubleArrayElements(inout, NULL);
	int length = env->GetArrayLength(inout);

	if (!length || length % 1024 != 0) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: length error in rceps");
		return;
	}

	rceps(data, length, data);

	env->SetDoubleArrayRegion(inout, 0, length, data);
}


JNIEXPORT jdoubleArray JNICALL Java_com_gobigg_cloud_common_util_iot_DspLib_dmcalc(JNIEnv* env, jobject object,
	jdoubleArray data, jfloat sampling_rate, jdouble max_freq, jdouble center_freq, jdouble filter_width, jfloat width_rate)
{
	jdouble *output;
	int out_length = 0;
	jdouble* input = NULL;

	if (data == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: data NULL in dmcalc");
		return NULL;
	}

	input = env->GetDoubleArrayElements(data, NULL);
	int in_length = env->GetArrayLength(data);

	dm_calc(input, in_length, sampling_rate, max_freq, center_freq, filter_width, width_rate, &output, &out_length);

	jdoubleArray outArray = env->NewDoubleArray((jsize)out_length);

	if (outArray == NULL) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: cannot allocate out array in dmcalc");
		return NULL;
	}

	if (out_length <= 0) {
		env->ThrowNew(env->FindClass("java/lang/Exception"), "exception from jni: out length is 0 in dmcalc");
		return NULL;
	}
	
	env->SetDoubleArrayRegion(outArray, 0, (jsize)out_length, output);

	free(output);

	return outArray;
}

