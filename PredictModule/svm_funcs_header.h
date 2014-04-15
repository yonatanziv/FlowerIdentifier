#ifndef SVM_FUNCS_HEADER_H
#define SVM_FUNCS_HEADER_H

#include "cv.h"
#include "ml.h"
#define PI 3.1415926f
using namespace std;
using namespace cv;

/*	struct Svm_Data
*	collects all samples we got so far (together with their labels)
*/
struct Svm_Data{ //collects all the samples we have 
	float* labels; //flower ids
	float* samples; //samples corresponding to labels
	int max_sample_amount;
	int current_sample_amount;
	int sample_dimension; //how many values in each sample
};

/*	struct Svm_Sample
*	represents one sample
*/
struct Svm_Sample{ 
	static const int dimension = 14; //number of fields, not including "label"
	float label;
	float min_max_radius_ratio;
	float num_points_max;
	float num_points_min;
	float dom_flower_color_red;
	float dom_flower_color_green;
	float dom_flower_color_blue;
	float dom_inner_color_red;
	float dom_inner_color_green;
	float dom_inner_color_blue;
	float angle_min;
	float angle_max;
	float min_max_flower_ratio;
	float length_area_ratio;
	float length_inner_length_flower_ratio;
};


/*	svmInitData
*	input: a struct which will contain data about samples we give the svm, maximum amount of samples, and the dimension (number of fields) of every sample
*	output: initializes the struct: allocates needed memory
*/
void svmInitData(Svm_Data* samples_data, int max_sample_amount, int sample_dimension);

/*	svmDataFree
*	input: a struct containing all samples data
*	output: releases the struct
*/
void svmDataFree(Svm_Data* samples_data);

/*	svmSampleStructToSampleArray
*	input: a struct representing one sample
*	output: an array representing the same sample (needed as input for train)
*/
void svmSampleStructToSampleArray(Svm_Sample* sample, float* sample_array);

/*	svmAddSample
*	input: a struct representing one sample (sample), and the struct representing all samples we have so far (sample_data)
*	output: adds sample to samples_data
*/
void svmAddSample(Svm_Data* samples_data, Svm_Sample* sample);

/*	svmTrainFromSvmData
*	input: SVM class, and a struct containing all samples data we have (samples_data)
*	output: gives the struct to the SVM as train data
*/
void svmTrainFromSvmData(CvSVM* svm, Svm_Data* samples_data);

/*	svmPredict
*	input: SVM class (which MUST be trained before passing it to this function), and a sample
*	output: returns the label that SVM thinks the sample fits most
*/
double svmPredict(CvSVM* svm, Svm_Sample* sample);

#endif