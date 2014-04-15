#include "svm_funcs_header.h"

/*	svmInitData
*	input: a struct which will contain data about samples we give the svm, maximum amount of samples, and the dimension (number of fields) of every sample
*	output: initializes the struct: allocates needed memory
*/
void svmInitData(Svm_Data* samples_data, int max_sample_amount, int sample_dimension){ //initializes Svm_Data struct
	samples_data->labels = (float*)malloc(max_sample_amount*sizeof(float));
	samples_data->samples = (float*)malloc(max_sample_amount * sample_dimension * sizeof(float));
	samples_data->max_sample_amount = max_sample_amount;
	samples_data->sample_dimension = sample_dimension;
	samples_data->current_sample_amount = 0;
}

/*	svmDataFree
*	input: a struct containing all samples data
*	output: releases the struct
*/
void svmDataFree(Svm_Data* samples_data){ //frees memory from pointers in samples_data struct
	free(samples_data->labels);
	samples_data->labels=NULL;
	free(samples_data->samples);
	samples_data->samples=NULL;
}

/*	svmSampleStructToSampleArray
*	input: a struct representing one sample
*	output: an array representing the same sample (needed as input for train)
*/
void svmSampleStructToSampleArray(Svm_Sample* sample, float* sample_array){ //"manually" transfers sample struct to an array
	//NOTE: no need to put sample->label in array
	
	sample_array[0]=sample->min_max_radius_ratio;
	sample_array[1]=sample->dom_flower_color_red/255;
	sample_array[2]=sample->dom_flower_color_green/255;
	sample_array[3]=sample->dom_flower_color_blue/255;
	sample_array[4]=sample->dom_inner_color_red/255;
	sample_array[5]=sample->dom_inner_color_green/255;
	sample_array[6]=sample->dom_inner_color_blue/255;
	sample_array[7]=sample->angle_min/(PI);
	sample_array[8]=sample->angle_max/(PI);
	sample_array[9]=sample->num_points_min/30;
	sample_array[10]=sample->num_points_max/30;
	sample_array[11]=sample->min_max_flower_ratio;
	sample_array[12]=sample->length_area_ratio;
	sample_array[13]=sample->length_inner_length_flower_ratio;

	for(int i=0; i<sample->dimension; i++){
		//printf("%.2f ", sample_array[i]);
	}
	printf("\n");
}

/*	svmAddSample
*	input: a struct representing one sample (sample), and the struct representing all samples we have so far (sample_data)
*	output: adds sample to samples_data
*/
void svmAddSample(Svm_Data* samples_data, Svm_Sample* sample){ //adds a sample to the Svm_Data struct
	float* sample_array=(float*)malloc(samples_data->sample_dimension * sizeof(float));
	svmSampleStructToSampleArray(sample,sample_array);

	if(samples_data->current_sample_amount < samples_data->max_sample_amount){
		samples_data->labels[samples_data->current_sample_amount]=sample->label;

		memcpy(samples_data->samples + samples_data->current_sample_amount * samples_data->sample_dimension, 
				sample_array, samples_data->sample_dimension * sizeof(float)); //copies array to our array

		samples_data->current_sample_amount++;
	}
	free(sample_array);
}

/*	svmTrainFromSvmData
*	input: SVM class, and a struct containing all samples data we have (samples_data)
*	output: gives the struct to the SVM as train data
*/
void svmTrainFromSvmData(CvSVM* svm, Svm_Data* samples_data){ //gets samples (with labels) and trains svm according to them
	CvMat labels,trainData;
	cvInitMatHeader(&labels,samples_data->current_sample_amount,1,CV_32FC1,samples_data->labels);
	cvInitMatHeader(&trainData,samples_data->current_sample_amount,samples_data->sample_dimension,CV_32FC1,samples_data->samples);
	//printf("trainData.cols: %d trainData.rows: %d ",trainData.cols,trainData.rows);
	svm->train(&trainData,&labels);
}

/*	svmPredict
*	input: SVM class (which MUST be trained before passing it to this function), and a sample
*	output: returns the label that SVM thinks the sample fits most
*/
double svmPredict(CvSVM* svm, Svm_Sample* sample){ //gets a sample, returns which label it thinks the sample should fit
	double res;
	float* predict_arr=(float*)malloc(sample->dimension * sizeof(float));
	CvMat predict_mat;
	svmSampleStructToSampleArray(sample,predict_arr);
	cvInitMatHeader(&predict_mat,sample->dimension,1,CV_32FC1,predict_arr);
	res = svm->predict(&predict_mat);
	free(predict_arr); //don't free before finished using the Mat!
	return res;
}