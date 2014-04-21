#ifndef SVM_FUNCS_HEADER_H
#define SVM_FUNCS_HEADER_H

#include "cv.h"
#include "ml.h"
#include <windows.h>
#include "Shlwapi.h"
#include "exceptions.h"
#include "consts.h"

using namespace std;
using namespace cv;

/*	struct Svm_Sample
*	represents one sample
*/
class DataSample
{ 
public:
	float m_label;
	float m_min_max_radius_ratio;
	float m_num_points_max;
	float m_num_points_min;
	float m_dom_flower_color_red;
	float m_dom_flower_color_green;
	float m_dom_flower_color_blue;
	float m_dom_inner_color_red;
	float m_dom_inner_color_green;
	float m_dom_inner_color_blue;
	float m_angle_min;
	float m_angle_max;
	float m_min_max_flower_ratio;
	float m_length_area_ratio;
	float m_length_inner_length_flower_ratio;
	static const int m_dimension = 14; //number of fields, not including "label"

	/*	svmSampleStructToSampleArray
	*	input: a struct representing one sample
	*	output: an array representing the same sample (needed as input for train)
	*/
	void toMatrix(Mat& sample_mat);

	void toCSV(string filepath);
};


/*	struct Svm_Data
*	collects all samples we got so far (together with their labels)
*/
//collects all the samples we  have
class TrainSet
{
public:
	Mat m_labels; //flower ids
	Mat m_samples; //samples corresponding to labels

	/*	svmAddSample
	*	input: a struct representing one sample (sample), and the struct representing all samples we have so far (sample_data)
	*	output: adds sample to samples_data
	*/
	void addSample(DataSample& sample);

	void toCSV(string filepath);
};


class SvmModel
{
private:
	SVM m_svm;

public:
	SvmModel() {};

	/*	svmPredict
	*	input: SVM class (which MUST be trained before passing it to this function), and a sample
	*	output: returns the label that SVM thinks the sample fits most
	*/
	double predict(DataSample& sample);

	/*	svmTrainFromSvmData
	*	input: SVM class, and a struct containing all samples data we have (samples_data)
	*	output: gives the struct to the SVM as train data
	*/
	void train(TrainSet& trainSet);

	void load(string filepath);
	void save(string filepath);
};

#endif