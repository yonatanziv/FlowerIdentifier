#ifndef SVM_FUNCS_H
#define SVM_FUNCS_H

#include "cv.h"
#include "ml.h"
#include <windows.h>
#include "Shlwapi.h"
#include "exceptions.h"
#include "consts.h"

using namespace std;
using namespace cv;

/*	DataSample class
*	represents one sample - contains all the flower's features.
*/
class DataSample
{ 
public:
	float m_label;
	float m_inner_outer_radius_ratio;
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
	float m_length_area_ratio;					//The ratio between the length of the closest min points(to center) to the the radius of the flower.
	float m_length_inner_length_flower_ratio;
	static const int m_dimension = 14;			//number of fields, not including "label"

	/*	toMatrix
	*	output: sample_mat - matrix representation of the sample.
	*/
	void toMatrix(Mat& sample_mat);

	void toCSV(string filepath);
};


/*	TrainSet class
*	This class represents out training set. 
*/
class TrainSet
{
public:
	Mat m_labels;	//flower ids
	Mat m_samples;	//samples corresponding to labels

	/*	addSample
	*	input: sample
	*	adds sample to the training set.
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

	/*	predict
	*	input: A sample
	*	output: Returns the label that SVM thinks the sample fits most
	*/
	double predict(DataSample& sample);

	void train(TrainSet& trainSet);

	void load(string filepath);

	void save(string filepath);
};

#endif