#include "svm_funcs_header.h"

/*	svmSampleStructToSampleArray
*	input: a struct representing one sample
*	output: an array representing the same sample (needed as input for train)
*/
void DataSample::toMatrix(Mat& sample_mat)
{ //"manually" transfers sample struct to an array
	//NOTE: no need to put sample->label in array
	sample_mat = Mat(1, m_dimension, CV_32FC1);

	sample_mat.at<float>(0,0) = m_min_max_radius_ratio;
	sample_mat.at<float>(0,1) = m_dom_flower_color_red/255;
	sample_mat.at<float>(0,2) = m_dom_flower_color_green/255;
	sample_mat.at<float>(0,3) = m_dom_flower_color_blue/255;
	sample_mat.at<float>(0,4) = m_dom_inner_color_red/255;
	sample_mat.at<float>(0,5) = m_dom_inner_color_green/255;
	sample_mat.at<float>(0,6) = m_dom_inner_color_blue/255;
	sample_mat.at<float>(0,7) = m_angle_min/(PI);
	sample_mat.at<float>(0,8) = m_angle_max/(PI);
	sample_mat.at<float>(0,9) = m_num_points_min/30;
	sample_mat.at<float>(0,10) = m_num_points_max/30;
	sample_mat.at<float>(0,11) = m_min_max_flower_ratio;
	sample_mat.at<float>(0,12) = m_length_area_ratio;
	sample_mat.at<float>(0,13) = m_length_inner_length_flower_ratio;
}

/*	svmAddSample
*	input: a struct representing one sample (sample), and the struct representing all samples we have so far (sample_data)
*	output: adds sample to samples_data
*/

void TrainSet::addSample(DataSample& sample)
{ //adds a sample to the Svm_Data struct
	Mat sample_without_label;
	sample.toMatrix(sample_without_label);

	m_samples.push_back(sample_without_label);
	m_labels.push_back<float>(sample.m_label);
}

void SvmModel::load(string filepath) 
{
	if(!PathFileExists(filepath.c_str())){
		throw PathNotFoundException();
	}

	m_svm.load(filepath.c_str());
}

void SvmModel::save(string filepath) 
{
	m_svm.save(filepath.c_str());
}

void SvmModel::train(TrainSet& trainSet) 
{
	m_svm.train(trainSet.m_samples, trainSet.m_labels);
}

double SvmModel::predict(DataSample& sample) 
{
	Mat sample_mat;
	sample.toMatrix(sample_mat);
	return m_svm.predict(sample_mat);
}

