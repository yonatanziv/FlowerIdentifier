#ifndef FLOWER_DETECTION_HEADER_H
#define FLOWER_DETECTION_HEADER_H
//version 58

#include <highgui.h>
//#include <cv.h>
#include <math.h>
#include "svm_funcs_header.h"
#include "exceptions.h"
#include "utils.h"

#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <direct.h>
#include  <io.h>
#include <iostream>
#include <exception>

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <math.h>
#include <time.h>
#include <limits.h>

#define PI 3.1415926f
#define OUTER_COLOR_DST_THRESHOLD_DENOMINATOR 16 //denominator of flower contour detection color thresholds. (determines step size)
#define OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR 8 //start looking for inner contour, at START/DENOMINATOR of the color distance between flower dominant color and background dominant color
#define OUTER_COLOR_DST_THRESHOLD_MAX_NUMERATOR 12 //finish looking for inner contour, at MAX/DENOMINATOR of the color distance between flower dominant color and background dominant color
#define COLOR_FLAG_ITERATIONS 4 //changes color threshold in outer flower contour detection (iterative algorithm which tries different thresholds)
#define COLOR_FLAG_STEP 0.05 //changes the step of color threshold in outer flower contour detection (iterative algorithm which tries different thresholds)
#define BACKGROUND_MASK_FLAG_ITERATIONS 4 //changes background mask radius in outer flower contour detection (iterative algorithm which tries different radiuses)
#define CONTOUR_FAIL_NO_POINTS -2 //failure return flag - no points on contour found
#define CONTOUR_FAIL_CENTER_POINT -3 //failure return flag - center point not in contour
#define CONTOUR_FAIL_BINARY_BACKGROUND -4 //failure return flag - binary background image contains more white than black
#define SUCCESS 0
#define NUM_BINS_FOR_COLOR 30 //number of histogram bins per color axis (actually there are NUM_BINS_FOR_COLOR^3 bins)
#define MAX_POINTS_THRESHOLD 7 //padding threshold for max local points algorithm (minimum cells between two max points)
#define MIN_POINTS_THRESHOLD 8 //padding threshold for min local points algorithm (minimum cells between two min points)
#define MEANSHIFT_LEVEL 1 //parameter for meanshift algorithm
#define MEANSHIFT_SPATIAL_RADIUS 20 //parameter for meanshift algorithm
#define MEANSHIFT_COLOR_RADIUS 40 //parameter for meanshift algorithm
#define CONTRAST_RATIO 1.2 //level of contrast we apply to image, as preprocessing
#define CONTRAST_RATIO_OUTER_CONTOUR 1.5 //level of contrast we apply to image, when finding outer contour with lightness fix
#define LIGHTNESS_THRESHOLD_WHITE 0.8 //if lightness is higher than this threshold, color is considered white
#define LIGHTNESS_THRESHOLD_BLACK 0.2 //if lightness is lower than this threshold, color is considered black
#define SATURATION_THRESHOLD_GRAY 0.15 //if saturation is lower than this threshold, color is considered gray
#define LIGHTNESS_THRESHOLD_FIX 0.2 //lightness fix raises lightness only if original lightness is higher than this threshold
#define LIGHTNESS_VALUE_AFTER_FIX 0.5 //lightness fix raises lightness to this value
#define BASE_SIZE 200 //resize image to this size
#define RADIUS_SMALL 3 //radius of small binary mask 
#define OPTIMAL_INNER_CONTOUR_START 10 //start looking for inner contour, at START/QUOTIENT of the color distance between inner dominant color and flower dominant color
#define OPTIMAL_INNER_CONTOUR_END 20 //end looking for inner contour, at END/QUOTIENT of the color distance between inner dominant color and flower dominant color
#define OPTIMAL_INNER_CONTOUR_QUOTIENT 30 //quotient of color threshold ratio, as described above. (determines step size)
#define FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG_ANGLES 50 //angles larger than this one are considered "big"
#define FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL_ANGLES 25 //angles smaller than this one are considered "small". angles larger than SMALL_ANGLES and smaller than BIG_ANGLES are considered "medium"
#define FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG 2 //flag - big threshold for fix bad min max points algorithm
#define FIX_BAD_MIN_MAX_POINTS_THRESHOLD_MEDIUM 3 //flag - medium threshold for fix bad min max points algorithm
#define FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL 4 //flag - small threshold for fix bad min max points algorithm
#define INNER_CONTOUR_EQUAL_FLOWER_CONTOUR_THRESHOLD 5 //This threshold tells that if the the color distance between the inner contour and the flower contour is smaller than its value than the inner contour is same as the flower contour.
#define MINIMUM_INNER_CONTOUR_POINTS 5 //This number is the minimum number of points that inner contour should contain.
#define INNER_RADIUS_TO_FLOWER_RADIUS_RATIO 0.5 //In order to eliminate bad circles contours around center point.
#define INNER_CONTUOR_MINIMUM_THRESHOLD 10 //Minimum threshold to find inner contour. 
#define SVM_SERIALIZE_PATH "samples\\serialization\\svmSerialize.txt"
#define MAX_POINTS_FLAG 0
#define MIN_POINTS_FLAG 1
#define DEBUG 1

class FlowerFeatureExtractor
{
public:
	DataSample m_sample; //Will hold the features.

	/*	getSampleDataFromImage
	*	input: image path, samples struct, center point
	*	output: updates samples struct with sample parameters of image in image_path. uses center point for calculations
	*	return: -1 on fail (image dropped because of problems with sample struct, 0 on success
	*/
	void extractFeaturesFromImage();

	FlowerFeatureExtractor(string image_path, Point& center);

private:

	Mat m_image;
	Mat m_draw_image;
	Point m_center;
	vector<Point> m_max_flower_contour;


	void localMinMaxPointsInner(Mat& image, Mat& outImage, int resolution, vector<Point>& min_max_points, int min_max_flag);

	/*	localMaxPoints
	*	input: an int array representing a binary image, its size, a point, an output array, resolution
	*	output: outImage will have 2 at the cells where there are local maximum points of the distance function between the given point and the array point
	*	not allows minimum points closer than <resolution> cells
	*/
	void localMaxPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& max_points);

	/*	localMinPoints
	*	input: an int array representing a binary image, its size, a point, an output array, resolution
	*	output: outImage will have 2 at the cells where there are local minimum points of the distance function between the given point and the array point
	*	not allows minimum points closer than <resolution> cells
	*/
	void localMinPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& min_points);


	/*	drawArrPointsOnImage
	*	input: an image, an array to draw, colors and such
	*	output: draws the extremum points from the array on the image, with specified color, thickness, etc
	*/
	void drawPointsOnImage(Mat& contour_matrix_helper, const Scalar& color);



	/* createCircleMask
	*	input: an image (for mask output), radius, center
	*	output: a binary image containing a filled circle, with specified radius, around center
	*/
	void createCircleMask(Mat& img_mask, int radius);


	/* flipMask
	*	input: a binary mask
	*	output: the binary mask, with 1's and 0's switched
	*/
	void flipMask(Mat& img_mask);


	/* toBinaryByDominantColorWithContour
	*	input: a color image (3 colors * 8 bits depth), a pointer to output image (1 color, 8 bits depth), color threshold, contour pointer
	*	output: img_output is a binary image, with white at all points inside contour (if contour is null, ignores contour) which are "closer" to dom_color than the color threshold
	*/
	void toBinaryByDominantColorWithContour(Mat& img_input, Mat& img_output, Scalar& dom_color, int threshold, vector<Point>& contour);


	/* findMaxContour
	*	input: a pointer to a list of contours (first_contour), a pointer to a pointer of a contour
	*	output: max_contour will contain the contour from the list with the largest amount of points on the contour
	*/
	void findMaxContour(vector<vector<Point> >& contours, vector<Point>& max_contour);


	/* findClosestToContour
	*	input: a pointer to a list of contours (first_contour), a pointer to a pointer of a contour
	*	output: contour will contain the contour from the list which is closest to the center point.
	*/
	void findClosestToCenterPointContour(vector<vector<Point> >& contours, vector<Point>& inner_contour);


	/* applyMeanShift
	*	input: an image
	*	output: the image, after applying the meanshift algorithm on it. contains a fix to the size of the image, to allow the meanshift to work
	*/
	void applyMeanShift();


	/* createAndCalcHistogram
	*	input: an image and a mask
	*	output: a pointer to an allocated histogram of the image area inside img_mask
	*/
	CvHistogram* createAndCalcHistogram(IplImage* img_hist, IplImage* img_mask);


	/* createDistanceArray
	*	input: an array of points, its size, a center point, a pointer to an allocated distances array
	*	output: updates distances array to contain all distances between a point from the points array and the center point
	*/
	void createDistanceArray(double* distances,CvPoint* points,int num_points,CvPoint* center);


	/*	getOptimalInnerContour
	*	input: max_flower_conotur (contour of outer flower), pointer to max_inner_contour pointer (to update pointer), dominant colors, 
	circle_around_center_ponit_flag - this flag will be activated (by the calling function) if no inner contour found the first time this function is run. if this flag is on, the inner contour will be a circle around the center point with radius of the average distance between the center point and the contour with the smallest variance.
	radius_max in order to elimnate bad contours - in case of the radius of the circle around center point is too big (more than half of the radius_max).
	*	output: runs on color thresholds (according to OPTIMAL_INNER_CONTOUR_START, OPTIMAL_INNER_CONTOUR_END, OPTIMAL_INNER_CONTOUR_QUOTIENT,
	*			paints the are which is within <threshold> distance from the center color, and within the outer flower contour, and tries to find contour there.
	picks the contour with smallest variance of distances between a point on the contour and the center point
	*	fail: max_inner_contour=NULL if no inner contour found.
	*/
	void getOptimalInnerContour(vector<Point>& max_inner_contour, Scalar& dom_flower_color, Scalar& dom_inner_color, int circle_around_center_ponit_flag, double radius_max);


	/*
	*	input:	- array of point (min or max points) in size greater than 2.
	*			- center point
	*	output:	The median of all angles between two close points and the center.
	*/
	double getMedianAngle(vector<Point>& min_max_points, double& outMinAngle, int* outMinAngleTwoPointsIndexes);


	/*	getFixBadMinMaxPointsThreshold
	*	input: median of angles between min or max points
	*	output: threshold, according to medAngle (big, medium, or small threshold)
	*/
	int getFixBadMinMaxPointsThreshold(double medAngle);


	/*
	*	This function removes bad max points - if the angle between two max points is smaller 
	*	than the medianAngle/fix_bad_min_max_points_threshold, it is removed.
	*	
	*	output:	The function update pointsArr and return its new length.
	*			- ouMedAngle - will contain the new median angle. 
	*/
	int fixBadMinMaxPoints(vector<Point>& min_max_points, double& outMedAngle, Mat& img_with_contours);


	/*	createCircleMasks
	*	input: binary mask (IplImage) pointers, mask_flag (changes inner mask radius, by steps same as COLOR_FLAG_STEP), background_mask_flag (changes background mask radius, by steps same as COLOR_FLAG_STEP)
	*	output: creates the required inner part, flower, and background binary masks according to flags
	*/
	void createCircleMasks(Mat& img_mask_small, Mat& img_mask_large, Mat& img_mask_background, int mask_flag, int background_mask_flag);


	/*	getDominantColors
	*	input: dominant color struct pointers, mask_flag (changes inner mask radius), background_mask_flag (changes background mask radius)
	*	output: dominant colors, according to the flags given
	*/
	void getDominantColors(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color, Mat& img, int mask_flag, int background_mask_flag);

	/*	getThresholdByFlag
	*	input: color_flag (changes threshold which seperates background and flower color )
	*	output: the threshold which seperates background and flower color, according to the flag
	*/
	int getThresholdByFlag(int color_flag,Scalar& dom_flower_color,Scalar& dom_background_color);


	/* checkBinaryBackground
	*	input: binary image, background_mask_flag (changes background mask radius)
	*	output: "masks" the image with background mask, and checks if there are more black pixels than white pixels in the image
	*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black pixels in that area (what would indicate a bad binarization - too much background)
	*/
	void checkBinaryBackground(Mat& img_binary, int background_mask_flag);


	/*	getOptimalFlowerContour
	*	input: pointer to pointer of outer flower contour (in order to update pointer), dominant colors, color flag (changes color threshold), background mask flag (changes background mask radius)
	*	output: gets optimal (largest) contour of outer flower (hopefully), according to the parameters.
	*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black in the area masked by background mask (fail of checkBinaryBackground)
	: CONTOUR_FAIL_CENTER_POINT if center point is not inside the contour
	*/
	int getOptimalFlowerContour(Mat& img, Scalar& dom_flower_color, Scalar& dom_background_color, int color_flag, int background_mask_flag);


	/*	getOuterFlowerProperties
	*	input: outer contour of flower, and parameters to update
	*	update: updates sample parameters related to outer part of flower: num of max and min points, median of angles between min points and max points, ratio between closest min point to the center and the flower radius
	*/
	void getOuterFlowerProperties(int& num_points_max, int& num_points_min, double& angle_max, double& angle_min, double& min_max_flower_ratio, double& radius_max);


	/*	getInnerPartProperties
	*	input: outer contour of the flower(in order to create a mask, to define where we look for the inner contour), radius_inner to update, colors and so on (radius_max in order to elimnate bad contours)
	*	output: updates sample properties related to the inner part of the flower. for now: radius_inner
	*/
	void getInnerPartProperties(double& radius_inner, double& inner_length, Scalar& dom_flower_color, Scalar& dom_inner_color, double radius_max);


	/*	getFlowerContourAndColor
	*	input: pointer to pointer of max_flower_contour (in order to update it), color struct pointers
	*	output: updates outer flower contour and dominant colors of the flower, inner, and background.
	*	may try several background and inner mask radiuses, and several color thresholds, if it gets bad contours (no points, or center point not inside contour)
	*/
	void getFlowerContourAndColor(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color);

};


#endif