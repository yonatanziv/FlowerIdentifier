#ifndef FLOWER_DETECTION_H
#define FLOWER_DETECTION_H

#include <highgui.h>
#include <math.h>
#include "svm_funcs.h"
#include "exceptions.h"
#include "utils.h"

#include <stdio.h>
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
#include "consts.h"

/*
 * FlowerFeatureExtractor class:
 * This class is responsible for extracting features out of image, given its center point.
 */
class FlowerFeatureExtractor
{
public:
	DataSample m_sample; //Will hold the flower's features.

	/*	extractFeatures
	 *	This is the main function which fills m_sample with the flower's features.
	 */
	void extractFeatures();

	/*	constructor
	 *	input: path of the image and center point of the flower.
	 */
	FlowerFeatureExtractor(string image_path, Point& center);

private:

	Mat m_image;							//The flower image.
	Mat m_draw_image;						//Image for drawing stuff (for debugging).
	Mat m_grabCutMask;
	Point m_center;							//Center point.
	vector<Point> m_max_flower_contour;		//Contour of the flower (Outer contour).
	


	/*	localMaxPoints
	*	input:	image - Binary image representing a contour (of the outer or inner flower).
	*			outImage - Same as image.
	*			resolution - The epsilon environment of the max points. (maximum points closer than <resolution> cells are not allowed)
	*	output: outImage - Will contain the value 2 at the cells where there are local maximum points.
    *			max_points - Vector of the local max points.
	*	
	*/
	void localMaxPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& max_points);

	/*	localMinPoints
	*	input:	image - Binary image representing a contour (of the outer or inner flower).
	*			outImage - Same as image.
	*			resolution - The epsilon environment of the min points. (minimum points closer than <resolution> cells are not allowed)
	*	output: outImage - Will contain the value 2 at the cells where there are local minimum points.
    *			min_points - Vector of the local min points.
	*	
	*/
	void localMinPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& min_points);

	/*
	 * localMinMaxPointsInner
	 * This function implements "localMaxPoints" and "localMinPoints" code.
	 * input: min_max_flag - flag which determines whether to find min or max points. 
	 */
	void localMinMaxPointsInner(Mat& image, Mat& outImage, int resolution, vector<Point>& min_max_points, int min_max_flag);


	/*	drawPointsOnImage
	*	input:	contour_matrix_helper - Output of "localMaxPoints" or "localMinPoints".
	*									This is binary contour matrix with value 2 at the extremum points.
	*			color - The wanted color of points.
	*	output: Draws the extremum points with specified color.
	*/
	void drawPointsOnImage(Mat& contour_matrix_helper, const Scalar& color);


	/* createCircleMask
	*	input: An image (for mask output), and a radius.
	*	output: A binary image containing a filled circle, with specified radius, around the center.
	*/
	void createCircleMask(Mat& img_mask, int radius);


	/* flipMask
	*	input: A binary mask
	*	output: The binary mask, with 1's and 0's switched
	*/
	void flipMask(Mat& img_mask);


	/* toBinaryByDominantColorWithContour
	*	input:	img_input - a color image (3 colors * 8 bits depth).
				dom_color - dominant color.
				threshold - threshold for binarization process.
				contour - (optional) - If its not empty, the binarization occures inside the contour area.
	*	output: img_output - Binary image, with white at all points inside contour (if contour is empty, ignores contour) which are "closer" to dom_color than the color threshold.
	*/
	void toBinaryByDominantColorWithContour(Mat& img_input, Mat& img_output, Scalar& dom_color, int threshold, vector<Point>& contour);


	/* findMaxContour
	*	input: contours - vector of contours.
	*	output: max_contour - will contain the contour from the vector with the largest amount of points.
	*/
	void findMaxContour(vector<vector<Point> >& contours, vector<Point>& max_contour);


	/* findClosestToContour
	*	input: contours - vector of contours.
	*	output: inner_contour - will contain the contour from the vector which is closest to the center point.
	*/
	void findClosestToCenterPointContour(vector<vector<Point> >& contours, vector<Point>& inner_contour);


	/* applyMeanShift
	 * Applying the meanshift algorithm on m_image. contains a fix to the size of the image, to allow the meanshift to work.
	 */
	void applyMeanShift();

	/*	grabCutImage
	*	* Applying the grabcut algorithm on m_image.
	*/
	void grabCutImage();

	/* createAndCalcHistogram
	*	input: an image and a mask (we use IplImage* in order to get 3-dimention histogram! (Impossible with matrix)).
	*	output: a pointer to an allocated histogram of the image area inside img_mask
	*/
	CvHistogram* createAndCalcHistogram(IplImage* img_hist, Mat& img_mask);


	/*	getOptimalInnerContour
	*	input:	
	*		dom_flower_color - dominant color of the flower.
	*		dom_inner_color - dominant color near the center point.
	*		circle_around_center_ponit_flag - this flag will be activated (by the calling function) if no inner contour found the first time this function is run. 
	*										  If this flag is on, the inner contour will be a circle around the center point with radius of the average distance 
	*										  between the center point and the contour with the smallest variance.
	*		outer_max_radius - radius of the flower. In order to eliminate bad contours - in case of the radius of the circle around center point is too big (more than half of the outer_max_radius).
	*	output: 
	*		runs on color thresholds (according to OPTIMAL_INNER_CONTOUR_START, OPTIMAL_INNER_CONTOUR_END, OPTIMAL_INNER_CONTOUR_QUOTIENT,
	*		paints the area which is within <threshold> distance from the center color, and within the outer flower contour, and tries to find contour there.
	*		picks the contour with smallest variance of distances between a point on the contour and the center point
	*	fail: max_inner_contour is empty if no inner contour found.
	*/
	void getOptimalInnerContour(vector<Point>& max_inner_contour, Scalar& dom_flower_color, Scalar& dom_inner_color, int circle_around_center_ponit_flag, double outer_max_radius);


	/*	getMedianAngle
	*	input:	min_max_points - vector of min or max points in size greater than 2.			
	*	output:	
	*		outMinAngle - The minimun angle between two points.
	*		outMinAngleTwoPointsIndexes - pointer to array of size two containing indexes of the twho points which create the min angle.
	*		return the median of all angles between two close points and the center.
	*/
	double getMedianAngle(vector<Point>& min_max_points, double& outMinAngle, int* outMinAngleTwoPointsIndexes);


	/*	getFixBadMinMaxPointsThreshold
	*	input: median of angles between min or max points
	*	output: threshold, according to medAngle (big, medium, or small threshold)
	*/
	int getFixBadMinMaxPointsThreshold(double medAngle);


	/*	fixBadMinMaxPoints
	*	This function removes bad min or max points - if the angle between two max points is smaller 
	*	than the medianAngle/fix_bad_min_max_points_threshold, it is removed.
	*	input: 
	*			min_max_points - vector of min or max points on the flower contour.	
	*			img_with_contours - The image with painted flower contour (for drawing reasons).
	*	output:	
	*			- ouMedAngle - will contain the new median angle.
	*			Returns the new number of min or max points after removing bad ones.
	*/
	int fixBadMinMaxPoints(vector<Point>& min_max_points, double& outMedAngle, Mat& img_with_contours);


	/*	createCircleMasks
	*	input:	
				mask_flag - changes inner mask radius, by steps same as COLOR_FLAG_STEP
	*			background_mask_flag - changes background mask radius, by steps same as COLOR_FLAG_STEP
	*	output: creates the required inner part, flower, and background binary masks according to flags
	*/
	void createCircleMasks(Mat& img_mask_small, Mat& img_mask_large, Mat& img_mask_background, int mask_flag, int background_mask_flag);


	/*	getDominantColors
	*	input:	
	*			mask_flag - changes inner mask radius 
	*			background_mask_flag - changes background mask radius
	*			img - The image we want to find dominant colors.
	*	output: dominant colors, according to the flags given
	*/
	void getDominantColors(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color, Mat& img, int mask_flag, int background_mask_flag);

	/*	getThresholdByFlag
	*	input:	color_flag - changes threshold which seperates background and flower color
	*			dom_flower_color - dominant color of the flower.
	*			dom_background_color - dominant color of the background.
	*	output: Returns the threshold which seperates background and flower color, according to the flag.
	*/
	int getThresholdByFlag(int color_flag, Scalar& dom_flower_color,Scalar& dom_background_color);


	/* checkBinaryBackground
	*	input: 
	*			img_binary - output of toBinaryByDominantColorWithContour.
	*			background_mask_flag - changes background mask radius.
	*	output: "masks" the image with background mask, and checks if there are more black pixels than white pixels in the image
	*	fail: throw ContourFailBinaryBackgroundException if there are more white pixels than black pixels in that area,
	*		  what would indicate a bad binarization - too much background.
	*/
	void checkBinaryBackground(Mat& img_binary, int background_mask_flag);


	/*	getOptimalFlowerContour
	*	input: 
	*			img - The image which we want to find the flower contour in (will be m_image or m_image after lightness fix).
	*			dom_flower_color - dominant color of the flower.
	*			dom_background_color - dominant color of the background.
	*			color_flag - changes color threshold.
	*			background_mask_flag - changes background mask radius.
	*	output: update m_max_flower_contour to be the optimal (largest) contour of outer flower (hopefully), according to the parameters.
	*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black in the area masked by background mask (fail of checkBinaryBackground).
	*		  CONTOUR_FAIL_CENTER_POINT if center point is not inside the contour.
	*/
	int getOptimalFlowerContour(Mat& img, Scalar& dom_flower_color, Scalar& dom_background_color, int color_flag, int background_mask_flag);


	/*	getOuterFlowerProperties
	*	output: updates sample parameters related to outer part of flower: 
	*	num of max and min points, median of angles between min points and max points, ratio between closest min point to the center and the flower radius.
	*/
	void getOuterFlowerProperties(int& num_points_max, int& num_points_min, double& angle_max, double& angle_min, double& min_max_flower_ratio, double& outer_max_radius);


	/*	getInnerPartProperties
	*	input: 
	*		dom_flower_color - dominant color of the flower.
	*		dom_inner_color - dominant color near center point.
	*		outer_max_radius - radius of the outer flower contour (in order to eliminate bad contours).
	*	output: updates sample properties related to the inner part of the flower: radius_inner and radius_length.
	*/
	void getInnerPartProperties(double& radius_inner, double& inner_length, Scalar& dom_flower_color, Scalar& dom_inner_color, double outer_max_radius);


	/*	getFlowerContourAndColor
	*	output: updates m_max_flower_contour to be the flower contour.
				updates the dominant colors of the flower, inner, and background.
	*	may try several background and inner mask radiuses, and several color thresholds, if it gets bad contours (no points, or center point not inside contour).
	*/
	void getFlowerContourAndColor(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color);

};


#endif