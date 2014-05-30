#ifndef CONSTS_H
#define CONSTS_H

#include <iostream>
using namespace std;

class Consts
{
public:
	static const float PI;
	static const int CONTOUR_FAIL_NO_POINTS;							//failure return flag - no points on contour found
	static const int CONTOUR_FAIL_CENTER_POINT;							//failure return flag - center point not in contour
	static const int CONTOUR_FAIL_BINARY_BACKGROUND;					//failure return flag - binary background image contains more white than black
	static const int SUCCESS;											//Success flag.
	static const int NUM_BINS_FOR_COLOR;								//number of histogram bins per color axis (actually there are NUM_BINS_FOR_COLOR^3 bins)
	static const int MAX_POINTS_THRESHOLD;								//padding threshold for max local points algorithm (minimum cells between two max points)
	static const int MIN_POINTS_THRESHOLD;								//padding threshold for min local points algorithm (minimum cells between two min points)
	static const int MEANSHIFT_LEVEL;									//parameter for meanshift algorithm
	static const int MEANSHIFT_SPATIAL_RADIUS;							//parameter for meanshift algorithm
	static const int MEANSHIFT_COLOR_RADIUS;							//parameter for meanshift algorithm
	static const float CONTRAST_RATIO;									//level of contrast we apply to image, as preprocessing
	static const float CONTRAST_RATIO_OUTER_CONTOUR;					//level of contrast we apply to image, when finding outer contour with lightness fix
	static const float LIGHTNESS_THRESHOLD_WHITE;						//if lightness is higher than this threshold, color is considered white
	static const float LIGHTNESS_THRESHOLD_BLACK;						//if lightness is lower than this threshold, color is considered black
	static const float SATURATION_THRESHOLD_GRAY;						//if saturation is lower than this threshold, color is considered gray
	static const float LIGHTNESS_THRESHOLD_FIX;							//lightness fix raises lightness only if original lightness is higher than this threshold
	static const float LIGHTNESS_VALUE_AFTER_FIX;						//lightness fix raises lightness to this value
	static const int BASE_SIZE;											//resize image to this size
	static const int RADIUS_SMALL;										//radius of small binary mask 
	static const int OUTER_COLOR_DST_THRESHOLD_DENOMINATOR;				//denominator of flower contour detection color thresholds. (determines step size)
	static const int OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR;			//start looking for inner contour, at START/DENOMINATOR of the color distance between flower dominant color and background dominant color
	static const int OUTER_COLOR_DST_THRESHOLD_MAX_NUMERATOR;			//finish looking for inner contour, at MAX/DENOMINATOR of the color distance between flower dominant color and background dominant color
	static const int COLOR_FLAG_ITERATIONS;								//changes color threshold in outer flower contour detection (iterative algorithm which tries different thresholds)
	static const float COLOR_FLAG_STEP;									//changes the step of color threshold in outer flower contour detection (iterative algorithm which tries different thresholds)
	static const int BACKGROUND_MASK_FLAG_ITERATIONS;					//changes background mask radius in outer flower contour detection (iterative algorithm which tries different radiuses)
	static const int OPTIMAL_INNER_CONTOUR_START;						//start looking for inner contour, at START/QUOTIENT of the color distance between inner dominant color and flower dominant color
	static const int OPTIMAL_INNER_CONTOUR_END;							//end looking for inner contour, at END/QUOTIENT of the color distance between inner dominant color and flower dominant color
	static const int OPTIMAL_INNER_CONTOUR_QUOTIENT;					//quotient of color threshold ratio, as described above. (determines step size)
	static const int FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG_ANGLES;		//angles larger than this one are considered "big"
	static const int FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL_ANGLES;		//angles smaller than this one are considered "small". angles larger than SMALL_ANGLES and smaller than BIG_ANGLES are considered "medium"
	static const int FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG;				//flag - big threshold for fix bad min max points algorithm
	static const int FIX_BAD_MIN_MAX_POINTS_THRESHOLD_MEDIUM;			//flag - medium threshold for fix bad min max points algorithm
	static const int FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL;			//flag - small threshold for fix bad min max points algorithm
	static const int INNER_CONTOUR_EQUAL_FLOWER_CONTOUR_THRESHOLD;		//This threshold tells that if the the color distance between the inner contour and the flower contour is smaller than its value than the inner contour is same as the flower contour.
	static const int MINIMUM_INNER_CONTOUR_POINTS;						//This number is the minimum number of points that inner contour should contain.
	static const float INNER_RADIUS_TO_FLOWER_RADIUS_RATIO;				//In order to eliminate bad circles contours around center point.
	static const int INNER_CONTUOR_MINIMUM_THRESHOLD;					//Minimum threshold to find inner contour. (Minimum color distance threshold used at the binarization process)
	static const int MAX_POINTS_FLAG;									//This flag is passed to localMinMaxPointsInner function if we want to find local max points.
	static const int MIN_POINTS_FLAG;									//This flag is passed to localMinMaxPointsInner function if we want to find local min points.
	static const int GRUBCUT_ITERATIONS;								//Number of iterations of grabcut algorithm.
	static const int DEBUG;												//Debug mode indicator.
	static const string SVM_SERIALIZE_PATH;
	static const string SAMPLES_DIR_PATH;
	static const string SVM_SERIALIZATION_DIR_PATH;
	static const string CENTERS_FILE_NAME;
	static const string SVM_SERIALIZE_FILE_NAME;
	static const string PREDICT_DIR_PATH;
	static const string ALL_CENTERS_FILE_NAME;
	
};

#endif