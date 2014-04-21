#ifndef UTILS_H
#define UTILS_H

#include "flower_detection_header.h"
#include "dirent.h"


class Utils
{
public:
	static DIR* MyChdir(char * path);

	static void toLowercase(char* str);


	/*	getCenter
	*	input: A center point pointer(to be update), a center points file, a number of a flower picture
	*	output: Updates center to contain the center of the flower, as written in the file. return 0 on success and -1 on fail.
	*	file format: each line contains:
	*				flowerPicNo,x,y
	*/
	static int getCenter(Point& center, FILE* centers,int flowerPicNo);

	static void resizeImage(Mat& image, Mat& out);

	/*show image if DEBUG MODE is on*/
	static void Utils::myShowImage(const char* name, const Mat& image);


	static void Utils::showAndWait(const char* name, Mat& image);


	static void Utils::convertRGBtoHSL (Scalar& rgb_color, Scalar& hsl_color);


	static double Utils::hue2rgb(double p, double q, double t); //helper for convertHSLtoRGB

	static void Utils::convertHSLtoRGB(Scalar& hsl_color, Scalar& rgb_color);

	static void Utils::bgr2rgb(Scalar& bgr, Scalar& rgb);

	/*return x^2*/
	static int square(int x);

	static bool isGrayscale(Scalar& color);

	/*	equalizeLightness
	*	input: an IplImage pointer
	*	output: sets lightness in all pixels to be 50 (shadow remover).
	*/
	static void equalizeLightness(Mat& img);

		/*	distance
	*	input: two points (x1,y1),(x2,y2)
	*	output: euclidean distance between them
	*/
	static double distance(int y1, int x1, int y2, int x2);


	/*	distance
	*	input: two CvPoint points
	*	output: euclidean distance between them
	*/
	static double point_distance(Point& point1, Point& point2);


	/*	padded_array_access
	*	input: an int array, its size, and a cell in it (row and col)
	*	output: returns image[row][col] if the indices are inside the array, 0 otherwise
	*	used to prevent falling out of array in localMaxPoints and localMinPoints
	*/
	static int padded_matrix_access(Mat& image, int x, int y);

		/*	colorDst
	*	input: two colors
	*	output: returns euclidean distance between the points (as three-dimensional points (R,G,B))
	*/
	static int colorDst(Scalar& color1, Scalar& color2);

	static void drawOneContour(Mat& draw_image, vector<Point>& contour,CvScalar& color);

	/*	angle
	*	input: gets three points
	*	output: calculates angle between the three points
	*/
	static double angle(Point& center, Point& p1, Point& p2);

	static double innerProduct(Point& v1, Point& v2);

	/*	toDegree
	*	input: angle in radians
	*	output: angle in degrees
	*/
	static double toDegree(double medAngle);

	/* contourVar
	*	input: a contour and flower's center
	*	output: returns variance of distances between a point on the contour and the center point
	*/
	static double contourVar(vector<Point>& contour, Point& center);

	/*avgDst
	*input: contour and center point.
	*output: The average distance between the center point and the contour. 
	*/
	static double avgDst(vector<Point>& contour, Point& center);

	/*
	*	input:	- one point of the array
	- array of points
	- isVisited - the index of points to be ignored. (chosen to be closest to other points already).
	- currentPointIndex - the index of the above point. its also need to be ignored.
	output:	- the index of the closest point
	*/
	static int getClosestPoint(int currentPointIndex,vector<Point>& min_max_points, vector<int>& isVisited);

	/*	median
	*	input: an array and its size
	*	output: median of the values in the array
	*/
	static double median(vector<double>& vec);

	/* radiusOutOfPoints
	input: an array of points, its size, a center point
	output: median of the distances between a point from the array and the center point
	*/
	static double radiusOutOfPoints(vector<Point>& max_points, Point& center);


	/*	getMinMaxFlowerRatio
	*	input: minimum points on outer flower contour, radius of outer flower contour
	*	output: returns ratio between closest minimum point to the center and the outer contour ratio
	*/
	static double getMinMaxFlowerRatio(vector<Point> min_points, double radius_max, Point& center);

	/* contourToArray
	*	input: a contour, a two-dim allocated array
	*	output: the contour points in a two-dim int array (which represents a binary image). contour points will have 1 in array
	*/
	static void contourToMatrix(vector<Point> contour, Mat& arr);

	static double contourToPointDst(vector<Point>& contour, Point& center);


};




#endif