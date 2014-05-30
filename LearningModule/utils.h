#ifndef UTILS_H
#define UTILS_H

#include "flower_detection.h"
#include "dirent.h"


class Utils
{
public:
	static DIR* MyChdir(string path);

	static void toLowercase(char* str);


	/*	getCenter
	*	input: A center point (to be update), a center points file, a number of a flower picture.
	*	output: Updates center to contain the center of the flower, as written in the file. 
	*			returns 0 on success and -1 on fail.
	*	file format: each line contains:
	*				flowerPicNo,x,y
	*/
	static int getCenter(Point& center, FILE* centers,int flowerPicNo);

	/*	resizeImage
	 *	input: an image
	 *	output: out - image after resize.
	 */
	static void resizeImage(Mat& image, Mat& out);

	/*show image if DEBUG MODE is 1*/
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
	*	input: an image
	*	output: sets lightness in all pixels to be 50 (shadow remover).
	*/
	static void equalizeLightness(Mat& img);

	/*	distance
	*	input: two points (x1,y1),(x2,y2)
	*	output: euclidean distance between them
	*/
	static double distance(int y1, int x1, int y2, int x2);

	/*	pointDistance
	*	input: two Point points
	*	output: euclidean distance between them
	*/
	static double pointDistance(Point& point1, Point& point2);


	/*	paddedMatrixAccess
	*	input: an image, and a cell in it (row and col)
	*	output: returns image[row][col] if the indices are inside the array, 0 otherwise
	*	used to prevent falling out of array in localMinMaxPointsInner.
	*/
	static int paddedMatrixAccess(Mat& image, int x, int y);

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

	/*	contourVar
	*	input: a contour and flower's center
	*	output: returns variance of distances between a point on the contour and the center point
	*/
	static double contourVar(vector<Point>& contour, Point& center);

	/*	avgDst
	 *	input: contour and center point.
	 *	output: The average distance between the center point and the contour. 
	 */
	static double avgDst(vector<Point>& contour, Point& center);

	/*
	 *	input:	
	 *			currentPointIndex - the index of the above point. its also need to be ignored.
	 *			min_max_points - vector of points
	 *			isVisited - vector of indexes of points to be ignored. (chosen to be closest to other points already).
	 *	output:	
	 *			The index of the closest point
	 */
	static int getClosestPoint(int currentPointIndex,vector<Point>& min_max_points, vector<int>& isVisited);

	/*	median
	*	input: a vector of numbers.
	*	output: median of the values in the vector.
	*/
	static double median(vector<double>& vec);

	/*	radiusOutOfPoints
	 *	input: a vector of points, and the center point.
	 *	output: median of the distances between a point from the vector and the center point.
	 */
	static double radiusOutOfPoints(vector<Point>& points, Point& center);


	/*	getMinMaxFlowerRatio
	*	input: minimum points on outer flower contour, radius of outer flower contour, and the center point.
	*	output: returns ratio between closest minimum point to the center and the outer contour radius.
	*/
	static double getMinMaxFlowerRatio(vector<Point> min_points, double outer_max_radius, Point& center);

	/*	contourToMatrix
	*	input: a contour, a two-dim allocated matrix.
	*	output: the contour points in a two-dim int matrix (which represents a binary image). contour points will have 1 in matrix.
	*/
	static void contourToMatrix(vector<Point> contour, Mat& arr);

	/*	contourToPointDst
	*	input: a contour and the center point.
	*	output: The distance between the contour and the center point -
	*			The distance between the closest point on the contour to the center point.
	*/
	static double contourToPointDst(vector<Point>& contour, Point& center);

	/*	intersectMasks
	*	input: Two masks.
	*	output: Intersect mask.
	*/
	static void intersectMasks(Mat& maskA, Mat& maskB, Mat& outMask);
};




#endif