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
};




#endif