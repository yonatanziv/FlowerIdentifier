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
};




#endif