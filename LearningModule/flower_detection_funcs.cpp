#include "flower_detection_header.h"
//version 58

void FlowerFeatureExtractor::convertRGBtoHSL (Scalar& rgb_color, Scalar& hsl_color){
	double r=rgb_color.val[0];
	double g=rgb_color.val[1];
	double b=rgb_color.val[2];
	double h,s,l, maxi, mini;
	double d;

	r /= 255.0, g /= 255.0, b /= 255.0;
	maxi = max(r, max(g, b));
	mini = min(r, min(g, b));

	l = (maxi + mini) / 2.0;

	if(maxi == mini){
		h = s = 0.0; // achromatic
	}else{
		d = maxi - mini;
		s = (l > 0.5) ? (d / (2.0 - maxi - mini)) : (d / (maxi + mini));
		if(maxi==r){h = (g - b) / d + ((g < b) ? 6.0 : 0.0);}
		else if(maxi==g){h = (b - r) / d + 2;}
		else if(maxi==b){h = (r - g) / d + 4;}
		h /= 6.0;
	}

	hsl_color.val[0] = h;
	hsl_color.val[1] = s;
	hsl_color.val[2] = l;

}

double FlowerFeatureExtractor::hue2rgb(double p, double q, double t){ //helper for convertHSLtoRGB
	double ret;
	if(t < 0) t += 1;
	if(t > 1) t -= 1;
	if(t < 1.0/6) return p + (q - p) * (6.0 * t);
	if(t < 1.0/2) return q;
	if(t < 2.0/3) return p + (q - p) * (2.0/3 - t) * 6.0;
	return p;
}

void FlowerFeatureExtractor::convertHSLtoRGB(Scalar& hsl_color, Scalar& rgb_color)
{
	double h = hsl_color.val[0];
	double s = hsl_color.val[1];
	double l = hsl_color.val[2];
	double r, g, b, q, p;

	if(s == 0){
		r = g = b = l; // achromatic
	}else{


		q = (l < 0.5) ? (l * (1 + s)) : (l + s - (l * s));
		p = 2 * l - q;
		r = hue2rgb(p, q, h + 1.0/3);
		g = hue2rgb(p, q, h);
		b = hue2rgb(p, q, h - 1.0/3);
	}
	rgb_color.val[0] = r*255;
	rgb_color.val[1] = g*255;
	rgb_color.val[2] = b*255;

}

int FlowerFeatureExtractor::square(int x)
{
	return x*x;
}

/*	distance
*	input: two points (x1,y1),(x2,y2)
*	output: euclidean distance between them
*/
double FlowerFeatureExtractor::distance(int y1, int x1, int y2, int x2)
{
	double dst = square(x1-x2)+square(y1-y2); 
	return sqrt(dst);
}

/*	distance
*	input: two CvPoint points
*	output: euclidean distance between them
*/
double FlowerFeatureExtractor::point_distance(Point& point1, Point& point2)
{
	double dst = square(point1.x-point2.x)+square(point1.y-point2.y);
	return sqrt(dst);
}

/*	padded_array_access
*	input: an int array, its size, and a cell in it (row and col)
*	output: returns image[row][col] if the indices are inside the array, 0 otherwise
*	used to prevent falling out of array in localMaxPoints and localMinPoints
*/
int FlowerFeatureExtractor::padded_matrix_access( Mat& image, int row, int col )
{
	if(row < 0 || col < 0 || col > image.cols-1 || row > image.rows-1)
		return 0;
	else
		return image.at<uchar>(row,col);

}

/*	localMinMaxPoints
*	input: an int array representing a binary image, its size, a point, an output array, resolution
*	output: outImage will have 2 at the cells where there are local maximum points of the distance function between the given point and the array point
*	not allows minimum points closer than <resolution> cells
*/
void FlowerFeatureExtractor::localMinMaxPointsInner(Mat& image, Point& center, Mat& outImage, int resolution, 
							vector<Point>& min_max_points, int min_max_flag)
{	
	double dst1, dst2;
	int dx,dy;
	int flag;
	for(int row = 0; row < image.rows; row++)
		for(int col = 0; col < image.cols; col++)
		{
			if(image.at<uchar>(row,col) == 1)
			{
				flag =0;
				dst1 = distance(row , col, center.y, center.x);
				for(dx = -resolution ; dx < resolution + 1; dx++)
					for(dy = -resolution; dy < resolution+1 ;dy++)
					{
						if(padded_matrix_access(image, row+dx, col+dy) == 1 && (dx!=0 || dy!=0))
						{
							dst2 = distance(row+dx, col+dy, center.y, center.x);
							if (min_max_flag == MAX_POINTS_FLAG) {
								if(dst2 > dst1)
									flag = 1;
							} else if (min_max_flag == MIN_POINTS_FLAG) {
								if(dst2 < dst1)
									flag =1;
							}

							if(padded_matrix_access(outImage, row+dx, col+dy) == 2) //in case of more than one max point at the same area.
								flag = 1;
						}
					}
					if(flag == 0) {
						outImage.at<uchar>(row,col) = 2;
						min_max_points.push_back(Point(col,row)); //col = x, row = y.
					}
			}				
		}
}

/*	localMaxPoints
*	input: an int array representing a binary image, its size, a point, an output array, resolution
*	output: outImage will have 2 at the cells where there are local maximum points of the distance function between the given point and the array point
*	not allows minimum points closer than <resolution> cells
*/
void FlowerFeatureExtractor::localMaxPoints(Mat& image, Point& center, Mat& outImage, int resolution, vector<Point>& max_points)
{	
	localMinMaxPointsInner(image, center, outImage, resolution, max_points, MAX_POINTS_FLAG);
}

/*	localMinPoints
*	input: an int array representing a binary image, its size, a point, an output array, resolution
*	output: outImage will have 2 at the cells where there are local minimum points of the distance function between the given point and the array point
*	not allows minimum points closer than <resolution> cells
*/
void FlowerFeatureExtractor::localMinPoints(Mat& image, Point& center, Mat& outImage, int resolution, vector<Point>& min_points)
{	
	localMinMaxPointsInner(image, center, outImage, resolution, min_points, MIN_POINTS_FLAG);
}

/*	drawArrPointsOnImage
*	input: an image, an array to draw, colors and such
*	output: draws the extremum points from the array on the image, with specified color, thickness, etc
*/
void FlowerFeatureExtractor::drawPointsOnImage(Mat& draw_img, Mat& contour_matrix_helper, int radius, 
					   const Scalar& color, int thickness, int line_type, int shift)
{
	Point now_point;
	for(int i=0; i<contour_matrix_helper.rows; i++){
		for(int j=0; j<contour_matrix_helper.cols; j++){
			if(contour_matrix_helper.at<uchar>(i,j) == 2){
				now_point.y=i;
				now_point.x=j;
				circle(draw_img, now_point, radius, color, thickness, line_type, shift);
			}
		}
	}
}

int FlowerFeatureExtractor::colorDst(Scalar& color1, Scalar& color2)
{	
	double dst = square(color1.val[0]-color2.val[0]) + square(color1.val[1]-color2.val[1]) + square(color1.val[2]-color2.val[2]);
	return sqrt(dst);
}

/* createCircleMask
*	input: an image (for mask output), radius, center
*	output: a binary image containing a filled circle, with specified radius, around center
*/
void FlowerFeatureExtractor::createCircleMask( Mat& img_mask, int radius, Point& center )
{
	int x,y;
	for(x=0; x<img_mask.cols; x++){
		for(y=0 ;y<img_mask.rows; y++){
			if(square(x-center.x) + square(y-center.y) > square(radius)){
				img_mask.at<uchar>(y,x) = 0;
			}
			else{
				img_mask.at<uchar>(y,x) = 1;
			}
		}
	}
}

/* flipMask
*	input: a binary mask
*	output: the binary mask, with 1's and 0's switched
*/
void FlowerFeatureExtractor::flipMask(Mat& img_mask )
{
	int x,y;
	for(x=0; x<img_mask.cols; x++){
		for(y=0; y<img_mask.rows; y++){
			if(img_mask.at<uchar>(y,x) == 1){
				img_mask.at<uchar>(y,x) = 0;
			}
			else{
				img_mask.at<uchar>(y,x) = 1;
			}
		}
	}
}

void FlowerFeatureExtractor::bgr2rgb(Scalar& bgr, Scalar& rgb) 
{
	rgb.val[0] = bgr.val[2];
	rgb.val[1] = bgr.val[1];
	rgb.val[2] = bgr.val[0];
}

/* toBinaryByDominantColorWithContour
*	input: a color image (3 colors * 8 bits depth), a pointer to output image (1 color, 8 bits depth), color threshold, contour pointer
*	output: img_output is a binary image, with white at all points inside contour (if contour is null, ignores contour) which are "closer" to dom_color than the color threshold
*/
void FlowerFeatureExtractor::toBinaryByDominantColorWithContour(Mat& img_input, Mat& img_output, Scalar& dom_color, int threshold, vector<Point>& contour)
{
	//colors only inside of contour. if contour is NULL, colors everywhere
	int x,y;
	Point2f point;
	Scalar rgb_color;

	for(x=0; x < img_input.cols; x++)
	{
		for(y=0; y < img_input.rows; y++)
		{
			Scalar bgr_color = img_input.at<Vec3b>(y,x);
			bgr2rgb(bgr_color, rgb_color);

			point.x=x;
			point.y=y;

			if(colorDst(rgb_color, dom_color) < threshold && (contour.empty() || pointPolygonTest(contour, point, false) > 0 )) {
				img_output.at<uchar>(y,x) = 255; //White to greater of threshold
			} else{
				img_output.at<uchar>(y,x) = 0; //Black other
			}
		}
	}
	//showAndWait("moo",img_output);
}

double FlowerFeatureExtractor::contourToPointDst(vector<Point>& contour, Point& center)
{
	double minDst = -1;
	double dst =0;

	for(int i=0; i < contour.size(); ++i ){	
		dst = point_distance(contour.at(i), center);
		if(minDst == -1 || dst < minDst)
			minDst = dst;
	}
	return minDst;
}

/* findMaxContour
*	input: a pointer to a list of contours (first_contour), a pointer to a pointer of a contour
*	output: max_contour will contain the contour from the list with the largest amount of points on the contour
*/
void FlowerFeatureExtractor::findMaxContour(vector<vector<Point> >& contours, vector<Point>& max_contour, Mat& img_draw_screen)
{
	int maxsize = 0;
	int maxIndex = 0;

	max_contour = vector<Point>(); //return NULL if no contour
	for(int i = 0; i < contours.size(); i++) { //find maximal contour
		drawContours(img_draw_screen, contours, i, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
		//showAndWait(argv[0],img_draw_screen);

		if(contours.at(i).size() > max_contour.size()) {
			max_contour = contours.at(i);
		}
	}
	//showAndWait("lol",img_draw_screen);
}

/* findClosestToContour
*	input: a pointer to a list of contours (first_contour), a pointer to a pointer of a contour
*	output: contour will contain the contour from the list which is closest to the center point.
*/
void FlowerFeatureExtractor::findClosestToCenterPointContour(vector<vector<Point> >& contours, vector<Point>& inner_contour, Point& center, Mat& img_draw_screen)
{
	double minDst = -1;
	double dst = 0;
	inner_contour = vector<Point>();

	for(int i = 0; i < contours.size(); i++) { //find maximal contour
		drawContours(img_draw_screen, contours, i, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
		//showAndWait(argv[0],img_draw_screen);

		if(contours.at(i).size() >= MINIMUM_INNER_CONTOUR_POINTS) {
			dst = contourToPointDst(contours.at(i), center);
			if(minDst == -1 || dst < minDst) {
				minDst = dst;
				inner_contour = contours.at(i);
			}
		}
	}
	//showAndWait("lol2",img_draw_screen);
}

/* contourToArray
*	input: a contour, a two-dim allocated array
*	output: the contour points in a two-dim int array (which represents a binary image). contour points will have 1 in array
*/
void FlowerFeatureExtractor::contourToMatrix(vector<Point> contour, Mat& arr)
{
	for( int i=0; i<contour.size(); ++i ) {
		arr.at<uchar>(contour.at(i).y, contour.at(i).x) = 1;
	}
}

/* applyMeanShift
*	input: an image
*	output: the image, after applying the meanshift algorithm on it. contains a fix to the size of the image, to allow the meanshift to work
*/
void FlowerFeatureExtractor::applyMeanShift( Mat& img_meanshift )
{
	//Here comes the thing (fix for meanshift)
	img_meanshift.cols &= -(1<<MEANSHIFT_LEVEL);
	img_meanshift.rows &= -(1<<MEANSHIFT_LEVEL);
	pyrMeanShiftFiltering(img_meanshift, img_meanshift,MEANSHIFT_SPATIAL_RADIUS,MEANSHIFT_COLOR_RADIUS,MEANSHIFT_LEVEL);
}

/* createAndCalcHistogram
*	input: an image and a mask
*	output: a pointer to an allocated histogram of the image area inside img_mask
*/
CvHistogram* FlowerFeatureExtractor::createAndCalcHistogram(IplImage* img_hist, IplImage* img_mask){
	CvHistogram* hist;
	int numBins[3] = {NUM_BINS_FOR_COLOR,NUM_BINS_FOR_COLOR,NUM_BINS_FOR_COLOR};
	float range[] = {0, 256};
	float *ranges[] = {range,range,range};
	hist = cvCreateHist(3, numBins, CV_HIST_ARRAY, ranges, 1);
	cvClearHist(hist);

	IplImage* imgRed = cvCreateImage(cvGetSize(img_hist), 8, 1);
	IplImage* imgGreen = cvCreateImage(cvGetSize(img_hist), 8, 1);
	IplImage* imgBlue = cvCreateImage(cvGetSize(img_hist), 8, 1);
	IplImage* img_channels[3];
	cvSplit(img_hist, imgBlue, imgGreen, imgRed, NULL);
	img_channels[0]=imgRed;
	img_channels[1]=imgGreen;
	img_channels[2]=imgBlue;
	cvCalcHist(img_channels, hist, 0, img_mask);

	cvReleaseImage(&imgRed);
	cvReleaseImage(&imgGreen);
	cvReleaseImage(&imgBlue);
	return hist;
}

/*	median
*	input: an array and its size
*	output: median of the values in the array
*/
double FlowerFeatureExtractor::median(vector<double>& vec)
{
	if(vec.empty()) 
		return 0;
	else {
		sort(vec.begin(), vec.end());
		if(vec.size() % 2 == 0)
			return (vec.at(vec.size()/2 - 1) + vec.at(vec.size()/2)) / 2;
		else
			return vec.at(vec.size()/2);
	}
}

/* radiusOutOfPoints
input: an array of points, its size, a center point
output: median of the distances between a point from the array and the center point
*/
double FlowerFeatureExtractor::radiusOutOfPoints(vector<Point>& max_points, Point& center)
{
	double result;
	vector<double> distances;
	for(int i =0; i < max_points.size(); i++) {
		distances.push_back(point_distance(max_points.at(i), center));
	}
	result = median(distances);
	return result;
}

/*	getMinMaxFlowerRatio
*	input: minimum points on outer flower contour, radius of outer flower contour
*	output: returns ratio between closest minimum point to the center and the outer contour ratio
*/
double FlowerFeatureExtractor::getMinMaxFlowerRatio(vector<Point> min_points, double radius_max, Point& center)
{
	double min_distance;
	vector<double> distances;
	for(int i =0; i < min_points.size(); i++) {
		distances.push_back(point_distance(min_points.at(i), center));
	}

	for(int i=0; i < min_points.size(); i++){
		if(i==0)
			min_distance = distances.at(i);
		else if(min_distance > distances.at(i))
			min_distance = distances.at(i);
	}

	return min_distance / radius_max;
}

/* contourVar
*	input: a contour and flower's center
*	output: returns variance of distances between a point on the contour and the center point
*/
double FlowerFeatureExtractor::contourVar(vector<Point>& contour, Point& center) {
	double Ex = 0; //E[x]
	double Ex2 = 0;//E[x^2]
	int i;

	for( i=0; i < contour.size(); ++i ){	
		Ex += point_distance(contour.at(i), center);
		Ex2 += square(contour.at(i).x-center.x) + square(contour.at(i).y-center.y);
	}

	Ex /= i;
	Ex2 /= i;

	return Ex2 - Ex*Ex;
}

/*avgDst
*input: contour and center point.
*output: The average distance between the center point and the contour. 
*/
double FlowerFeatureExtractor::avgDst(vector<Point>& contour, Point& center){
	double radius = 0;
	for(int i=0; i < contour.size(); ++i ){	
		radius += point_distance(contour.at(i), center);
	}
	radius /= contour.size();
	return radius;
}


/*	getOptimalInnerContour
*	input: max_flower_conotur (contour of outer flower), pointer to max_inner_contour pointer (to update pointer), dominant colors, 
circle_around_center_ponit_flag - this flag will be activated (by the calling function) if no inner contour found the first time this function is run. if this flag is on, the inner contour will be a circle around the center point with radius of the average distance between the center point and the contour with the smallest variance.
radius_max in order to elimnate bad contours - in case of the radius of the circle around center point is too big (more than half of the radius_max).
*	output: runs on color thresholds (according to OPTIMAL_INNER_CONTOUR_START, OPTIMAL_INNER_CONTOUR_END, OPTIMAL_INNER_CONTOUR_QUOTIENT,
*			paints the are which is within <threshold> distance from the center color, and within the outer flower contour, and tries to find contour there.
picks the contour with smallest variance of distances between a point on the contour and the center point
*	fail: max_inner_contour=NULL if no inner contour found.
*/
void FlowerFeatureExtractor::getOptimalInnerContour(vector<Point>& max_flower_contour, vector<Point>& max_inner_contour, Scalar& dom_flower_color,
							Scalar& dom_inner_color, int circle_around_center_ponit_flag, Mat& img_meanshift, 
							Mat& img_draw_screen, Point& center, double radius_max )
{

	Mat img_binary_small = Mat(img_meanshift.size(), CV_8UC1);
	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	vector<Point> optimal_contour, current_inner_contour;
	double min_variance, current_variance;
	Point2f center2f = Point2f(center.x, center.y);
	int colorDistance;
	int threshold_small;
	char name[100];

	for(int i=OPTIMAL_INNER_CONTOUR_START; i<OPTIMAL_INNER_CONTOUR_END; i++){
		colorDistance = colorDst(dom_flower_color, dom_inner_color);
		threshold_small = colorDistance * ((i*1.0) / OPTIMAL_INNER_CONTOUR_QUOTIENT);
		threshold_small = threshold_small < INNER_CONTUOR_MINIMUM_THRESHOLD ? INNER_CONTUOR_MINIMUM_THRESHOLD : threshold_small;
		//printf("SMALL THRESHOLD - %d",threshold_small);

		toBinaryByDominantColorWithContour(img_meanshift, img_binary_small, dom_inner_color, threshold_small, max_flower_contour);
		sprintf(name,"smallbw %d",i);
		Utils::myShowImage(name, img_binary_small);

		findContours(img_binary_small, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
		if(circle_around_center_ponit_flag == 0)
			findMaxContour(contours, current_inner_contour, img_draw_screen);
		else
			findClosestToCenterPointContour(contours, current_inner_contour, center, img_draw_screen);

		if(current_inner_contour.empty()){	
			if(colorDistance <= INNER_CONTOUR_EQUAL_FLOWER_CONTOUR_THRESHOLD){	//if there is no inner contour the inner and outer contours are identical.
				max_inner_contour=max_flower_contour;
				return;
			}
			continue;
		}

		if(circle_around_center_ponit_flag == 0 && pointPolygonTest(current_inner_contour, center2f, 0) <= 0) { //the last expression checks that the center point is inside the contour. 
			continue;
		}

		if(circle_around_center_ponit_flag == 1 && current_inner_contour.size() <= MINIMUM_INNER_CONTOUR_POINTS){
			continue;
		}

		current_variance = contourVar(current_inner_contour, center);
		if(optimal_contour.empty() || min_variance > current_variance) {
			min_variance = current_variance;
			optimal_contour = current_inner_contour;
		}
		img_meanshift.copyTo(img_draw_screen);
		vector<vector<Point> > current_inner_contour_vec;
		current_inner_contour_vec.push_back(current_inner_contour);
		drawContours(img_draw_screen, current_inner_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
		//showAndWait("Opencv_test2.exe", img_draw_screen );
	}

	if(circle_around_center_ponit_flag == 1 && !optimal_contour.empty()) {

		Mat circle_image = Mat(img_meanshift.size(), CV_8UC1);
		double radius = avgDst(optimal_contour, center);

		if(radius <= INNER_RADIUS_TO_FLOWER_RADIUS_RATIO * radius_max) {
			createCircleMask(circle_image,radius,center);
			findContours(circle_image, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
			findMaxContour(contours, optimal_contour, img_draw_screen);
		}
		else {
			optimal_contour = vector<Point>();
		}
	}

	max_inner_contour = optimal_contour;
}

double FlowerFeatureExtractor::innerProduct(Point& v1, Point& v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}

/*	angle
*	input: gets three points
*	output: calculates angle between the three points
*/
double FlowerFeatureExtractor::angle(Point& center, Point& p1, Point& p2)
{
	Point v1,v2;
	v1.x = p1.x-center.x;
	v1.y = p1.y-center.y;

	v2.x = p2.x-center.x;
	v2.y = p2.y-center.y;

	return abs(acos(innerProduct(v1,v2)/sqrt(innerProduct(v1,v1)*innerProduct(v2,v2))));
}

/*
*	input:	- one point of the array
- array of points
- isVisited - the index of points to be ignored. (chosen to be closest to other points already).
- currentPointIndex - the index of the above point. its also need to be ignored.
output:	- the index of the closest point
*/
int FlowerFeatureExtractor::getClosestPoint(int currentPointIndex,vector<Point>& min_max_points, vector<int>& isVisited)
{
	double dst, min = -1;
	int minIndex;
	Point minPoint;
	Point point = min_max_points.at(currentPointIndex);

	for(int i = 0; i <min_max_points.size(); i++)
	{		
		if(i == currentPointIndex || isVisited.at(i)==1)
			continue;
		dst = point_distance(point, min_max_points.at(i));
		if(min == -1 || dst < min) {
			min = dst;
			minPoint = min_max_points.at(i);
			minIndex = i;
		}
	}
	return minIndex;
}

/*
*	input:	- array of point (min or max points) in size greater than 2.
*			- center point
*	output:	The median of all angles between two close points and the center.
*/
double FlowerFeatureExtractor::getMedianAngle(vector<Point>& min_max_points, Point& center, Mat& draw_img, 
					  double& outMinAngle, int* outMinAngleTwoPointsIndexes)
{
	if(min_max_points.size() < 2)
		throw MinMaxPointsNumberException();

	double result;
	vector<double> angles(min_max_points.size(),0);
	vector<int> isVisited(min_max_points.size(),0); //contain ponits which have been chosen to be closest to other point.		
	int closestPointIndex;
	int currentPointIndex;
	int i;
	double minAngle=361;	//if outMinAnglePointIndex!=NULL it will contain the min angle. notice:we work with radian so 361 is arbitrary big number.

	currentPointIndex = 0;
	for(i=0; i < min_max_points.size()-1; i++)	//we run until min_max_points.size()-1 because for the last run isVisited is full. we will take care of the last run seperately.
	{	
		closestPointIndex = getClosestPoint(currentPointIndex, min_max_points, isVisited);
		isVisited.at(currentPointIndex) = 1;
		angles.at(i) = angle(center, min_max_points.at(currentPointIndex), min_max_points.at(closestPointIndex));

		if(outMinAngleTwoPointsIndexes != NULL){
			if(angles.at(i) < minAngle){
				minAngle = angles.at(i);
				outMinAngle = minAngle;
				outMinAngleTwoPointsIndexes[0] = currentPointIndex;
				outMinAngleTwoPointsIndexes[1] = closestPointIndex;
			}
		}

		circle(draw_img, min_max_points.at(currentPointIndex), 2, CV_RGB(0xff,0x00,0x00), -1, 8, 0);
		//showAndWait("test",draw_img);
		//printf("%lf \n",angles[i]*180/3.14159);*/
		currentPointIndex = closestPointIndex;
	}

	if(i == min_max_points.size() - 1)
	{
		angles.at(i) = angle(center, min_max_points.at(currentPointIndex), min_max_points.at(0));

		if(outMinAngleTwoPointsIndexes != NULL){
			if(angles.at(i) < minAngle){
				minAngle = angles.at(i);
				outMinAngle = minAngle;				
				outMinAngleTwoPointsIndexes[0] = currentPointIndex;
				outMinAngleTwoPointsIndexes[1] = 0;
			}
		}

		//printf("%lf \n",angles[i]*180/3.14159);
		circle(draw_img, min_max_points.at(currentPointIndex), 2, CV_RGB(0xff,0x00,0x00), -1, 8, 0);
		//showAndWait("test",draw_img);
	}

	result = median(angles);

	return result;
}

/*	toDegree
*	input: angle in radians
*	output: angle in degrees
*/
double FlowerFeatureExtractor::toDegree(double medAngle)
{
	return medAngle*180.0/PI;
}

/*	getFixBadMinMaxPointsThreshold
*	input: median of angles between min or max points
*	output: threshold, according to medAngle (big, medium, or small threshold)
*/
int FlowerFeatureExtractor::getFixBadMinMaxPointsThreshold(double medAngle)
{
	double angle = toDegree(medAngle);
	if(angle <= FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL_ANGLES)
		return FIX_BAD_MIN_MAX_POINTS_THRESHOLD_SMALL;
	if(angle <= FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG_ANGLES)
		return FIX_BAD_MIN_MAX_POINTS_THRESHOLD_MEDIUM;
	return FIX_BAD_MIN_MAX_POINTS_THRESHOLD_BIG;

}

/*
*	This function removes bad max points - if the angle between two max points is smaller 
*	than the medianAngle/fix_bad_min_max_points_threshold, it is removed.
*	
*	output:	The function update pointsArr and return its new length.
*			- ouMedAngle - will contain the new median angle. 
*/
int FlowerFeatureExtractor::fixBadMinMaxPoints(vector<Point>& min_max_points, Point& center, Mat& draw_img, double& outMedAngle, Mat& img_meanshift)
{
	double medAngle;
	double minAngle;
	int minPointsIndex[2];	//outMinAngleTwoPointsIndexes parameter for getMedianAngle.
	vector<int> isVisited(min_max_points.size(), 0);
	CvPoint* tempPointsArr;
	int secondClosestPoinIndex;
	double dst1,dst2;
	int i,j;
	int removePointIndex;
	int fix_bad_min_max_points_threshold;

	medAngle = getMedianAngle(min_max_points, center, draw_img, minAngle, minPointsIndex);
	if(medAngle <= 0) {
		throw MedianAngleException();
	}

	fix_bad_min_max_points_threshold = getFixBadMinMaxPointsThreshold(medAngle);

	while(minAngle < medAngle/fix_bad_min_max_points_threshold)
	{		
		for(i=0; i < min_max_points.size(); i++)
			isVisited.at(i)=0;

		isVisited.at(minPointsIndex[1]) = 1;

		secondClosestPoinIndex = getClosestPoint(minPointsIndex[0], min_max_points, isVisited);
		dst1 = point_distance(min_max_points.at(minPointsIndex[0]), min_max_points.at(secondClosestPoinIndex));

		// printf("\n dst1 == %lf\n",dst1);
		// printf("\n point1 - x=%d y=%d\n",(*pointsArr)[minPointsIndex[0]].x,(*pointsArr)[minPointsIndex[0]].y);

		isVisited.at(minPointsIndex[1]) = 0;
		isVisited.at(minPointsIndex[0]) = 1;
		secondClosestPoinIndex = getClosestPoint(minPointsIndex[1],min_max_points, isVisited);
		dst2 = point_distance(min_max_points.at(minPointsIndex[1]), min_max_points.at(secondClosestPoinIndex));

		//printf("\n dst2 == %lf\n",dst2);
		// printf("\n point2 - x=%d y=%d\n",(*pointsArr)[minPointsIndex[1]].x,(*pointsArr)[minPointsIndex[1]].y);

		removePointIndex = dst1<dst2?	minPointsIndex[0] :	minPointsIndex[1];
		min_max_points.erase(min_max_points.begin() + removePointIndex);

		img_meanshift.copyTo(draw_img);

		medAngle = getMedianAngle(min_max_points, center, draw_img, minAngle, minPointsIndex);
		if(medAngle <=0) {
			throw MedianAngleException();
		}
		fix_bad_min_max_points_threshold = getFixBadMinMaxPointsThreshold(medAngle);
	}

	outMedAngle = medAngle;
	return min_max_points.size();
}


/*	createCircleMasks
*	input: binary mask (IplImage) pointers, mask_flag (changes inner mask radius, by steps same as COLOR_FLAG_STEP), background_mask_flag (changes background mask radius, by steps same as COLOR_FLAG_STEP)
*	output: creates the required inner part, flower, and background binary masks according to flags
*/
void FlowerFeatureExtractor::createCircleMasks( Mat& img_mask_small, Mat& img_mask_large, Mat& img_mask_background,Point& center,int mask_flag,int background_mask_flag )
{
	int radius_large = center.x > center.y ? center.y : center.x;	
	int radius_background = radius_large;

	radius_large*=1 - (COLOR_FLAG_STEP * mask_flag);
	radius_background*= 1 + (COLOR_FLAG_STEP * background_mask_flag);

	createCircleMask(img_mask_small, RADIUS_SMALL, center);
	createCircleMask(img_mask_large, radius_large, center);
	createCircleMask(img_mask_background, radius_background, center);	
	flipMask(img_mask_background);
}

/*	getDominantColors
*	input: dominant color struct pointers, mask_flag (changes inner mask radius), background_mask_flag (changes background mask radius)
*	output: dominant colors, according to the flags given
*/
void FlowerFeatureExtractor::getDominantColors(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color, 
					   Mat& img, Point& center,int mask_flag, int background_mask_flag)
{

	/* CREATE MASK */
	Mat img_mask_small = Mat(img.size(), CV_8UC1);//inner circle mask
	Mat img_mask_large = Mat(img.size(), CV_8UC1); //whole flower mask 
	Mat img_mask_background = Mat(img.size(), CV_8UC1); //background mask

	createCircleMasks(img_mask_small,img_mask_large,img_mask_background,center,mask_flag, background_mask_flag);

	/* END MASK */
	
	CvHistogram *hist_small = createAndCalcHistogram(&(IplImage)img, &(IplImage)img_mask_small);	//inner circle histogram
	CvHistogram *hist_large = createAndCalcHistogram(&(IplImage)img, &(IplImage)img_mask_large);	//whole flower histogram
	CvHistogram *hist_background = createAndCalcHistogram(&(IplImage)img, &(IplImage)img_mask_background); //background histogram

	/* GET DOMINANT COLOR */

	float histMaxInner = 0;
	int maxIndexInner[3];
	cvGetMinMaxHistValue(hist_small, 0, &histMaxInner, 0, maxIndexInner);

	float histMaxFlower = 0;
	int maxIndexFlower[3];
	cvGetMinMaxHistValue(hist_large, 0, &histMaxFlower, 0, maxIndexFlower);

	float histMaxBackground = 0;
	int maxIndexBackground[3];
	cvGetMinMaxHistValue(hist_background, 0, &histMaxBackground, 0, maxIndexBackground);	

	/* END DOMINANT COLOR */

	dom_inner_color=Scalar((maxIndexInner[0]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexInner[1]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexInner[2]*256)/NUM_BINS_FOR_COLOR);

	dom_flower_color=Scalar((maxIndexFlower[0]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexFlower[1]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexFlower[2]*256)/NUM_BINS_FOR_COLOR);

	dom_background_color=Scalar((maxIndexBackground[0]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexBackground[1]*256)/NUM_BINS_FOR_COLOR,
		(maxIndexBackground[2]*256)/NUM_BINS_FOR_COLOR);

	cvReleaseHist(&hist_small);
	cvReleaseHist(&hist_large);
	cvReleaseHist(&hist_background);
}



/*	getThresholdByFlag
*	input: color_flag (changes threshold which seperates background and flower color )
*	output: the threshold which seperates background and flower color, according to the flag
*/
int FlowerFeatureExtractor::getThresholdByFlag( int color_flag, Scalar& dom_flower_color, Scalar& dom_background_color )
{
	return colorDst(dom_flower_color,dom_background_color)*(OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR + color_flag)/OUTER_COLOR_DST_THRESHOLD_DENOMINATOR;
}

/* checkBinaryBackground
*	input: binary image, background_mask_flag (changes background mask radius)
*	output: "masks" the image with background mask, and checks if there are more black pixels than white pixels in the image
*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black pixels in that area (what would indicate a bad binarization - too much background)
*/
void FlowerFeatureExtractor::checkBinaryBackground( Mat& img_binary, int background_mask_flag, Point& center )
{
	Mat img_mask_background = Mat(img_binary.size(), CV_8UC1);
	Mat img_binary_after_mask = Mat(img_binary.size(), CV_8UC1);

	int radius_background = center.x > center.y ? center.y : center.x;	

	radius_background*= 1 + (COLOR_FLAG_STEP * background_mask_flag);

	createCircleMask(img_mask_background, radius_background, center);
	flipMask(img_mask_background);

	img_binary.copyTo(img_binary_after_mask, img_mask_background);

	int x, y, black_count = 0, white_count = 0;
	//showAndWait("bin mask",img_binary_after_mask);
	for(x=0; x < img_binary_after_mask.cols; x++){
		for(y=0; y < img_binary_after_mask.rows; y++){
			if(square(x-center.x) + square(y-center.y) <= square(radius_background))
				continue;
			if(img_binary_after_mask.at<uchar>(y,x) == 0) {
				black_count++;

			}
			else{
				white_count++;
			}
		}
	}

	if(white_count > black_count)
		throw ContourFailBinaryBackgroundException();
}

/*	getOptimalFlowerContour
*	input: pointer to pointer of outer flower contour (in order to update pointer), dominant colors, color flag (changes color threshold), background mask flag (changes background mask radius)
*	output: gets optimal (largest) contour of outer flower (hopefully), according to the parameters.
*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black in the area masked by background mask (fail of checkBinaryBackground)
: CONTOUR_FAIL_CENTER_POINT if center point is not inside the contour
*/
int FlowerFeatureExtractor::getOptimalFlowerContour(vector<Point>& max_flower_contour, Mat& img, Scalar& dom_flower_color, Scalar& dom_background_color, 
							Point& center, int color_flag, int background_mask_flag )
{
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	Mat img_binary = Mat(img.size(), CV_8UC1);
	Mat img_draw_screen;

	int threshold_large=getThresholdByFlag(color_flag,dom_flower_color,dom_background_color);

	//printf("using threshold_large: %d\n",threshold_large);
	try 
	{
		toBinaryByDominantColorWithContour(img, img_binary, dom_flower_color, threshold_large, vector<Point>());
		checkBinaryBackground(img_binary, background_mask_flag, center);

		Utils::myShowImage("screen5", img_binary);

		findContours(img_binary, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
		if (contours.size() == 0) {
			return CONTOUR_FAIL_NO_POINTS;
		}

		img.copyTo(img_draw_screen);

		findMaxContour(contours, max_flower_contour, img_draw_screen);

		Point2f center2f;
		center2f.x=center.x;
		center2f.y=center.y;
		if(pointPolygonTest(max_flower_contour, center2f, 0) <= 0)
			return CONTOUR_FAIL_CENTER_POINT;

		return SUCCESS;
		
	} catch (ContourFailBinaryBackgroundException& e) {
		return CONTOUR_FAIL_BINARY_BACKGROUND;
	} 
	
}

/*	getOuterFlowerProperties
*	input: outer contour of flower, and parameters to update
*	update: updates sample parameters related to outer part of flower: num of max and min points, median of angles between min points and max points, ratio between closest min point to the center and the flower radius
*/
void FlowerFeatureExtractor::getOuterFlowerProperties(int& num_points_max, int& num_points_min, double& angle_max, double& angle_min, 
							  double& min_max_flower_ratio, double& radius_max, Point& center, 
							  vector<Point> max_flower_contour, Mat& img_meanshift, Mat& img_draw_screen)
{
	vector<Point> max_points;		//Will hold the flower max points.
	vector<Point> min_points;		//Will hold the flower min points.
	int cols = img_meanshift.cols;
	int rows = img_meanshift.rows;
	Mat contour_matrix = Mat::zeros(img_meanshift.size(), CV_8UC1);
	Mat contour_matrix_helper = Mat::zeros(img_meanshift.size(), CV_8UC1);

	contourToMatrix(max_flower_contour, contour_matrix);
	contour_matrix.copyTo(contour_matrix_helper);

	localMaxPoints(contour_matrix, center, contour_matrix_helper, MAX_POINTS_THRESHOLD, max_points);
	drawPointsOnImage(img_draw_screen, contour_matrix_helper, 2, CV_RGB(0x00,0x00,0x00), -1, 8, 0);
	if(max_points.size() == 0)
		throw NoOuterMaxPointsFoundException();

	radius_max = radiusOutOfPoints(max_points, center);

	Mat img_with_contours = Mat(img_meanshift.size(), img_meanshift.type());
	img_draw_screen.copyTo(img_with_contours);
	
	num_points_max = fixBadMinMaxPoints(max_points, center, img_draw_screen, angle_max, img_with_contours);
	contour_matrix.copyTo(contour_matrix_helper);
	localMinPoints(contour_matrix, center, contour_matrix_helper, MIN_POINTS_THRESHOLD, min_points);
	drawPointsOnImage(img_draw_screen, contour_matrix_helper, 2, CV_RGB(0x00,0x00,0xff), -1, 8, 0);
	Utils::myShowImage("screen6", img_draw_screen);

	if(min_points.size() == 0)
		throw NoOuterMinPointsFoundException(); // NO_OUTER_MIN_POINTS_FOUND_ERROR;


	min_max_flower_ratio = getMinMaxFlowerRatio(min_points, radius_max, center);

	img_meanshift.copyTo(img_draw_screen);
	vector<vector<Point> > max_flower_contour_vec;
	max_flower_contour_vec.push_back(max_flower_contour);

	drawContours(img_draw_screen, max_flower_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
	drawPointsOnImage(img_draw_screen, contour_matrix_helper, 2, CV_RGB(0x00,0x00,0x00), -1, 8, 0);

	img_draw_screen.copyTo(img_with_contours);

	num_points_min = fixBadMinMaxPoints(min_points, center, img_draw_screen, angle_min, img_with_contours); //updates angle_min
}

/*	getInnerPartProperties
*	input: outer contour of the flower(in order to create a mask, to define where we look for the inner contour), radius_inner to update, colors and so on (radius_max in order to elimnate bad contours)
*	output: updates sample properties related to the inner part of the flower. for now: radius_inner
*/
void FlowerFeatureExtractor::getInnerPartProperties(double& radius_inner, double& inner_length, vector<Point>& max_flower_contour, 
							Scalar& dom_flower_color, Scalar& dom_inner_color, Point& center, Mat& img, 
							Mat& img_draw_screen, double radius_max)
{
	vector<Point> max_inner_contour;
	int cols = img.cols;
	int rows = img.rows;
	int num_points_max_inner, num_points_min_inner;

	getOptimalInnerContour(max_flower_contour, max_inner_contour, dom_flower_color, dom_inner_color, 0, 
						   img, img_draw_screen, center, radius_max);
	if(max_inner_contour.empty())
		getOptimalInnerContour(max_flower_contour, max_inner_contour, dom_flower_color,dom_inner_color, 1,
							   img,img_draw_screen,center,radius_max);
	if(max_inner_contour.empty())
		throw MaxInnerContourException();// NO_INNER_CONTOUR_FOUND_ERROR;

	/* ALLOCATE ARRS */
	Mat out_arr_max = Mat(img.size(), CV_8UC1);
	Mat out_arr_min = Mat(img.size(), CV_8UC1);
	Mat contour_arr = Mat(img.size(), CV_8UC1);
	vector<Point> max_points_inner, min_points_inner, combined_points;

	/* END ALLOC ARRS */
	inner_length = arcLength(max_inner_contour, true);
	//printf("Inner contour Length: %lf, Area: %lf, Ratio: %lf\n", cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED), cvContourArea(max_inner_contour),(cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED)*cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED))/ (cvContourArea(max_inner_contour)));
	img.copyTo(img_draw_screen);
	vector<vector<Point> > max_inner_contour_vec;
	max_inner_contour_vec.push_back(max_inner_contour);
	drawContours(img_draw_screen, max_inner_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens

	contourToMatrix(max_inner_contour, contour_arr);

	contour_arr.copyTo(out_arr_max);
	contour_arr.copyTo(out_arr_min);

	localMaxPoints(contour_arr, center, out_arr_max, MAX_POINTS_THRESHOLD, max_points_inner);
	drawPointsOnImage(img_draw_screen, out_arr_max, 2, CV_RGB(0xff,0x00,0x00), -1, 8, 0);

	localMinPoints(contour_arr, center, out_arr_min, MIN_POINTS_THRESHOLD, min_points_inner);
	drawPointsOnImage(img_draw_screen, out_arr_min, 2, CV_RGB(0x00,0x00,0xff), -1, 8, 0);
	Utils::myShowImage("screen7", img_draw_screen);

	if(max_points_inner.empty() || min_points_inner.empty())
		throw NoInnerMinMaxPointsException(); //NO_INNER_MIN_MAX_POINTS_FOUND_ERROR;

	for(int i = 0; i < max_points_inner.size(); i ++) {
		combined_points.push_back(max_points_inner.at(i));
	}

	for(int i = 0; i < min_points_inner.size(); i ++) {
		combined_points.push_back(min_points_inner.at(i));
	}

	radius_inner = radiusOutOfPoints(combined_points, center);

}

/*	getFlowerContourAndColor
*	input: pointer to pointer of max_flower_contour (in order to update it), color struct pointers
*	output: updates outer flower contour and dominant colors of the flower, inner, and background.
*	may try several background and inner mask radiuses, and several color thresholds, if it gets bad contours (no points, or center point not inside contour)
*/
void FlowerFeatureExtractor::getFlowerContourAndColor(vector<Point>& max_flower_contour, Mat& img, Scalar& dom_inner_color, Scalar& dom_flower_color, 
							  Scalar& dom_background_color, Point& center)
{
	int mask_flag, color_flag, contour_fail, background_mask_flag;
	int mask_end = OUTER_COLOR_DST_THRESHOLD_MAX_NUMERATOR - OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR;
	Scalar dom_color_temp_hsl;
	/*LIGHTNESS FIX*/
	Scalar dom_inner_color_outer_contour, dom_flower_color_outer_contour, dom_background_color_outer_contour;
	Mat img_outer_contour;
	Mat img_draw_screen;
	vector<vector<Point> > max_flower_contour_vec;

	img.copyTo(img_outer_contour);
	equalizeLightness(img_outer_contour);
	img_outer_contour.convertTo(img_outer_contour, -1, CONTRAST_RATIO_OUTER_CONTOUR); //contrast
	Utils::myShowImage("lightness",img_outer_contour);

	for(background_mask_flag=0; background_mask_flag < BACKGROUND_MASK_FLAG_ITERATIONS; background_mask_flag++) { //if no contour at all is found, probably background color and flower color are the same. then, try to increase background mask radius to fix this
		getDominantColors(dom_inner_color,dom_flower_color,dom_background_color,img,center,0, background_mask_flag);
		
		getDominantColors(dom_inner_color_outer_contour, dom_flower_color_outer_contour, dom_background_color_outer_contour, img_outer_contour, center, 0, background_mask_flag);
		convertRGBtoHSL(dom_flower_color, dom_color_temp_hsl); //for printing color
		
		if(isGrayscale(dom_flower_color)){ //don't fix lightness
			contour_fail = getOptimalFlowerContour(max_flower_contour, img, dom_flower_color, dom_background_color,
				center, 0, 0);
		}

		else{ //lightness fix
			
			contour_fail = getOptimalFlowerContour(max_flower_contour, img_outer_contour, dom_flower_color_outer_contour,
				dom_background_color_outer_contour, center, 0, 0);
			//printf("DOM COLOR1: H:%.2lf S:%.2lf L:%.2lf ",dom_color_temp_hsl.val[0],dom_color_temp_hsl.val[1],dom_color_temp_hsl.val[2]);
			//printf("first, lightness fix APPLIED\n");
		}

		img.copyTo(img_draw_screen);
		max_flower_contour_vec.push_back(max_flower_contour);
		drawContours(img_draw_screen, max_flower_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
		
		//char txt[20];
		//sprintf(txt,"flower contour %d",background_mask_flag);
		//myShowImage(txt, img_draw_screen);
		if(contour_fail != CONTOUR_FAIL_NO_POINTS)	//Here we care only about that error, because we want to make sure that we got the right background dominant color. Therefore we dont check CONTOUR_FAIL_BINARY_BACKGROUND
			break;
	}
	for(color_flag=0; color_flag < COLOR_FLAG_ITERATIONS; color_flag++) { //runs on color thresholds
		for(mask_flag=0; mask_flag < mask_end; mask_flag++){ //runs on mask radius	
			getDominantColors(dom_inner_color, dom_flower_color, dom_background_color, 
				img, center, mask_flag, background_mask_flag); //mask_flag to set radius for mask

			getDominantColors(dom_inner_color_outer_contour, dom_flower_color_outer_contour, 
				dom_background_color_outer_contour, img_outer_contour, center, 0, background_mask_flag);

			convertRGBtoHSL(dom_flower_color, dom_color_temp_hsl);
			if(isGrayscale(dom_flower_color)) { //don't fix lightness
				contour_fail = getOptimalFlowerContour(max_flower_contour, img, dom_flower_color, dom_background_color, 
													   center, color_flag, background_mask_flag);
			}
			else {//lightness fix
				contour_fail = getOptimalFlowerContour(max_flower_contour, img_outer_contour, dom_flower_color_outer_contour, 
													   dom_background_color_outer_contour, center, color_flag, 
													   background_mask_flag);
				//printf("DOM COLOR2: H:%.2lf S:%.2lf L:%.2lf ",dom_color_temp_hsl.val[0],dom_color_temp_hsl.val[1],dom_color_temp_hsl.val[2]);
				//printf("second, lightness fix APPLIED\n");
			}

			img.copyTo(img_draw_screen);
			drawContours(img_draw_screen, max_flower_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0);
			
			//printf("mask_flag: %d, color_flag: %d\n",mask_flag,color_flag);
			if(contour_fail==CONTOUR_FAIL_BINARY_BACKGROUND)
				printf("CONTOUR_FAIL_BINARY_BACKGROUND\n");
			//char txt[20];
			//sprintf(txt,"flower contour %d",color_flag*20+mask_flag*50);
			//mycvShowImage(txt, img_draw_screen);
			if(contour_fail == SUCCESS){ //if we found the contour
				return;
			}
		}
	}

	throw ContourFailCenterPointException(); //failed to find contour
}

/*	equalizeLightness
*	input: an Image matrix
*	output: sets lightness in all pixels to be 50 
*/
void FlowerFeatureExtractor::equalizeLightness(Mat& img)
{
	Scalar rgb_pixel, hsl_pixel;
	Mat bgr[3];
	split(img, bgr);

	for(int i = 0; i < img.rows; i++) {
		for(int j = 0; j < img.cols; j++) {
			rgb_pixel = Scalar(bgr[2].at<uchar>(i,j),bgr[1].at<uchar>(i,j),bgr[0].at<uchar>(i,j));
			convertRGBtoHSL(rgb_pixel, hsl_pixel);

			if(hsl_pixel.val[2]>=LIGHTNESS_THRESHOLD_FIX) {
				hsl_pixel.val[2] = LIGHTNESS_VALUE_AFTER_FIX;

				convertHSLtoRGB(hsl_pixel,rgb_pixel);
				bgr[0].at<uchar>(i,j) = rgb_pixel.val[2];
				bgr[1].at<uchar>(i,j) = rgb_pixel.val[1];
				bgr[2].at<uchar>(i,j) = rgb_pixel.val[0];
			}
		}
	}
	merge(bgr, 3, img);
}

bool FlowerFeatureExtractor::isGrayscale(Scalar& color)
{
	Scalar hsl_color;
	convertRGBtoHSL(color, hsl_color);
	if(hsl_color.val[2]>=LIGHTNESS_THRESHOLD_WHITE || hsl_color.val[2]<=LIGHTNESS_THRESHOLD_BLACK || hsl_color.val[1]<=SATURATION_THRESHOLD_GRAY)
		return true;
	return false;
}

void FlowerFeatureExtractor::extractFeaturesFromImage(string image_path, Point& center) {
	cout << image_path << "\n";
	Scalar dom_inner_color, dom_flower_color, dom_background_color;	

	Mat img_opened;
	Mat img_resized;

	if(!PathFileExists(image_path.c_str())){
		throw PathNotFoundException();
	}

	img_opened = imread( image_path, CV_LOAD_IMAGE_COLOR );

	if(img_opened.empty()){
		throw EmptyImageException();
	}

	Utils::resizeImage(img_opened,img_resized);


	center.x=((center.x)*(img_resized.cols))/(img_opened.cols);
	center.y=((center.y)*(img_resized.rows))/(img_opened.rows);

	img_opened.release();

	Mat img_meanshift;
	//img_resized.copyTo(img_meanshift);

	Mat img_draw_screen;
	img_resized.copyTo(img_draw_screen);
	circle(img_draw_screen, center, 5, CV_RGB(0x00,0x00,0x00), -1, 8, 0);
	Utils::myShowImage("screen1", img_draw_screen);

	//equalizeLightness(img_meanshift);
	//mycvShowImage("screen1.3 lightness", img_meanshift);


	img_resized.convertTo(img_meanshift, -1, CONTRAST_RATIO);

	//cvSmooth(img_meanshift, img_meanshift,CV_MEDIAN,7,3);

	Utils::myShowImage("screen1.5 contrast", img_meanshift);
	applyMeanShift(img_meanshift); //changes img size!
	Utils::myShowImage("screen2", img_meanshift);

	
	img_meanshift.copyTo(img_draw_screen);


	/* OUTER PART */
	vector<Point> max_flower_contour;
	vector<vector<Point> > max_flower_contour_vec;

	getFlowerContourAndColor(max_flower_contour, img_meanshift, dom_inner_color, dom_flower_color, 
		dom_background_color, center);

	Scalar temp_hsl;
	convertRGBtoHSL(dom_flower_color, temp_hsl);
	//printf("DOM FL COL: R: %.2lf G: %.2lf B: %.2lf, H: %.2lf S: %.2lf L: %.2lf\n",dom_flower_color.val[0],dom_flower_color.val[1],dom_flower_color.val[2],temp_hsl.val[0],temp_hsl.val[1],temp_hsl.val[2]);

	//printf("Contour Length: %lf, Area: %lf, Ratio: %lf\n", cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED), cvContourArea(max_flower_contour),(cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED)*cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED))/ (cvContourArea(max_flower_contour)));

	/* UPDATE SAMPLE STRUCT */
	m_sample.m_dom_flower_color_red=dom_flower_color.val[0];
	m_sample.m_dom_flower_color_green=dom_flower_color.val[1];
	m_sample.m_dom_flower_color_blue=dom_flower_color.val[2];
	m_sample.m_dom_inner_color_red=dom_inner_color.val[0];
	m_sample.m_dom_inner_color_green=dom_inner_color.val[1];
	m_sample.m_dom_inner_color_blue=dom_inner_color.val[2];

	double flower_length = arcLength(max_flower_contour, true);
	double flower_area = contourArea(max_flower_contour);
	m_sample.m_length_area_ratio = flower_area*(4*PI)/(flower_length*flower_length);
	/* END UPDATE SAMPLE STRUCT */

	img_meanshift.copyTo(img_draw_screen);
	max_flower_contour_vec.push_back(max_flower_contour);
	drawContours(img_draw_screen, max_flower_contour_vec, 0, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens

	int num_points_max;
	int num_points_min;
	double radius_max;			//The radius of the flower.
	double min_max_flower_ratio;	//The ratio between the length of the closest min points(to center) to the the radius of the flower.
	double angle_max;			//The median angle between two min points.
	double angle_min;			//The median angle between two max points.


	getOuterFlowerProperties(num_points_max, num_points_min, angle_max, angle_min, 
		min_max_flower_ratio, radius_max, center, 
		max_flower_contour, img_meanshift, img_draw_screen);

	//printf("angle of min: %f\n", angle_min);

	/* UPDATE SAMPLE STRUCT */
	m_sample.m_num_points_max = num_points_max;
	m_sample.m_angle_max=angle_max;
	m_sample.m_min_max_flower_ratio=min_max_flower_ratio;
	m_sample.m_num_points_min = num_points_min;
	m_sample.m_angle_min = angle_min;
	/* END UPDATE SAMPLE STRUCT */

	/* END OUTER PART */

	/* INNER PART */

	double radius_inner;
	double inner_length;			//The length of the inner contour.

	getInnerPartProperties(radius_inner, inner_length, max_flower_contour, dom_flower_color, dom_inner_color,
		center, img_meanshift, img_draw_screen, radius_max);

	/* UPDATE SAMPLE STRUCT */
	m_sample.m_min_max_radius_ratio = radius_inner / radius_max;
	m_sample.m_length_inner_length_flower_ratio = inner_length / flower_length;
	/* END UPDATE SAMPLE STRUCT */

	/* END INNER PART */
}