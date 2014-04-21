#include "flower_detection_header.h"
//version 58

/*	localMinMaxPoints
*	input: an int array representing a binary image, its size, a point, an output array, resolution
*	output: outImage will have 2 at the cells where there are local maximum points of the distance function between the given point and the array point
*	not allows minimum points closer than <resolution> cells
*/
void FlowerFeatureExtractor::localMinMaxPointsInner(Mat& image, Mat& outImage, int resolution, 
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
				dst1 = Utils::distance(row , col, m_center.y, m_center.x);
				for(dx = -resolution ; dx < resolution + 1; dx++)
					for(dy = -resolution; dy < resolution+1 ;dy++)
					{
						if(Utils::padded_matrix_access(image, row+dx, col+dy) == 1 && (dx!=0 || dy!=0))
						{
							dst2 = Utils::distance(row+dx, col+dy, m_center.y, m_center.x);
							if (min_max_flag == MAX_POINTS_FLAG) {
								if(dst2 > dst1)
									flag = 1;
							} else if (min_max_flag == MIN_POINTS_FLAG) {
								if(dst2 < dst1)
									flag =1;
							}

							if(Utils::padded_matrix_access(outImage, row+dx, col+dy) == 2) //in case of more than one max point at the same area.
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
void FlowerFeatureExtractor::localMaxPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& max_points)
{	
	localMinMaxPointsInner(image, outImage, resolution, max_points, MAX_POINTS_FLAG);
}

/*	localMinPoints
*	input: an int array representing a binary image, its size, a point, an output array, resolution
*	output: outImage will have 2 at the cells where there are local minimum points of the distance function between the given point and the array point
*	not allows minimum points closer than <resolution> cells
*/
void FlowerFeatureExtractor::localMinPoints(Mat& image, Mat& outImage, int resolution, vector<Point>& min_points)
{	
	localMinMaxPointsInner(image, outImage, resolution, min_points, MIN_POINTS_FLAG);
}

/*	drawArrPointsOnImage
*	input: an image, an array to draw, colors and such
*	output: draws the extremum points from the array on the image, with specified color, thickness, etc
*/
void FlowerFeatureExtractor::drawPointsOnImage(Mat& contour_matrix_helper, const Scalar& color)
{
	Point now_point;
	for(int i=0; i<contour_matrix_helper.rows; i++){
		for(int j=0; j<contour_matrix_helper.cols; j++){
			if(contour_matrix_helper.at<uchar>(i,j) == 2){
				now_point.y=i;
				now_point.x=j;
				circle(m_draw_image, now_point, 2, color, -1, 8, 0);
			}
		}
	}
}

/* createCircleMask
*	input: an image (for mask output), radius, center
*	output: a binary image containing a filled circle, with specified radius, around center
*/
void FlowerFeatureExtractor::createCircleMask( Mat& img_mask, int radius)
{
	int x,y;
	for(x=0; x<img_mask.cols; x++){
		for(y=0 ;y<img_mask.rows; y++){
			if(Utils::square(x-m_center.x) + Utils::square(y-m_center.y) > Utils::square(radius)){
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
			Utils::bgr2rgb(bgr_color, rgb_color);

			point.x=x;
			point.y=y;

			if(Utils::colorDst(rgb_color, dom_color) < threshold && (contour.empty() || pointPolygonTest(contour, point, false) > 0 )) {
				img_output.at<uchar>(y,x) = 255; //White to greater of threshold
			} else{
				img_output.at<uchar>(y,x) = 0; //Black other
			}
		}
	}
	//showAndWait("moo",img_output);
}

/* findMaxContour
*	input: a pointer to a list of contours (first_contour), a pointer to a pointer of a contour
*	output: max_contour will contain the contour from the list with the largest amount of points on the contour
*/
void FlowerFeatureExtractor::findMaxContour(vector<vector<Point> >& contours, vector<Point>& max_contour)
{
	int maxsize = 0;
	int maxIndex = 0;

	max_contour = vector<Point>(); //return NULL if no contour
	for(int i = 0; i < contours.size(); i++) { //find maximal contour
		drawContours(m_draw_image, contours, i, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
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
void FlowerFeatureExtractor::findClosestToCenterPointContour(vector<vector<Point> >& contours, vector<Point>& inner_contour)
{
	double minDst = -1;
	double dst = 0;
	inner_contour = vector<Point>();

	for(int i = 0; i < contours.size(); i++) { //find maximal contour
		drawContours(m_draw_image, contours, i, CV_RGB(0x00,0xff,0x00), 2, 8, noArray(), 0); // Try different values of max_level, and see what happens

		//showAndWait(argv[0],img_draw_screen);

		if(contours.at(i).size() >= MINIMUM_INNER_CONTOUR_POINTS) {
			dst = Utils::contourToPointDst(contours.at(i), m_center);
			if(minDst == -1 || dst < minDst) {
				minDst = dst;
				inner_contour = contours.at(i);
			}
		}
	}
	//showAndWait("lol2",img_draw_screen);
}

/* applyMeanShift
*	input: an image
*	output: the image, after applying the meanshift algorithm on it. contains a fix to the size of the image, to allow the meanshift to work
*/
void FlowerFeatureExtractor::applyMeanShift()
{
	//Here comes the thing (fix for meanshift)
	m_image.cols &= -(1<<MEANSHIFT_LEVEL);
	m_image.rows &= -(1<<MEANSHIFT_LEVEL);
	pyrMeanShiftFiltering(m_image, m_image, MEANSHIFT_SPATIAL_RADIUS, MEANSHIFT_COLOR_RADIUS, MEANSHIFT_LEVEL);
}

/* createAndCalcHistogram
*	input: an image and a mask
*	output: a pointer to an allocated histogram of the image area inside img_mask
*/
CvHistogram* FlowerFeatureExtractor::createAndCalcHistogram(IplImage* img_hist, IplImage* img_mask){
	CvHistogram* hist;
	int numBins[3] = {NUM_BINS_FOR_COLOR, NUM_BINS_FOR_COLOR, NUM_BINS_FOR_COLOR};
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

/*	getOptimalInnerContour
*	input: max_flower_conotur (contour of outer flower), pointer to max_inner_contour pointer (to update pointer), dominant colors, 
circle_around_center_ponit_flag - this flag will be activated (by the calling function) if no inner contour found the first time this function is run. if this flag is on, the inner contour will be a circle around the center point with radius of the average distance between the center point and the contour with the smallest variance.
radius_max in order to elimnate bad contours - in case of the radius of the circle around center point is too big (more than half of the radius_max).
*	output: runs on color thresholds (according to OPTIMAL_INNER_CONTOUR_START, OPTIMAL_INNER_CONTOUR_END, OPTIMAL_INNER_CONTOUR_QUOTIENT,
*			paints the are which is within <threshold> distance from the center color, and within the outer flower contour, and tries to find contour there.
picks the contour with smallest variance of distances between a point on the contour and the center point
*	fail: max_inner_contour=NULL if no inner contour found.
*/
void FlowerFeatureExtractor::getOptimalInnerContour(vector<Point>& max_inner_contour, Scalar& dom_flower_color,
							Scalar& dom_inner_color, int circle_around_center_ponit_flag, double radius_max )
{

	Mat img_binary_small = Mat(m_image.size(), CV_8UC1);
	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	vector<Point> optimal_contour, current_inner_contour;
	double min_variance, current_variance;
	Point2f center2f = Point2f(m_center.x, m_center.y);
	int colorDistance;
	int threshold_small;
	char name[100];

	for(int i = OPTIMAL_INNER_CONTOUR_START; i < OPTIMAL_INNER_CONTOUR_END; i++){
		colorDistance = Utils::colorDst(dom_flower_color, dom_inner_color);
		threshold_small = colorDistance * ((i*1.0) / OPTIMAL_INNER_CONTOUR_QUOTIENT);
		threshold_small = threshold_small < INNER_CONTUOR_MINIMUM_THRESHOLD ? INNER_CONTUOR_MINIMUM_THRESHOLD : threshold_small;
		//printf("SMALL THRESHOLD - %d",threshold_small);

		toBinaryByDominantColorWithContour(m_image, img_binary_small, dom_inner_color, threshold_small, m_max_flower_contour);
		sprintf(name,"smallbw %d",i);
		Utils::myShowImage(name, img_binary_small);

		findContours(img_binary_small, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
		if(circle_around_center_ponit_flag == 0)
			findMaxContour(contours, current_inner_contour);
		else
			findClosestToCenterPointContour(contours, current_inner_contour);

		if(current_inner_contour.empty()){	
			if(colorDistance <= INNER_CONTOUR_EQUAL_FLOWER_CONTOUR_THRESHOLD){	//if there is no inner contour the inner and outer contours are identical.
				max_inner_contour = m_max_flower_contour;
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

		current_variance = Utils::contourVar(current_inner_contour, m_center);
		if(optimal_contour.empty() || min_variance > current_variance) {
			min_variance = current_variance;
			optimal_contour = current_inner_contour;
		}
		m_image.copyTo(m_draw_image);
		Utils::drawOneContour(m_draw_image, current_inner_contour, CV_RGB(0x00,0xff,0x00));
		//showAndWait("Opencv_test2.exe", img_draw_screen );
	}

	if(circle_around_center_ponit_flag == 1 && !optimal_contour.empty()) {

		Mat circle_image = Mat(m_image.size(), CV_8UC1);
		double radius = Utils::avgDst(optimal_contour, m_center);

		if(radius <= INNER_RADIUS_TO_FLOWER_RADIUS_RATIO * radius_max) {
			createCircleMask(circle_image, radius);
			findContours(circle_image, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
			findMaxContour(contours, optimal_contour);
		}
		else {
			optimal_contour = vector<Point>();
		}
	}

	max_inner_contour = optimal_contour;
}

/*
*	input:	- array of point (min or max points) in size greater than 2.
*			- center point
*	output:	The median of all angles between two close points and the center.
*/
double FlowerFeatureExtractor::getMedianAngle(vector<Point>& min_max_points, double& outMinAngle, 
											  int* outMinAngleTwoPointsIndexes)
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
		closestPointIndex = Utils::getClosestPoint(currentPointIndex, min_max_points, isVisited);
		isVisited.at(currentPointIndex) = 1;
		angles.at(i) = Utils::angle(m_center, min_max_points.at(currentPointIndex), min_max_points.at(closestPointIndex));

		if(outMinAngleTwoPointsIndexes != NULL){
			if(angles.at(i) < minAngle){
				minAngle = angles.at(i);
				outMinAngle = minAngle;
				outMinAngleTwoPointsIndexes[0] = currentPointIndex;
				outMinAngleTwoPointsIndexes[1] = closestPointIndex;
			}
		}

		circle(m_draw_image, min_max_points.at(currentPointIndex), 2, CV_RGB(0xff,0x00,0x00), -1, 8, 0);
		//showAndWait("test",draw_img);
		//printf("%lf \n",angles[i]*180/3.14159);*/
		currentPointIndex = closestPointIndex;
	}

	if(i == min_max_points.size() - 1)
	{
		angles.at(i) = Utils::angle(m_center, min_max_points.at(currentPointIndex), min_max_points.at(0));

		if(outMinAngleTwoPointsIndexes != NULL){
			if(angles.at(i) < minAngle){
				minAngle = angles.at(i);
				outMinAngle = minAngle;				
				outMinAngleTwoPointsIndexes[0] = currentPointIndex;
				outMinAngleTwoPointsIndexes[1] = 0;
			}
		}

		//printf("%lf \n",angles[i]*180/3.14159);
		circle(m_draw_image, min_max_points.at(currentPointIndex), 2, CV_RGB(0xff,0x00,0x00), -1, 8, 0);
		//showAndWait("test",draw_img);
	}

	result = Utils::median(angles);

	return result;
}

/*	getFixBadMinMaxPointsThreshold
*	input: median of angles between min or max points
*	output: threshold, according to medAngle (big, medium, or small threshold)
*/
int FlowerFeatureExtractor::getFixBadMinMaxPointsThreshold(double medAngle)
{
	double angle = Utils::toDegree(medAngle);
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
int FlowerFeatureExtractor::fixBadMinMaxPoints(vector<Point>& min_max_points, double& outMedAngle, Mat& img_with_contours)
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

	medAngle = getMedianAngle(min_max_points, minAngle, minPointsIndex);
	if(medAngle <= 0) {
		throw MedianAngleException();
	}

	fix_bad_min_max_points_threshold = getFixBadMinMaxPointsThreshold(medAngle);

	while(minAngle < medAngle/fix_bad_min_max_points_threshold)
	{		
		for(i=0; i < min_max_points.size(); i++)
			isVisited.at(i)=0;

		isVisited.at(minPointsIndex[1]) = 1;

		secondClosestPoinIndex = Utils::getClosestPoint(minPointsIndex[0], min_max_points, isVisited);
		dst1 = Utils::point_distance(min_max_points.at(minPointsIndex[0]), min_max_points.at(secondClosestPoinIndex));

		// printf("\n dst1 == %lf\n",dst1);
		// printf("\n point1 - x=%d y=%d\n",(*pointsArr)[minPointsIndex[0]].x,(*pointsArr)[minPointsIndex[0]].y);

		isVisited.at(minPointsIndex[1]) = 0;
		isVisited.at(minPointsIndex[0]) = 1;
		secondClosestPoinIndex = Utils::getClosestPoint(minPointsIndex[1],min_max_points, isVisited);
		dst2 = Utils::point_distance(min_max_points.at(minPointsIndex[1]), min_max_points.at(secondClosestPoinIndex));

		//printf("\n dst2 == %lf\n",dst2);
		// printf("\n point2 - x=%d y=%d\n",(*pointsArr)[minPointsIndex[1]].x,(*pointsArr)[minPointsIndex[1]].y);

		removePointIndex = dst1<dst2?	minPointsIndex[0] :	minPointsIndex[1];
		min_max_points.erase(min_max_points.begin() + removePointIndex);

		img_with_contours.copyTo(m_draw_image);

		medAngle = getMedianAngle(min_max_points, minAngle, minPointsIndex);
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
void FlowerFeatureExtractor::createCircleMasks( Mat& img_mask_small, Mat& img_mask_large, Mat& img_mask_background, int mask_flag, int background_mask_flag)
{
	int radius_large = m_center.x > m_center.y ? m_center.y : m_center.x;	
	int radius_background = radius_large;

	radius_large*=1 - (COLOR_FLAG_STEP * mask_flag);
	radius_background*= 1 + (COLOR_FLAG_STEP * background_mask_flag);

	createCircleMask(img_mask_small, RADIUS_SMALL);
	createCircleMask(img_mask_large, radius_large);
	createCircleMask(img_mask_background, radius_background);	
	flipMask(img_mask_background);
}

/*	getDominantColors
*	input: dominant color struct pointers, mask_flag (changes inner mask radius), background_mask_flag (changes background mask radius)
*	output: dominant colors, according to the flags given
*/
void FlowerFeatureExtractor::getDominantColors(Scalar& dom_inner_color, Scalar& dom_flower_color, Scalar& dom_background_color, 
					   Mat& img, int mask_flag, int background_mask_flag)
{

	/* CREATE MASK */
	Mat img_mask_small = Mat(img.size(), CV_8UC1);//inner circle mask
	Mat img_mask_large = Mat(img.size(), CV_8UC1); //whole flower mask 
	Mat img_mask_background = Mat(img.size(), CV_8UC1); //background mask

	createCircleMasks(img_mask_small, img_mask_large, img_mask_background, mask_flag, background_mask_flag);

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
	return Utils::colorDst(dom_flower_color,dom_background_color)*(OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR + color_flag)/OUTER_COLOR_DST_THRESHOLD_DENOMINATOR;
}

/* checkBinaryBackground
*	input: binary image, background_mask_flag (changes background mask radius)
*	output: "masks" the image with background mask, and checks if there are more black pixels than white pixels in the image
*	fail: CONTOUR_FAIL_BINARY_BACKGROUND if there are more white pixels than black pixels in that area (what would indicate a bad binarization - too much background)
*/
void FlowerFeatureExtractor::checkBinaryBackground(Mat& img_binary, int background_mask_flag)
{
	Mat img_mask_background = Mat(img_binary.size(), CV_8UC1);
	Mat img_binary_after_mask = Mat(img_binary.size(), CV_8UC1);

	int radius_background = m_center.x > m_center.y ? m_center.y : m_center.x;	

	radius_background*= 1 + (COLOR_FLAG_STEP * background_mask_flag);

	createCircleMask(img_mask_background, radius_background);
	flipMask(img_mask_background);

	img_binary.copyTo(img_binary_after_mask, img_mask_background);

	int x, y, black_count = 0, white_count = 0;
	//showAndWait("bin mask",img_binary_after_mask);
	for(x=0; x < img_binary_after_mask.cols; x++){
		for(y=0; y < img_binary_after_mask.rows; y++){
			if(Utils::square(x-m_center.x) + Utils::square(y-m_center.y) <= Utils::square(radius_background))
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
int FlowerFeatureExtractor::getOptimalFlowerContour(Mat& img, Scalar& dom_flower_color, Scalar& dom_background_color, 
													int color_flag, int background_mask_flag)
{
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	Mat img_binary = Mat(img.size(), CV_8UC1);

	int threshold_large = getThresholdByFlag(color_flag, dom_flower_color, dom_background_color);

	//printf("using threshold_large: %d\n",threshold_large);
	try 
	{
		toBinaryByDominantColorWithContour(img, img_binary, dom_flower_color, threshold_large, vector<Point>());
		checkBinaryBackground(img_binary, background_mask_flag);

		Utils::myShowImage("screen5", img_binary);

		findContours(img_binary, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);
		if (contours.size() == 0) {
			return CONTOUR_FAIL_NO_POINTS;
		}

		img.copyTo(m_draw_image);

		findMaxContour(contours, m_max_flower_contour);

		Point2f center2f;
		center2f.x = m_center.x;
		center2f.y = m_center.y;
		if(pointPolygonTest(m_max_flower_contour, center2f, 0) <= 0)
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
void FlowerFeatureExtractor::getOuterFlowerProperties(int& num_points_max, int& num_points_min, double& angle_max, 
													  double& angle_min, double& min_max_flower_ratio, double& radius_max)
{
	vector<Point> max_points;		//Will hold the flower max points.
	vector<Point> min_points;		//Will hold the flower min points.

	Mat contour_matrix = Mat::zeros(m_image.size(), CV_8UC1);
	Mat contour_matrix_helper = Mat::zeros(m_image.size(), CV_8UC1);

	Utils::contourToMatrix(m_max_flower_contour, contour_matrix);
	contour_matrix.copyTo(contour_matrix_helper);

	localMaxPoints(contour_matrix, contour_matrix_helper, MAX_POINTS_THRESHOLD, max_points);
	drawPointsOnImage(contour_matrix_helper, CV_RGB(0x00,0x00,0x00));
	if(max_points.size() == 0)
		throw NoOuterMaxPointsFoundException();

	radius_max = Utils::radiusOutOfPoints(max_points, m_center);

	Mat img_with_contours = Mat(m_image.size(), m_image.type());
	m_draw_image.copyTo(img_with_contours);
	
	num_points_max = fixBadMinMaxPoints(max_points, angle_max, img_with_contours);
	contour_matrix.copyTo(contour_matrix_helper);
	localMinPoints(contour_matrix, contour_matrix_helper, MIN_POINTS_THRESHOLD, min_points);
	drawPointsOnImage(contour_matrix_helper, CV_RGB(0x00,0x00,0xff));
	Utils::myShowImage("screen6", m_draw_image);

	if(min_points.size() == 0)
		throw NoOuterMinPointsFoundException(); // NO_OUTER_MIN_POINTS_FOUND_ERROR;


	min_max_flower_ratio = Utils::getMinMaxFlowerRatio(min_points, radius_max, m_center);

	m_image.copyTo(m_draw_image);

	Utils::drawOneContour(m_draw_image, m_max_flower_contour, CV_RGB(0x00,0xff,0x00));
	drawPointsOnImage(contour_matrix_helper, CV_RGB(0x00,0x00,0x00));

	m_draw_image.copyTo(img_with_contours);

	num_points_min = fixBadMinMaxPoints(min_points, angle_min, img_with_contours); //updates angle_min
}

/*	getInnerPartProperties
*	input: outer contour of the flower(in order to create a mask, to define where we look for the inner contour), radius_inner to update, colors and so on (radius_max in order to elimnate bad contours)
*	output: updates sample properties related to the inner part of the flower. for now: radius_inner
*/
void FlowerFeatureExtractor::getInnerPartProperties(double& radius_inner, double& inner_length, 
													Scalar& dom_flower_color, Scalar& dom_inner_color, double radius_max)
{
	vector<Point> max_inner_contour;

	int num_points_max_inner, num_points_min_inner;

	getOptimalInnerContour(max_inner_contour, dom_flower_color, dom_inner_color, 0, radius_max);
	if(max_inner_contour.empty())
		getOptimalInnerContour(max_inner_contour, dom_flower_color,dom_inner_color, 1, radius_max);
	if(max_inner_contour.empty())
		throw MaxInnerContourException();// NO_INNER_CONTOUR_FOUND_ERROR;

	Mat out_arr_max = Mat(m_image.size(), CV_8UC1);
	Mat out_arr_min = Mat(m_image.size(), CV_8UC1);
	Mat contour_arr = Mat(m_image.size(), CV_8UC1);
	vector<Point> max_points_inner, min_points_inner, combined_points;

	/* END ALLOC ARRS */
	inner_length = arcLength(max_inner_contour, true);
	//printf("Inner contour Length: %lf, Area: %lf, Ratio: %lf\n", cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED), cvContourArea(max_inner_contour),(cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED)*cvArcLength(max_inner_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED))/ (cvContourArea(max_inner_contour)));
	m_image.copyTo(m_draw_image);
	Utils::drawOneContour(m_draw_image, max_inner_contour, CV_RGB(0x00,0xff,0x00));

	Utils::contourToMatrix(max_inner_contour, contour_arr);

	contour_arr.copyTo(out_arr_max);
	contour_arr.copyTo(out_arr_min);

	localMaxPoints(contour_arr, out_arr_max, MAX_POINTS_THRESHOLD, max_points_inner);
	drawPointsOnImage(out_arr_max, CV_RGB(0xff,0x00,0x00));

	localMinPoints(contour_arr, out_arr_min, MIN_POINTS_THRESHOLD, min_points_inner);
	drawPointsOnImage(out_arr_min, CV_RGB(0x00,0x00,0xff));
	Utils::myShowImage("screen7", m_draw_image);

	if(max_points_inner.empty() || min_points_inner.empty())
		throw NoInnerMinMaxPointsException(); //NO_INNER_MIN_MAX_POINTS_FOUND_ERROR;

	for(int i = 0; i < max_points_inner.size(); i ++) {
		combined_points.push_back(max_points_inner.at(i));
	}

	for(int i = 0; i < min_points_inner.size(); i ++) {
		combined_points.push_back(min_points_inner.at(i));
	}

	radius_inner = Utils::radiusOutOfPoints(combined_points, m_center);

}

/*	getFlowerContourAndColor
*	input: pointer to pointer of max_flower_contour (in order to update it), color struct pointers
*	output: updates outer flower contour and dominant colors of the flower, inner, and background.
*	may try several background and inner mask radiuses, and several color thresholds, if it gets bad contours (no points, or center point not inside contour)
*/
void FlowerFeatureExtractor::getFlowerContourAndColor(Scalar& dom_inner_color, Scalar& dom_flower_color, 
													  Scalar& dom_background_color)
{
	int mask_flag, color_flag, contour_fail, background_mask_flag;
	int mask_end = OUTER_COLOR_DST_THRESHOLD_MAX_NUMERATOR - OUTER_COLOR_DST_THRESHOLD_START_NUMERATOR;
	Scalar dom_color_temp_hsl;
	/*LIGHTNESS FIX*/
	Scalar dom_inner_color_outer_contour, dom_flower_color_outer_contour, dom_background_color_outer_contour;
	Mat img_outer_contour;

	m_image.copyTo(img_outer_contour);
	Utils::equalizeLightness(img_outer_contour);
	img_outer_contour.convertTo(img_outer_contour, -1, CONTRAST_RATIO_OUTER_CONTOUR); //contrast
	Utils::myShowImage("lightness",img_outer_contour);

	for(background_mask_flag=0; background_mask_flag < BACKGROUND_MASK_FLAG_ITERATIONS; background_mask_flag++) { //if no contour at all is found, probably background color and flower color are the same. then, try to increase background mask radius to fix this
		getDominantColors(dom_inner_color, dom_flower_color, dom_background_color, m_image, 0, background_mask_flag);
		Utils::convertRGBtoHSL(dom_flower_color, dom_color_temp_hsl); //for printing color
		
		if(Utils::isGrayscale(dom_flower_color)) { //don't fix lightness
			contour_fail = getOptimalFlowerContour(m_image, dom_flower_color, dom_background_color, 0, 0);
		}

		else { //lightness fix

			getDominantColors(dom_inner_color_outer_contour, dom_flower_color_outer_contour, 
				dom_background_color_outer_contour, img_outer_contour, 0, background_mask_flag);

			contour_fail = getOptimalFlowerContour(img_outer_contour, dom_flower_color_outer_contour, 
												   dom_background_color_outer_contour, 0, 0);
			//printf("DOM COLOR1: H:%.2lf S:%.2lf L:%.2lf ",dom_color_temp_hsl.val[0],dom_color_temp_hsl.val[1],dom_color_temp_hsl.val[2]);
			//printf("first, lightness fix APPLIED\n");
		}

		m_image.copyTo(m_draw_image);
		Utils::drawOneContour(m_draw_image, m_max_flower_contour, CV_RGB(0x00,0xff,0x00));
		
		if(contour_fail != CONTOUR_FAIL_NO_POINTS)	//Here we care only about that error, because we want to make sure that we got the right background dominant color. Therefore we dont check CONTOUR_FAIL_BINARY_BACKGROUND
			break;
	}
	for(color_flag=0; color_flag < COLOR_FLAG_ITERATIONS; color_flag++) { //runs on color thresholds
		for(mask_flag=0; mask_flag < mask_end; mask_flag++) { //runs on mask radius	
			getDominantColors(dom_inner_color, dom_flower_color, dom_background_color, 
				m_image, mask_flag, background_mask_flag); //mask_flag to set radius for mask

			Utils::convertRGBtoHSL(dom_flower_color, dom_color_temp_hsl);

			if(Utils::isGrayscale(dom_flower_color)) { //don't fix lightness
				contour_fail = getOptimalFlowerContour(m_image, dom_flower_color, dom_background_color, 
													   color_flag, background_mask_flag);
			}
			else {//lightness fix
				getDominantColors(dom_inner_color_outer_contour, dom_flower_color_outer_contour, 
					dom_background_color_outer_contour, img_outer_contour, 0, background_mask_flag);

				contour_fail = getOptimalFlowerContour(img_outer_contour, dom_flower_color_outer_contour, 
													   dom_background_color_outer_contour, color_flag, 
													   background_mask_flag);
				//printf("DOM COLOR2: H:%.2lf S:%.2lf L:%.2lf ",dom_color_temp_hsl.val[0],dom_color_temp_hsl.val[1],dom_color_temp_hsl.val[2]);
				//printf("second, lightness fix APPLIED\n");
			}

			m_image.copyTo(m_draw_image);
			Utils::drawOneContour(m_draw_image, m_max_flower_contour, CV_RGB(0x00,0xff,0x00));
			
			//printf("mask_flag: %d, color_flag: %d\n",mask_flag,color_flag);
			if(contour_fail==CONTOUR_FAIL_BINARY_BACKGROUND)
				printf("CONTOUR_FAIL_BINARY_BACKGROUND\n");

			if(contour_fail == SUCCESS){ //if we found the contour
				return;
			}
		}
	}

	throw ContourFailCenterPointException(); //failed to find contour
}


void FlowerFeatureExtractor::extractFeaturesFromImage() {	
	Scalar dom_inner_color, dom_flower_color, dom_background_color;	

	Mat img_resized;

	Utils::resizeImage(m_image, img_resized);


	m_center.x=((m_center.x)*(img_resized.cols))/(m_image.cols);
	m_center.y=((m_center.y)*(img_resized.rows))/(m_image.rows);

	m_image = img_resized;
	m_image.copyTo(m_draw_image);
	circle(m_draw_image, m_center, 5, CV_RGB(0x00,0x00,0x00), -1, 8, 0);
	Utils::myShowImage("screen1", m_draw_image);

	//equalizeLightness(img_meanshift);
	//mycvShowImage("screen1.3 lightness", img_meanshift);

	m_image.convertTo(m_image, -1, CONTRAST_RATIO);

	//cvSmooth(img_meanshift, img_meanshift,CV_MEDIAN,7,3);

	Utils::myShowImage("screen1.5 contrast", m_image);
	applyMeanShift(); //changes img size!
	Utils::myShowImage("screen2", m_image);

	
	m_image.copyTo(m_draw_image);


	/* OUTER PART */	
	

	getFlowerContourAndColor(dom_inner_color, dom_flower_color, dom_background_color);

	Scalar temp_hsl;
	Utils::convertRGBtoHSL(dom_flower_color, temp_hsl);
	//printf("DOM FL COL: R: %.2lf G: %.2lf B: %.2lf, H: %.2lf S: %.2lf L: %.2lf\n",dom_flower_color.val[0],dom_flower_color.val[1],dom_flower_color.val[2],temp_hsl.val[0],temp_hsl.val[1],temp_hsl.val[2]);

	//printf("Contour Length: %lf, Area: %lf, Ratio: %lf\n", cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED), cvContourArea(max_flower_contour),(cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED)*cvArcLength(max_flower_contour,CV_WHOLE_SEQ,CV_SEQ_FLAG_CLOSED))/ (cvContourArea(max_flower_contour)));

	/* UPDATE SAMPLE STRUCT */
	m_sample.m_dom_flower_color_red=dom_flower_color.val[0];
	m_sample.m_dom_flower_color_green=dom_flower_color.val[1];
	m_sample.m_dom_flower_color_blue=dom_flower_color.val[2];
	m_sample.m_dom_inner_color_red=dom_inner_color.val[0];
	m_sample.m_dom_inner_color_green=dom_inner_color.val[1];
	m_sample.m_dom_inner_color_blue=dom_inner_color.val[2];

	double flower_length = arcLength(m_max_flower_contour, true);
	double flower_area = contourArea(m_max_flower_contour);
	m_sample.m_length_area_ratio = flower_area*(4*PI)/(flower_length*flower_length);
	/* END UPDATE SAMPLE STRUCT */

	m_image.copyTo(m_draw_image);	
	Utils::drawOneContour(m_draw_image, m_max_flower_contour, CV_RGB(0x00,0xff,0x00));

	int num_points_max;
	int num_points_min;
	double radius_max;			//The radius of the flower.
	double min_max_flower_ratio;	//The ratio between the length of the closest min points(to center) to the the radius of the flower.
	double angle_max;			//The median angle between two min points.
	double angle_min;			//The median angle between two max points.


	getOuterFlowerProperties(num_points_max, num_points_min, angle_max, angle_min, min_max_flower_ratio, radius_max);

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

	getInnerPartProperties(radius_inner, inner_length, dom_flower_color, dom_inner_color, radius_max);

	/* UPDATE SAMPLE STRUCT */
	m_sample.m_min_max_radius_ratio = radius_inner / radius_max;
	m_sample.m_length_inner_length_flower_ratio = inner_length / flower_length;
	/* END UPDATE SAMPLE STRUCT */

	/* END INNER PART */
}

FlowerFeatureExtractor::FlowerFeatureExtractor(string image_path, Point& center) 
{
	cout << image_path << "\n";
	if(!PathFileExists(image_path.c_str())){
		throw PathNotFoundException();
	}

	m_image = imread(image_path);

	if(m_image.empty()){
		throw EmptyImageException();
	}

	m_center = Point(center);
}