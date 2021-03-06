#include "utils.h"

DIR* Utils::MyChdir(string path )
{
	if(!PathFileExists(path.c_str())){
		printf("ERROR: directory - %s does not exists\n",path);
		return NULL;
	}
	DIR* dir = opendir(path.c_str());
	if(_chdir(path.c_str())!=0){
		printf("ERROR: couldnt change working directory - %s\n",path);
		return NULL;
	}
	return dir;
}

void Utils::toLowercase(char* str){
	for(int i=0; i < strlen(str); i++){
		str[i]=tolower(str[i]);
	}
}


int Utils::getCenter(Point& center, FILE* centers,int flowerPicNo){
	int lineCounter=0;
	char line[20];
	int x,y,id;
	int flag = 1;

	rewind(centers);

	while(flag && fgets(line,20,centers)!=NULL){
		sscanf(line,"%d",&id);
		if(id==flowerPicNo)
			flag=0;
	}
	if(flag==1)
		return -1;
	sscanf(line,"%d,%d,%d",&id,&x,&y);
	center.x=x;
	center.y=y;
	return 0;	
}

void Utils::resizeImage(Mat& image, Mat& resized_image) {
	double ratio = (image.rows*1.0)/image.cols;
	int height, width;

	if(ratio < 1){
		height = Consts::BASE_SIZE;
		width = (int)(height*1.0/ratio);
	}
	else{
		width = Consts::BASE_SIZE;
		height = (int)width*ratio;
	}

	resize(image, resized_image, Size(width,height));
}

void Utils::myShowImage(const char* name, const Mat& image){
	if(Consts::DEBUG){
		imshow(name, image);
		cvWaitKey(10);
	}
}

void Utils::showAndWait(const char* name, Mat& image)
{
	imshow(name, image);
	cvWaitKey(0);
}

void Utils::convertRGBtoHSL (Scalar& rgb_color, Scalar& hsl_color){
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

double Utils::hue2rgb(double p, double q, double t){ //helper for convertHSLtoRGB
	if(t < 0) t += 1;
	if(t > 1) t -= 1;
	if(t < 1.0/6) return p + (q - p) * (6.0 * t);
	if(t < 1.0/2) return q;
	if(t < 2.0/3) return p + (q - p) * (2.0/3 - t) * 6.0;
	return p;
}

void Utils::convertHSLtoRGB(Scalar& hsl_color, Scalar& rgb_color)
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
		r = Utils::hue2rgb(p, q, h + 1.0/3);
		g = Utils::hue2rgb(p, q, h);
		b = Utils::hue2rgb(p, q, h - 1.0/3);
	}
	rgb_color.val[0] = r*255;
	rgb_color.val[1] = g*255;
	rgb_color.val[2] = b*255;

}

void Utils::bgr2rgb(Scalar& bgr, Scalar& rgb) 
{
	rgb.val[0] = bgr.val[2];
	rgb.val[1] = bgr.val[1];
	rgb.val[2] = bgr.val[0];
}

int Utils::square(int x)
{
	return x*x;
}

bool Utils::isGrayscale(Scalar& color)
{
	Scalar hsl_color;
	convertRGBtoHSL(color, hsl_color);
	if(hsl_color.val[2] >= Consts::LIGHTNESS_THRESHOLD_WHITE || hsl_color.val[2] <= Consts::LIGHTNESS_THRESHOLD_BLACK || hsl_color.val[1] <= Consts::SATURATION_THRESHOLD_GRAY)
		return true;
	return false;
}


void Utils::equalizeLightness(Mat& img)
{
	Scalar rgb_pixel, hsl_pixel;
	Mat bgr[3];
	split(img, bgr);

	for(int i = 0; i < img.rows; i++) {
		for(int j = 0; j < img.cols; j++) {
			rgb_pixel = Scalar(bgr[2].at<uchar>(i,j),bgr[1].at<uchar>(i,j),bgr[0].at<uchar>(i,j));
			convertRGBtoHSL(rgb_pixel, hsl_pixel);

			if(hsl_pixel.val[2] >= Consts::LIGHTNESS_THRESHOLD_FIX) {
				hsl_pixel.val[2] = Consts::LIGHTNESS_VALUE_AFTER_FIX;

				convertHSLtoRGB(hsl_pixel,rgb_pixel);
				bgr[0].at<uchar>(i,j) = rgb_pixel.val[2];
				bgr[1].at<uchar>(i,j) = rgb_pixel.val[1];
				bgr[2].at<uchar>(i,j) = rgb_pixel.val[0];
			}
		}
	}
	merge(bgr, 3, img);
}


double Utils::distance(int y1, int x1, int y2, int x2)
{
	double dst = square(x1-x2) + square(y1-y2); 
	return sqrt(dst);
}


double Utils::pointDistance(Point& point1, Point& point2)
{
	double dst = square(point1.x-point2.x) + square(point1.y-point2.y);
	return sqrt(dst);
}


int Utils::paddedMatrixAccess( Mat& image, int row, int col )
{
	if(row < 0 || col < 0 || col > image.cols-1 || row > image.rows-1)
		return 0;
	else
		return image.at<uchar>(row,col);

}

int Utils::colorDst(Scalar& color1, Scalar& color2)
{	
	double dst = square(color1.val[0]-color2.val[0]) + square(color1.val[1]-color2.val[1]) + square(color1.val[2]-color2.val[2]);
	return sqrt(dst);
}

void Utils::drawOneContour(Mat& draw_image, vector<Point>& contour,CvScalar& color)
{
	vector<vector<Point> > contour_vec;
	contour_vec.push_back(contour);
	drawContours(draw_image, contour_vec, 0, color, 2, 8, noArray(), 0); // Try different values of max_level, and see what happens
}


double Utils::angle(Point& center, Point& p1, Point& p2)
{
	Point v1,v2;
	v1.x = p1.x-center.x;
	v1.y = p1.y-center.y;

	v2.x = p2.x-center.x;
	v2.y = p2.y-center.y;

	return abs(acos(innerProduct(v1,v2)/sqrt(innerProduct(v1,v1)*innerProduct(v2,v2))));
}

double Utils::innerProduct(Point& v1, Point& v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}


double Utils::toDegree(double medAngle)
{
	return medAngle*180.0 / Consts::PI;
}


double Utils::contourVar(vector<Point>& contour, Point& center) {
	double Ex = 0; //E[x]
	double Ex2 = 0;//E[x^2]
	int i;

	for( i=0; i < contour.size(); ++i ){	
		Ex += pointDistance(contour.at(i), center);
		Ex2 += square(contour.at(i).x-center.x) + square(contour.at(i).y-center.y);
	}

	Ex /= i;
	Ex2 /= i;

	return Ex2 - Ex*Ex;
}


double Utils::avgDst(vector<Point>& contour, Point& center){
	double radius = 0;
	for(int i=0; i < contour.size(); ++i ){	
		radius += pointDistance(contour.at(i), center);
	}
	radius /= contour.size();
	return radius;
}


int Utils::getClosestPoint(int currentPointIndex, vector<Point>& min_max_points, vector<int>& isVisited)
{
	double dst, min = -1;
	int minIndex;
	Point minPoint;
	Point point = min_max_points.at(currentPointIndex);

	for(int i = 0; i <min_max_points.size(); i++)
	{		
		if(i == currentPointIndex || isVisited.at(i)==1)
			continue;
		dst = pointDistance(point, min_max_points.at(i));
		if(min == -1 || dst < min) {
			min = dst;
			minPoint = min_max_points.at(i);
			minIndex = i;
		}
	}
	return minIndex;
}


double Utils::median(vector<double>& vec)
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


double Utils::radiusOutOfPoints(vector<Point>& points, Point& center)
{
	double result;
	vector<double> distances;
	for(int i =0; i < points.size(); i++) {
		distances.push_back(pointDistance(points.at(i), center));
	}
	result = median(distances);
	return result;
}


double Utils::getMinMaxFlowerRatio(vector<Point> min_points, double outer_max_radius, Point& center)
{
	double min_distance;
	vector<double> distances;
	for(int i =0; i < min_points.size(); i++) {
		distances.push_back(pointDistance(min_points.at(i), center));
	}

	for(int i=0; i < min_points.size(); i++){
		if(i==0)
			min_distance = distances.at(i);
		else if(min_distance > distances.at(i))
			min_distance = distances.at(i);
	}

	return min_distance / outer_max_radius;
}


void Utils::contourToMatrix(vector<Point> contour, Mat& arr)
{
	for( int i=0; i<contour.size(); ++i ) {
		arr.at<uchar>(contour.at(i).y, contour.at(i).x) = 1;
	}
}

double Utils::contourToPointDst(vector<Point>& contour, Point& center)
{
	double minDst = -1;
	double dst =0;

	for(int i=0; i < contour.size(); ++i ){	
		dst = Utils::pointDistance(contour.at(i), center);
		if(minDst == -1 || dst < minDst)
			minDst = dst;
	}
	return minDst;
}

void Utils::intersectMasks(Mat& maskA, Mat& maskB, Mat& intersectMask)
{
	int x,y;
	for(x=0; x < maskA.cols; x++)
	{
		for(y=0; y < maskA.rows; y++)
		{


			if(maskA.at<uchar>(y,x) != 0 && maskB.at<uchar>(y,x) != 0) {
				intersectMask.at<uchar>(y,x) = 255; //White to greater of threshold
			} else{
				intersectMask.at<uchar>(y,x) = 0; //Black other
			}
		}
	}
}