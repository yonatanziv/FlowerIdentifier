#include "utils.h"

DIR* Utils::MyChdir(char * path){
	if(!PathFileExists(path)){
		printf("ERROR: directory - %s does not exists\n",path);
		return NULL;
	}
	DIR* dir = opendir(path);
	if(_chdir(path)!=0){
		printf("ERROR: couldnt change working directory - %s\n",path);
		return NULL;
	}
	return dir;
}

void Utils::toLowercase(char* str){
	for(int i=0; i<strlen(str); i++){
		str[i]=tolower(str[i]);
	}
}

/*	getCenter
*	input: A center point pointer(to be update), a center points file, a number of a flower picture
*	output: Updates center to contain the center of the flower, as written in the file. return 0 on success and -1 on fail.
*	file format: each line contains:
*				flowerPicNo,x,y
*/
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
		height = BASE_SIZE;
		width = (int)(height*1.0/ratio);
	}
	else{
		width = BASE_SIZE;
		height = (int)width*ratio;
	}

	resize(image, resized_image, Size(width,height));
}

void Utils::myShowImage(const char* name, const Mat& image){
	if(DEBUG){
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
	double ret;
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
	if(hsl_color.val[2]>=LIGHTNESS_THRESHOLD_WHITE || hsl_color.val[2]<=LIGHTNESS_THRESHOLD_BLACK || hsl_color.val[1]<=SATURATION_THRESHOLD_GRAY)
		return true;
	return false;
}

/*	equalizeLightness
*	input: an Image matrix
*	output: sets lightness in all pixels to be 50 
*/
void Utils::equalizeLightness(Mat& img)
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

/*	distance
*	input: two points (x1,y1),(x2,y2)
*	output: euclidean distance between them
*/
double Utils::distance(int y1, int x1, int y2, int x2)
{
	double dst = square(x1-x2) + square(y1-y2); 
	return sqrt(dst);
}

/*	distance
*	input: two CvPoint points
*	output: euclidean distance between them
*/
double Utils::point_distance(Point& point1, Point& point2)
{
	double dst = square(point1.x-point2.x) + square(point1.y-point2.y);
	return sqrt(dst);
}

/*	padded_array_access
*	input: an int array, its size, and a cell in it (row and col)
*	output: returns image[row][col] if the indices are inside the array, 0 otherwise
*	used to prevent falling out of array in localMaxPoints and localMinPoints
*/
int Utils::padded_matrix_access( Mat& image, int row, int col )
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