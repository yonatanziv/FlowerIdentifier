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


