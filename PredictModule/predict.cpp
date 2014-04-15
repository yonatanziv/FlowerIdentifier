#include "flower_detection_header.h"

#define SVM_SERIALIZATION_DIR_PATH "samples\\serialization"



/* command line input: argv[1]=predict_img_path, argv[2]=center.x, argv[3]= center.y
 * This program uses the samples\serialization\svmSerialize.txt file created by the train module, 
 * takes the file given as a parameter and predicts which flower it is.
 * return -1 on error, and the id of the input flower otherwise
*/
int main(int argc, char* argv[]) {

	Svm_Sample predict_sample;
	CvSVM SVM;
	double res =-1;
	int ret;
	Point center;
	int flowerID=0;
	int flowerPicNo=0;

	if(argc != 4){
		printf("ERROR: Bad Argument Number - %d.",argc);
		cvWaitKey(0);
		return -1;
	}
			
	printf("\n\n Predict:\n");
	
	if(_chdir(SVM_SERIALIZATION_DIR_PATH)!=0){
		printf("ERROR: couldnt change working directory - %s",SVM_SERIALIZATION_DIR_PATH);
		return -1;
	}
	if(!exists(SVM_SERIALIZE_FILE_NAME)){
		printf("ERROR: file is not exists - %s",SVM_SERIALIZE_FILE_NAME);
		return -1;
	}
	SVM.load(SVM_SERIALIZE_FILE_NAME);
	_chdir("..");
	_chdir("..");

	center.x=atoi(argv[2]);
	center.y=atoi(argv[3]);
	ret=getSampleDataFromImage(argv[1], &predict_sample,center);
	if(ret!=-1){
		res=svmPredict(&SVM,&predict_sample);
		printf("** Predicting: %s is %lf **\n\n", argv[1], res);
	}
	else
		printf("Image %s dropped!\n", argv[1]);	
	waitKey(0);
	return (int)res;
}

