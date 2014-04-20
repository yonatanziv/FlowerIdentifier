#include "flower_detection_header.h"

#define SVM_SERIALIZATION_DIR_PATH "samples\\serialization"

/* command line input: argv[1]=predict_img_path, argv[2]=center.x, argv[3]= center.y
 * This program uses the samples\serialization\svmSerialize.txt file created by the train module, 
 * takes the file given as a parameter and predicts which flower it is.
 * return -1 on error, and the id of the input flower otherwise
*/
int main(int argc, char* argv[]) {

	FlowerFeatureExtractor flower_featurs;
	SvmModel svm;
	double res =-1;
	Point center;

	if(argc != 4){
		cout << "Usage: <exe> <flower_image_path> <center_x> <center_y>";
		return -1;
	}

	try 
	{
		cout << "\n\n Predict:\n";

		svm.load(SVM_SERIALIZE_PATH);

		center.x=atoi(argv[2]);
		center.y=atoi(argv[3]);

		flower_featurs.extractFeaturesFromImage(argv[1], center);
		res = svm.predict(flower_featurs.m_sample);
		printf("** Predicting: %s is %lf **\n\n", argv[1], res);

	}  catch(exception& e) {
		cout << e.what() << "\n";
	}

	waitKey(0);
	return (int)res;
} 

