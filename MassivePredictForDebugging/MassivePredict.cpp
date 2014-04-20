#include "utils.h"

#define SVM_SERIALIZATION_DIR_PATH "samples\\serialization"
#define PREDICT_DIR_PATH "predict"
#define CENTERS_FILE_NAME "centers.txt"


/* command line input: argv[1]=predict_img_path, argv[2]=center.x, argv[3]= center.y
   return -1 on error, 0 if couldnt recognize the image, and the id of the input flower otherwise*/
int main(int argc, char* argv[]) {
	FlowerFeatureExtractor flower_featurs;
	SvmModel svm;
	double res =-1;
	Point center;
	DataSample current_sample;
	TrainSet trainSet;
	char temp_filename[257];

	FILE* centers=NULL;
	int center_res;

	try 
	{
		svm.load(SVM_SERIALIZE_PATH);

		printf("\n\n Predict:\n");
		//loop over files

		/******** DIR SCANNING CODE *************/

		int flowerID=0;
		int flowerPicNo=0;

		DIR *dir,*subdir, *temp;
		struct dirent *ent;	

		/* print contents of directories listed in command line */

		/* open directory stream */
		dir = Utils::MyChdir(PREDICT_DIR_PATH);
		if (dir != NULL) {

			centers = fopen(CENTERS_FILE_NAME,"r");

			/* print all the files and directories within directory */
			while ((ent = readdir (dir)) != NULL) {
				switch (ent->d_type) {
				case DT_REG:


					/******************** YOUR CODE HERE *************************/

					sscanf(ent->d_name, "%d.jpg", &flowerPicNo);
					//printf ("ID:%d, No:%d:\n",flowerID, flowerPicNo);
					Utils::toLowercase(ent->d_name + strlen(ent->d_name)-4);
					if(strcmp(".jpg",ent->d_name + strlen(ent->d_name)-4)!=0){
						break;
					}
					strcpy(temp_filename,PREDICT_DIR_PATH);
					strcat(temp_filename,"\\");
					strcat(temp_filename,ent->d_name);
					Utils::getCenter(center,centers,flowerPicNo);
					if(center.x!=-1 && center.y!=-1)
					{
						try {
							flower_featurs.extractFeaturesFromImage(ent->d_name, center);
							res = svm.predict(flower_featurs.m_sample);
							printf("** Predicting: %s is %lf **\n\n", ent->d_name, res);	
						}   catch(exception& e) {
							cout << e.what() << "\n";
						}
										
					}else
						printf("Image %s dropped!\n", temp_filename);
					

					/**************************************************************/
					break;

				case DT_DIR:
					//printf ("%s (dir)\n", ent->d_name);
					break;

					/*default:
					break;*/
				}
			}

			closedir (dir);
			fclose(centers);
		} else {
			/* could not open directory */
			perror ("");
		}

		/******** END DIR SCANNING CODE FOR PREDICTION *************/

	
		char c;
		scanf("%c",&c);

	}  catch(exception& e) {
		cout << e.what() << "\n";
	}

	return 0;
}