#include "flower_detection.h"

/*
*	This program trains our algorithm.
*	It uses samples directory, which contains images for training, in the following format:
*	<TRAIN_EXE_PATH>/samples/<FLOWER_ID>/<FLOWER_PIC_NUM>.jpg
*	Each folder must contain a centers.txt file in the following format:
*	for each flower photo;
*	<FLOWER_PIC_NUM>,<CENTER_POINT_X>,<CENTER_POINT_Y>
*	program output is stored under serialization/svmSerialize.txt file.
*/
int main(int argc, char* argv[]) 
{
	SvmModel svm;
	TrainSet trainSet;	

	FILE* centers=NULL;
	Point center;
	int center_res;
	cout << "\n\n Samples:\n";
	
	//loop over files
	/******** DIR SCANNING CODE *************/
	int flowerID=0;
	int flowerPicNo=0;
	DIR *dir,*subdir, *temp;
	struct dirent *ent;

	/* open directory stream */
	dir = Utils::MyChdir(Consts::SAMPLES_DIR_PATH);

	if (dir == NULL)
		return -1;
	
	ent = readdir(dir);
	/* loop all over the files and directories within directory */
	while (ent != NULL) {
		switch (ent->d_type) {

		case DT_REG:
			Utils::toLowercase(ent->d_name + strlen(ent->d_name)-4);
			if(strcmp(".jpg",ent->d_name + strlen(ent->d_name)-4)!=0){
				break;
			}
			sscanf(ent->d_name, "%d.jpg", &flowerPicNo);
			center_res = Utils::getCenter(center,centers,flowerPicNo);
			//printf ("No = %d, Filename: %s\n", flowerPicNo, temp_filename);
			if(center_res==0)
			{
				try 
				{
					FlowerFeatureExtractor flower_featurs(ent->d_name, center);
					flower_featurs.extractFeatures();
					flower_featurs.m_sample.m_label = flowerID;//CHANGE TO FLOWER LABEL
					trainSet.addSample(flower_featurs.m_sample);

				} catch(exception& e) {
					printf("Image %s dropped!\n", ent->d_name);					
				}
					
			}
			else
				printf("Image %s dropped!\n", ent->d_name);
			break;

		case DT_DIR:
			if(strcmp(ent->d_name,".")==0 || strcmp(ent->d_name,"..")==0)
				break;

			subdir = Utils::MyChdir(ent->d_name);

			temp = subdir;
			subdir = dir;
			dir = temp;

			if(!PathFileExists(Consts::CENTERS_FILE_NAME.c_str())){
				printf("ERROR: centers file does not exist in directory - %s. Please add it before continue training.\n",ent->d_name);
				return -1;
			}
			centers = fopen(Consts::CENTERS_FILE_NAME.c_str(),"r");
			if(centers == NULL){
				printf("ERROR: couldnt open centers file in directory - %s.\n",ent->d_name);
				return -1;
			}
			sscanf(ent->d_name, "%d", &flowerID);
			//printf ("ID = %d\n", flowerID);
			break;
		}

		if((ent = readdir(dir)) == NULL){
			dir = subdir;
			fclose(centers);
			_chdir("..");
			ent = readdir (dir);
			printf("\n");
		}
	}
	closedir (dir);

	/******** END DIR SCANNING CODE *************/
	//trainSet.toCSV("trainSetCsv.txt");
	svm.train(trainSet);
	
	if(_chdir(Consts::SVM_SERIALIZATION_DIR_PATH.c_str())!=0){
		printf("ERROR: directory didn't exist, creating it - %s\n", Consts::SVM_SERIALIZATION_DIR_PATH);
		_mkdir(Consts::SVM_SERIALIZATION_DIR_PATH.c_str());
		if(_chdir(Consts::SVM_SERIALIZATION_DIR_PATH.c_str())!=0){
			printf("ERROR: couldnt create directory - %s\n",Consts::SVM_SERIALIZATION_DIR_PATH);
			return -1;
		}
	}
	svm.save(Consts::SVM_SERIALIZE_FILE_NAME);

	return 0;
}

