﻿//version 54
#include "flower_detection_header.h"
#include "utils.h"

#define SVM_NUM_SAMPLES 5000
#define SAMPLES_DIR_PATH "samples"
#define SVM_SERIALIZATION_DIR_PATH "serialization"
#define CENTERS_FILE_NAME "centers.txt"
#define SVM_SERIALIZE_FILE_NAME "svmSerialize.txt"

/*
*	This program trains our algorithm.
*	it uses samples directory, which contains images for training, in the following format:
*	<TRAIN_EXE_PATH>/samples/<FLOWER_ID>/<FLOWER_PIC_NUM>.jpg
*	and each flower, its folder must contain a centers.txt file in the following format:
*	for each flower photo;
*	<FLOWER_PIC_NUM>,<CENTER_POINT_X>,<CENTER_POINT_Y>
*	program output is stored under serialization/svmSerialize.txt file.
*/
int main(int argc, char* argv[]) 
{
	DataSample current_sample;
	SvmModel svm;
	TrainSet trainSet;
	FlowerFeatureExtractor flower_featurs;

	FILE* centers=NULL;
	Point center;
	int center_res;
	printf("\n\n Samples:\n");
	
	//loop over files
	/******** DIR SCANNING CODE *************/
	int flowerID=0;
	int flowerPicNo=0;
	DIR *dir,*subdir, *temp;
	struct dirent *ent;

	/* open directory stream */
	dir = Utils::MyChdir(SAMPLES_DIR_PATH);

	if (dir == NULL)
		return -1;
	
	ent = readdir(dir);
	/* print all the files and directories within directory */
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
					current_sample.m_label=flowerID; //CHANGE TO FLOWER LABEL
					flower_featurs.extractFeaturesFromImage(ent->d_name, center);
					trainSet.addSample(current_sample);

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

			if(!PathFileExists(CENTERS_FILE_NAME)){
				printf("ERROR: centers file does not exist in directory - %s. Please add it before continue training.\n",ent->d_name);
				return -1;
			}
			centers = fopen(CENTERS_FILE_NAME,"r");
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

	svm.train(trainSet);
	
	if(_chdir(SVM_SERIALIZATION_DIR_PATH)!=0){
		printf("ERROR: directory didn't exist, creating it - %s\n",SVM_SERIALIZATION_DIR_PATH);
		_mkdir(SVM_SERIALIZATION_DIR_PATH);
		if(_chdir(SVM_SERIALIZATION_DIR_PATH)!=0){
			printf("ERROR: couldnt create directory - %s\n",SVM_SERIALIZATION_DIR_PATH);
			return -1;
		}
	}
	svm.save(SVM_SERIALIZE_FILE_NAME);

	return 0;
}
