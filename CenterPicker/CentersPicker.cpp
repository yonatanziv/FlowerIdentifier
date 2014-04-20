#include "utils.h"

#define SVM_NUM_SAMPLES 100
#define SAMPLES_DIR_PATH "flowers"
#define SVM_SERIALIZATION_DIR_PATH "serialization"
#define CENTERS_FILE_NAME "centers.txt"
#define ALL_CENTERS_FILE_NAME "All_Centers.txt"


FILE * centersFile;
FILE* allCentersFile;
double ratio;
int flowerPicNo =0;
int counter = 70;

void my_mouse_callback( int event, int x, int y, int flags, void* param ){

	if(event == CV_EVENT_LBUTTONUP){
		printf("(%d,%d)\n",(int)(x*ratio),(int)(y*ratio));

		if(flowerPicNo<1000) 
		{
			fprintf(centersFile,"%0d,%d,%d\n",flowerPicNo,(int)(x*ratio),(int)(y*ratio));
			fprintf(allCentersFile,"%0d,%d,%d\n",flowerPicNo,(int)(x*ratio),(int)(y*ratio));
		}
		else 
		{
			fprintf(centersFile,"%d,%d,%d\n",flowerPicNo,(int)(x*ratio),(int)(y*ratio));
			fprintf(allCentersFile,"%d,%d,%d\n",flowerPicNo,(int)(x*ratio),(int)(y*ratio));
		}
	}

}

int main(int argc, char* argv[]) 
{		
	//loop over files
	/******** DIR SCANNING CODE *************/
	int flowerID=0;
	DIR *dir,*subdir, *temp;
	struct dirent *ent;

	allCentersFile = fopen(ALL_CENTERS_FILE_NAME,"w+");

	/* open directory stream */
	dir = Utils::MyChdir(SAMPLES_DIR_PATH);

	if (dir == NULL)
		return -1;
	
	const char* name = "Box Example";

	Mat image;
	Mat img_resized;

	ent = readdir(dir);
	/* print all the files and directories within directory */
	while (ent != NULL) {
		switch (ent->d_type) {

		case DT_REG:
			Utils::toLowercase(ent->d_name + strlen(ent->d_name)-4);
			if(strcmp(".jpg",ent->d_name + strlen(ent->d_name)-4)!=0){
				break;
			}

			char newName[10];
			int index;
			for(index = 0; index<10;index++)
				newName[index] = '\0';
			
			index = 0;
			if(flowerID<10)
				newName[index++] = '0';

			itoa(flowerID,newName+index,10);
			index += (int)(log10((double)flowerID));
			index++;

			if(counter<10)
				newName[index++] = '0';

			itoa(counter,newName+index,10);
			index += log10((double)flowerID);
			index++;

			strcat(newName,".jpg");

			if(rename( ent->d_name, newName )!=0){
				perror("can't rename file");
				return -1;
			}
			
			counter++;

			sscanf(newName, "%d.jpg", &flowerPicNo);

			image = imread(newName);
			Utils::resizeImage(image, img_resized);


			ratio = image.rows*1.0/img_resized.rows;
			
			namedWindow(name);

			Utils::myShowImage(name, img_resized );
			printf("%d: ",counter);

			// Set up the callback
			setMouseCallback(name, my_mouse_callback);

			cvWaitKey(0);
			break;

		case DT_DIR:
			if(strcmp(ent->d_name,".")==0 || strcmp(ent->d_name,"..")==0)
				break;

			subdir = Utils::MyChdir(ent->d_name);

			temp = subdir;
			subdir = dir;
			dir = temp;

			centersFile = fopen(CENTERS_FILE_NAME,"w+");
			sscanf(ent->d_name, "%d", &flowerID);
			//printf ("ID = %d\n", flowerID);
			break;
		}

		if((ent = readdir(dir)) == NULL){
			dir = subdir;
			fclose(centersFile);
			counter = 70;
			chdir("..");
			ent = readdir (dir);
			printf("\n");
		}
	}
	closedir (dir);


	char c;
	scanf("%c",&c);

	return 0;

}

