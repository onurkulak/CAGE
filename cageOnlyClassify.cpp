#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <queue>
/* Total number of dimensions allowed after dimension reduction is sqrt(number of areas)*10. In the LP the summation of weights of dimensions should be less than half of this value*/

FILE *fpOut;

typedef struct samples
{
    char sampleName[50];
    int caseControl; // 1: case, 0: control
    int trainOrTest; //1: test, 0: train
    double* values;
    double distanceToCenter;
    int filterSample;
    bool covered;
} samples;

samples* listSamples;
int totalSamples;
int numberGene;
double** distanceMatrix;



int maxCloRow;
double* weightsCalculated;



typedef struct areas
{
    int* caseIDs;
    int numCases;
    int centerId;
} areas;


int compare (const void * a, const void * b)
{

    double xx = *(double*)a;
    double yy = *(double*)b;
    if (xx<yy) return -1;
    if (xx>yy) return 1;
    return 0;

}

void updateDistanceMatrix(double *geneWeights)
{
    for (int count=0; count<totalSamples; count++)
    {
        for (int count2=count+1; count2<totalSamples; count2++)
        {

            double distance=0;
            for (int i=0; i<numberGene; i++)
            {
                distance=distance+geneWeights[i]*fabs(listSamples[count].values[i]-listSamples[count2].values[i]);
            }
            distanceMatrix[count][count2]=distance;
            distanceMatrix[count2][count]=distance;
        }
    }
}


double* getRadiuses()
{
    double* radisusControls = new double[totalSamples];
    double* radiusOfArea = new double[totalSamples];
    for (int count=0; count<totalSamples; count++)
    {
        int numControls=0;
        if ((listSamples[count].caseControl==1) && (listSamples[count].filterSample==0) &&
                (listSamples[count].trainOrTest==0))
        {
            for (int tempC=0; tempC<totalSamples; tempC++)
            {
                radisusControls[tempC]=0;
            }

            for (int count2=0; count2<totalSamples; count2++)
            {
                if ((listSamples[count2].caseControl==0) && (listSamples[count2].filterSample==0) && (listSamples[count2].trainOrTest==0))
                {
                    radisusControls[numControls]=distanceMatrix[count][count2];
                    numControls++;
                }
            }
            qsort(radisusControls, numControls, sizeof(double), compare);
            radiusOfArea[count]=radisusControls[0];
        }
    }
    delete[] radisusControls;
    return radiusOfArea;
}


int* getCaseNeighbors(double* radiusOfArea, std::vector<int>* &neighborPts, bool epsilonIsMinRadius){
    int* casesInAreaCenter = new int[totalSamples];


    neighborPts = new std::vector<int>[totalSamples];

    for (int count=0; count<totalSamples; count++)
            casesInAreaCenter[count]=0;

    for (int count=0; count<totalSamples; count++)
        {
            if ((listSamples[count].caseControl==1) && (listSamples[count].filterSample==0) &&
                    (listSamples[count].trainOrTest==0))
            {
                for (int count2=0; count2<=count; count2++)
                {
                    if ((listSamples[count2].caseControl==1)  &&
                            (listSamples[count2].filterSample==0) &&
                            (listSamples[count2].trainOrTest==0))
                    {
                        if(epsilonIsMinRadius){
                            double distanceMin = std::min(radiusOfArea[count], radiusOfArea[count2]);
                            if(distanceMatrix[count][count2]<distanceMin){
                                casesInAreaCenter[count]++;
                                casesInAreaCenter[count2]++;
                                neighborPts[count].push_back(count2);
                                neighborPts[count2].push_back(count);
                            }
                        }
                        else{
                            if(distanceMatrix[count][count2]<radiusOfArea[count]){
                                casesInAreaCenter[count]++;
                                neighborPts[count].push_back(count2);
                            }
                            if(distanceMatrix[count][count2]<radiusOfArea[count2]){
                                casesInAreaCenter[count2]++;
                                neighborPts[count2].push_back(count);
                            }
                        }
                    }
                }
            }
        }
    return casesInAreaCenter;
}

double** getCenterWeights(int numberOfAreas,bool averagedCenters, int* numCase, int* centerIds, int** casesToConsider){
    double** centerWeights = new double*[numberOfAreas];
    if(averagedCenters){
        for(int i = 0; i < numberOfAreas; i++)
            {
                if(numCase[i]==1)
                    centerWeights[i] = listSamples[centerIds[i]].values;
                else{
                    centerWeights[i] = new double[numberGene];
                    for(int j = 0; j < numberGene; j++){
                        double sum = 0;
                        for(int k = 0; k < numCase[i]; k++)
                            sum+=listSamples[casesToConsider[i][k]].values[j];
                        centerWeights[i][j] = sum/numCase[i];
                    }
                }
            }
    }
    else{
        for(int i = 0; i < numberOfAreas; i++)
            centerWeights[i] = listSamples[centerIds[i]].values;
    }
    return centerWeights;
}



std::vector<areas> getClusters(int minpts, double *geneWeights, bool epsilonIsMinRadius, int* &casesInAreaCenter)
{
    /*method explanation:
    if epsilonIsMinRadius, epsilon value between two cases is small one of their radiuses, both should cover eachother to extend cluster
    else, to expand a cluster from a point 'P' it is enough that its radius covers another point, and number of points covered this way is >=minpts
    */
    std::vector<areas> listAreas;
    updateDistanceMatrix(geneWeights);
    double* radiusOfArea =  getRadiuses();
    std::vector<int>* neighborPts;

    casesInAreaCenter = getCaseNeighbors(radiusOfArea, neighborPts, epsilonIsMinRadius);
    for(int i=0; i<totalSamples; i++)
        listSamples[i].covered = false;


    for(int i = 0; i < totalSamples; i++){
        if(casesInAreaCenter[i]>=minpts && !listSamples[i].covered){
            int maxPoints = 0;
            int centerIndex = -1;
            std::queue<int> addingOrder;
            std::vector<int> ids;

            addingOrder.push(i);
            while(!addingOrder.empty()){
                int id = addingOrder.front();
                addingOrder.pop();

                if(!listSamples[id].covered){
                    listSamples[id].covered = true;
                    ids.push_back(id);
                    if(casesInAreaCenter[id]>=minpts){
                        for(int j = 0; j < casesInAreaCenter[id]; j++)
                            addingOrder.push(neighborPts[id][j]);
                        if(casesInAreaCenter[id]>maxPoints){
                            maxPoints = casesInAreaCenter[id];
                            centerIndex = id;
                        }
                    }
                }
            }

            areas a;
            a.numCases = ids.size();
            a.caseIDs = new int[ids.size()];
            a.centerId = centerIndex;
            std::copy(ids.begin(), ids.end(), a.caseIDs);
            listAreas.push_back(a);
        }

    }

    delete[] neighborPts;
    delete[] radiusOfArea;
    delete[] casesInAreaCenter;
    return listAreas;


}



int main(int argv, char *argc[])
{
    FILE *fp=fopen(argc[1],"r");
    numberGene=atoi(argc[2]);
    totalSamples=atoi(argc[3]);
    FILE *fp2=fopen(argc[4],"r");// Every samples name grouped as training (0) or testing sample (1)
    FILE *fp4=fopen(argc[5],"r");//Samples to remove
    FILE *preparedWeights = fopen(argc[6],"r");
    int minpts = atoi(argc[7]);

    bool averagedCenters = atoi(argc[8])==1;
    bool epsilonIsMinRadius = atoi(argc[9])==1;
    int lpType = atoi(argc[10]);

    int numberOfAreas = 0;
    maxCloRow = numberGene+1;
    weightsCalculated = new double[maxCloRow];


    distanceMatrix = new double*[totalSamples];
    for (int count=0; count<totalSamples; count++)
    {
        distanceMatrix[count] = new double[totalSamples];
    }



    int sampleId;
    int trainOrTest;
    char sampleName[50];

    
    listSamples = new samples[totalSamples];
    for (int count=0; count<totalSamples; count++)
    {
        fscanf(fp, "%s\t%i", listSamples[count].sampleName, &(listSamples[count].caseControl));
        listSamples[count].filterSample=0;
        listSamples[count].values = new double[numberGene];
        for (int count2=0; count2<numberGene; count2++)
        {
            fscanf(fp,"\t%lf", &(listSamples[count].values[count2]));
            //printf("%lf\n", listSamples[count].values[count2]);
        }
        fscanf(fp,"\n");

    }




    while(fscanf(fp2,"%s\t%i\n", sampleName, &trainOrTest)!=EOF)
    {
        for (int count=0; count<totalSamples; count++)
        {
            if (strcmp(listSamples[count].sampleName, sampleName)==0)
            {
                listSamples[count].trainOrTest=trainOrTest;
            }
        }
    }



    while(fscanf(fp4,"%s\n", sampleName)!=EOF)
    {
        for (int count=0; count<totalSamples; count++)
        {
            if (strcmp(listSamples[count].sampleName, sampleName)==0)
            {
                listSamples[count].filterSample=1;
            }
        }
    }

    int areaCount;
    fscanf(preparedWeights,"%i\n",&areaCount);
    for (int i=0; i<numberGene; i++)
    {
        fscanf(preparedWeights,"%lf ", weightsCalculated+i);
    }
    int* coreN;
    std::vector<areas> bestAreas = getClusters(minpts ,weightsCalculated, epsilonIsMinRadius,coreN);
    int bestNumberOfAreas = bestAreas.size();

    updateDistanceMatrix(weightsCalculated);
    double* oneLastFinal = getRadiuses();
    double* radisusControls = new double[totalSamples];


    int totalPredictions = 0, trueCase = 0, trueControl = 0, falseCase = 0, falseControl = 0; 



    for(int i = 0; i < totalSamples; i++)
        if(listSamples[i].trainOrTest==1 && listSamples[i].filterSample==0){
            totalPredictions++;
            double radius = 0;
            if(epsilonIsMinRadius)
            {
            int numControls = 0;
            for (int tempC=0; tempC<totalSamples; tempC++)
            {
                radisusControls[tempC]=0;
            }

            for (int count2=0; count2<totalSamples; count2++)
            {
                if ((listSamples[count2].caseControl==0) && (listSamples[count2].filterSample==0) && (listSamples[count2].trainOrTest==0))
                {
                    radisusControls[numControls]=distanceMatrix[i][count2];
                    numControls++;
                }
            }
            qsort(radisusControls, numControls, sizeof(double), compare);
            radius=radisusControls[0];
            }

            bool gotcha = false;
            for(int j = 0; j < bestNumberOfAreas; j++){
                for(int k = 0; k < bestAreas[j].numCases; k++)
                {
                    if(coreN[bestAreas[j].caseIDs[k]]>=minpts)
                    if(epsilonIsMinRadius)
                        gotcha = std::min(oneLastFinal[bestAreas[j].caseIDs[k]], radius) > 
                            distanceMatrix[i][bestAreas[j].caseIDs[k]];
                    gotcha = oneLastFinal[bestAreas[j].caseIDs[k]] > distanceMatrix[i][bestAreas[j].caseIDs[k]];

                    if(gotcha){
                        printf("%s\tis guessed to be a case around the area:\t%d\n", listSamples[i].sampleName, j);
                        break;
                    }
                }
                if(gotcha)
                    break;
            }
            if(!gotcha)
                printf("%s\tis guessed to be a control\n", listSamples[i].sampleName);
            
            if(listSamples[i].caseControl==0){
                printf("%s\tis a control\n", listSamples[i].sampleName);
                if(gotcha)
                    falseCase++;
                else trueControl++;
            }
            else {
                printf("%s\tis a case\n", listSamples[i].sampleName);
                if(gotcha)
                    trueCase++;
                else falseControl++;
            }
        }

        delete[] radisusControls;
        delete[] oneLastFinal;

        printf("Out of total %d predictions:\n", totalPredictions);
        printf("%d is guessed to be case, %d was control\n",  falseCase+trueCase, falseCase);
        printf("%d is guessed to be control, %d was case\n", falseControl+trueControl, falseControl);
        printf("false prediction: %f", falseCase/((double)falseCase+trueCase));

}




