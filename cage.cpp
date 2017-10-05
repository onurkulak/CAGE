#include <ilcplex/cplex.h>
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

void deallocateCenterWeights(double** centerWeights, int numberOfAreas, bool averagedCenters, int* numCase){
    if(averagedCenters){
        for(int i = 0; i < numberOfAreas; i++)
            if(numCase[i] > 1)
                delete[] centerWeights[i];
    }
    delete[] centerWeights;
}


void createLPEquation(double* obj, int numberOfAreas,int controlCount,int* numCase,int** casesToConsider, int lpType) {
    int totalCasesInAllAreas = 0;
    for (int i=0; i<numberOfAreas; i++)
        totalCasesInAllAreas+=numCase[i];
    
    for (int w=0; w<numberGene; w++)
        obj[w] = 0;
    double coefficient1, coefficient2, coefficient3;

    for(int sph = 0; sph<numberOfAreas; sph++)
    {
        int caseTimesControl = controlCount * numCase[sph];

        for (int w=0; w<numberGene; w++)
        {

            /*
            lp according to lpType:
            if lpType is 0 (log(1+dist))/count-1 applied for each term
            if lpType is 1 arithmetic means are used, like in the original odin
            if lpType is 2 geometric means are used of: (1+dist) to prevent 0 multipliers
            if lpType is 3 for area squeezing arithmetic mean is used, for control and other area's distancing geometric means are used
            */
            coefficient1=0;
            coefficient2=0;
            coefficient3=0;


            for (int cases=0; cases<numCase[sph]; cases++)
            {
                if(lpType == 1 || lpType == 3)
                    for(int cases2 = 0; cases2 < cases; cases2++)
                        coefficient1+=fabs(listSamples[casesToConsider[sph][cases]].values[w]-listSamples[casesToConsider[sph][cases2]].values[w]);
                else 
                    for(int cases2 = 0; cases2 < cases; cases2++)
                        coefficient1+=log(1+fabs(listSamples[casesToConsider[sph][cases]].values[w]-listSamples[casesToConsider[sph][cases2]].values[w]));
            }

            for (int samples=0; samples<totalSamples; samples++)
            {
                if ((listSamples[samples].caseControl==0) && (listSamples[samples].trainOrTest==0) && (listSamples[samples].filterSample==0))
                {
                    if(lpType == 1)
                        for (int cases=0; cases<numCase[sph]; cases++)
                            coefficient3+=fabs(listSamples[casesToConsider[sph][cases]].values[w]-listSamples[samples].values[w]);
                    else
                        for (int cases=0; cases<numCase[sph]; cases++)
                            coefficient3+=log(1+fabs(listSamples[casesToConsider[sph][cases]].values[w]-listSamples[samples].values[w]));
                }
            }

            int totalCasesToCasesCount = (numCase[sph]) * (totalCasesInAllAreas-numCase[sph]);

            for(int caseCounter = 0; caseCounter < numCase[sph]; caseCounter++)
                for(int i = 0; i<numberOfAreas; i++)
                    if(sph!=i)
                    {
                        if(lpType == 1)
                        for(int j = 0; j<numCase[i]; j++)
                                coefficient2+=fabs(listSamples[casesToConsider[sph][caseCounter]].values[w]-
                                              listSamples[casesToConsider[i][j]].values[w]);
                        else 
                            for(int j = 0; j<numCase[i]; j++)
                                coefficient2+=log(1+fabs(listSamples[casesToConsider[sph][caseCounter]].values[w]-
                                              listSamples[casesToConsider[i][j]].values[w]));
                    }


            
            int casesTimesCases = numCase[sph]*(numCase[sph]-1)/2;

            double m = numCase[sph];
            double v = 0;
            if(lpType==1 || lpType == 0){
                if(casesTimesCases!=0)
                    v+=coefficient1/casesTimesCases;
                if(caseTimesControl!=0)
                    v-=coefficient3/caseTimesControl;
                if(totalCasesToCasesCount!=0)
                    v-=(coefficient2/totalCasesToCasesCount);
            }
            else if(lpType == 2){
                if(casesTimesCases!=0)
                    v+=exp(coefficient1/casesTimesCases);
                if(caseTimesControl!=0)
                    v-=exp(coefficient3/caseTimesControl);
                if(totalCasesToCasesCount!=0)
                    v-=exp((coefficient2/totalCasesToCasesCount));
            }
            else{
                if(casesTimesCases!=0)
                    v+=coefficient1/casesTimesCases;
                if(caseTimesControl!=0)
                    v-=exp(coefficient3/caseTimesControl);
                if(totalCasesToCasesCount!=0)
                    v-=exp((coefficient2/totalCasesToCasesCount));
            }
            obj[w]+=m*v;
        }
    }
}

int calculateWeights(int numberOfAreas, int **casesToConsider, int *numCase, int dimensionCount, int* centerIds, 
    bool averagedCenters, int lpType, bool noConstraints)
{

    int totalCasesInAllAreas = 0;
    for (int i=0; i<numberOfAreas; i++)
        totalCasesInAllAreas+=numCase[i];

    int solstat;
    double objval;
    double *x = NULL;
    double *pi = NULL;
    double *slack = NULL;
    double *dj = NULL;
    int cur_numrows, cur_numcols;
    int status=0;
    int NumRows, NumCols, NumNZ;

    //if averaged center option is chosen, instead of choosing the case with most neigbors for center weights
    //takes average of all cases for a cluster and these averages forms new weight center, it's not an actual point!
    double** centerWeights = getCenterWeights(numberOfAreas, averagedCenters, numCase, centerIds, casesToConsider);
    
    double *obj = new double[maxCloRow];
    double *lb= new double[maxCloRow];
    double *ub= new double[maxCloRow];
    char **colname= new char*[maxCloRow];

    char *strings;

    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    env = CPXopenCPLEX (&status);
    status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
    status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
    lp = CPXcreateprob (env, &status, "lpex1");


    status = CPXchgobjsen (env, lp, CPX_MIN);

    int controlCount = 0;

    for (int samples=0; samples<totalSamples; samples++)
    {
        if ((listSamples[samples].caseControl==0) && (listSamples[samples].trainOrTest==0) && (listSamples[samples].filterSample==0))
            controlCount++;
    }

    int maxConstraintCount = numberOfAreas*controlCount+totalCasesInAllAreas+1;
    int *rmatbeg = new int[maxConstraintCount];
    int *rmatind = new int[maxConstraintCount*maxCloRow-1];
    double *rmatval = new double[maxConstraintCount*maxCloRow-1];
    double *rhs = new double[maxConstraintCount];
    char *sense = new char[maxConstraintCount];
    char **rowname = new char*[maxConstraintCount];

    createLPEquation(obj, numberOfAreas, controlCount, numCase, casesToConsider, lpType);

    for (int w=0; w<numberGene; w++)
    {
        strings=(char *)malloc(100*sizeof(char));
        
        obj[w] = obj[w]/((double)totalCasesInAllAreas);
        lb[w]=0;
        ub[w]=1;

        sprintf(strings, "w%i\0", w);
        colname[w]=strings;

    }


    obj[numberGene]=0;
    lb[numberGene]=0;
    ub[numberGene]=CPX_INFBOUND;

    colname[numberGene]="r\0";
    status = CPXnewcols (env, lp, numberGene+1, obj, lb, ub, NULL, colname);


    int rmatbegCount=0;
    int constraintCount=0;
    if(!noConstraints){
        for(int sph = 0; sph<numberOfAreas; sph++)
            for (int cases=0; cases<numCase[sph]; cases++)
            {
                rmatbeg[constraintCount]=rmatbegCount;
                strings=(char *)malloc(100*sizeof(char));
                sprintf(strings, "c%i\0", constraintCount);
                rowname[constraintCount]=strings;
                for (int w=0; w<numberGene; w++)
                {
                    rmatind[rmatbegCount]=w;
                    rmatval[rmatbegCount]=fabs(listSamples[casesToConsider[sph][cases]].values[w]-centerWeights[sph][w]);
                    rmatbegCount++;
                }
    
                rmatind[rmatbegCount]=numberGene;
                rmatval[rmatbegCount]=-1;
                rmatbegCount++;
    
                sense[constraintCount]='L';
                rhs[constraintCount]=0;
    
    
                constraintCount++;
            }
    
        
            for (int allSamples=0; allSamples<totalSamples; allSamples++)
            {
                if ((listSamples[allSamples].caseControl==0) && (listSamples[allSamples].trainOrTest==0))// it is a controls
                    for(int sph = 0; sph<numberOfAreas; sph++){
                        {
                            //printf("ConstCount %i %i\n", constraintCount, rmatbegCount);
                            rmatbeg[constraintCount]=rmatbegCount;
                            strings=(char *)malloc(100*sizeof(char));
                            sprintf(strings, "d%i\0", constraintCount);
        
                            rowname[constraintCount]=strings;
                            for (int w=0; w<numberGene; w++)
                            {
                                rmatind[rmatbegCount]=w;
                                rmatval[rmatbegCount]=-1*fabs(listSamples[allSamples].values[w]-centerWeights[sph][w]);
                                rmatbegCount++;
                            }
        
                            rmatind[rmatbegCount]=numberGene;
                            rmatval[rmatbegCount]=1;
                            rmatbegCount++;
        
                            sense[constraintCount]='L';
                            rhs[constraintCount]=0;
                            constraintCount++;
                        }
                }
                
            }
        }

    rmatbeg[constraintCount]=rmatbegCount;
    rowname[constraintCount]="t";
    for (int w=0; w<numberGene; w++)
    {
        rmatind[rmatbegCount]=w;
        rmatval[rmatbegCount]=1;
        rmatbegCount++;
    }
    sense[constraintCount]='L';
    rhs[constraintCount]=dimensionCount/2;
    constraintCount++;


    status = CPXaddrows (env, lp, 0, constraintCount, rmatbegCount, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);

    status = CPXlpopt (env, lp);

    cur_numrows = CPXgetnumrows (env, lp);
    cur_numcols = CPXgetnumcols (env, lp);


    x = (double *) malloc (cur_numcols * sizeof(double));
    slack = (double *) malloc (cur_numrows * sizeof(double));
    dj = (double *) malloc (cur_numcols * sizeof(double));
    pi = (double *) malloc (cur_numrows * sizeof(double));
    status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);


    for (int count=0; count<cur_numcols; count++)
    {
        if (x[count]>0)
            printf("Col %i %lf\n", count, x[count]);
        weightsCalculated[count] = x[count];
    }
    

    delete[] obj;
    delete[] lb;
    delete[] ub;
    delete[] colname;
    delete[] rmatbeg;
    delete[] rmatind;
    delete[] rmatval;
    delete[] sense;
    delete[] rowname;
    delete[] rhs;

    free(x);
    free(slack);
    free(dj);
    free(pi);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX (&env);
    deallocateCenterWeights(centerWeights, numberOfAreas, averagedCenters, numCase);
}


std::vector<areas> getClusters(int minpts, double *geneWeights, bool epsilonIsMinRadius)
{
    /*method explanation:
    if epsilonIsMinRadius, epsilon value between two cases is small one of their radiuses, so both should cover eachother to extend cluster.
    else, to expand a cluster from a point 'P' it is enough that its radius covers another point, and number of points covered this way is >=minpts
    */
    std::vector<areas> listAreas;
    updateDistanceMatrix(geneWeights);
    double* radiusOfArea =  getRadiuses();
    std::vector<int>* neighborPts;

    int* casesInAreaCenter = getCaseNeighbors(radiusOfArea, neighborPts, epsilonIsMinRadius);
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

/*bool continueLoop(int c, int cc, int cp, int ac, int ap, int method, double cw, double pw){
    if(c>10)
        return false;
    if(c<3)
        return true;
    if(0==method)
        return cc>=cp||ac>=ap;
    else if(method==1) 
        return cw>=pw;
    else return cc>=cp||ac>=ap||cw>=pw;
}*/
double getDavies_BouldinVal(std::vector<areas> baloon, double* weights){
    int numberOfAreas = baloon.size();
    int** caseIDsInAreas = new int*[numberOfAreas];
    int* caseCounts = new int[numberOfAreas];
    int* centerIds = new int[numberOfAreas];
    for (int i=0; i<numberOfAreas; i++)
    {
        centerIds[i] = baloon[i].centerId;
        caseIDsInAreas[i] = baloon[i].caseIDs;
        caseCounts[i] = baloon[i].numCases;
    }
    double** pCenters = getCenterWeights(numberOfAreas, true, caseCounts, centerIds, caseIDsInAreas);

    double* centroidClusterDistanceAverages = new double[numberOfAreas]();
    double** centroidDistances = new double*[numberOfAreas];

    for(int i = 0; i < numberOfAreas; i++){
        double* avgDistArray = new double[numberGene]();
        for(int j = 0; j < numberGene; j++){
            for(int k = 0; k < caseCounts[i]; k++)
                avgDistArray[j]+=fabs(pCenters[i][j]-listSamples[caseIDsInAreas[i][k]].values[j]);
            centroidClusterDistanceAverages[i] += avgDistArray[j]*weights[j];
        }
        delete[] avgDistArray;
        centroidClusterDistanceAverages[i]/=caseCounts[i];
        centroidDistances[i] = new double[i]();
        for(int j = 0; j < i; j++)
            for(int k = 0; k < numberGene; k++)
                centroidDistances[i][j] += fabs(pCenters[i][k]-pCenters[j][k])*weights[k];
    }

    double valSum = 0;
    for(int i = 0 ; i < numberOfAreas; i++){
        double maxVal = 0;
        for(int j = 0; j < numberOfAreas; j++)
        {
            if(i==j)
                continue;
            double val = centroidClusterDistanceAverages[i] + centroidClusterDistanceAverages[j];
            if(j < i)
                val/=centroidDistances[i][j];
            else 
                val/= centroidDistances[j][i];
            if(maxVal<val)
                maxVal = val;
        }
        valSum += maxVal;
    }
    deallocateCenterWeights(pCenters, numberOfAreas, true, caseCounts);
    delete[] centroidClusterDistanceAverages;
    delete[] centerIds;
    delete[] caseIDsInAreas;
    delete[] caseCounts;
    for(int i = 0; i < numberOfAreas; i++)
        delete[] centroidDistances[i];
    delete[] centroidDistances;
    return valSum/numberOfAreas;
}


bool continueLoop(int c, std::vector<areas> preList, std::vector<areas> curList, double* preWeights, bool increaseInCoverage){
    if(c>10)
        return false;
    if(c<2)
        return true;
    if(curList.size()>preList.size()||increaseInCoverage)
        return true;

    return getDavies_BouldinVal(curList,weightsCalculated) < getDavies_BouldinVal(preList,preWeights);
}


int main(int argv, char *argc[])
{
    FILE *fp=fopen(argc[1],"r");
    numberGene=atoi(argc[2]);
    totalSamples=atoi(argc[3]);
    FILE *fp2=fopen(argc[4],"r");// Every samples name grouped as training (0) or testing sample (1)
    FILE *fp4=fopen(argc[5],"r");//Samples to remove
    fpOut = fopen(argc[6],"w");
    int minpts = atoi(argc[7]);

    bool averagedCenters = atoi(argc[8])==1;
    bool epsilonIsMinRadius = atoi(argc[9])==1;
    int lpType = atoi(argc[10]);
    int dimensionParameter = atoi(argc[11]);
    int iterationUntilWeight = atoi(argc[12]);
    bool noConstraints = atoi(argc[13])==1;
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



    for (int i=0; i<numberGene; i++)
    {
        weightsCalculated[i]=1;
    }


    double *tempWeights = new double[numberGene];

    int numCasesTotalCur=1;
    int numCasesTotalPre=0;
    int numCasesAverageCur=1;
    int numCasesAveragePre=0;
    double previousWeight = 0;
    double currentWeight = 0;

    int countX=0;

    std::vector<areas> listAreas;
    std::vector<areas> bestAreas;
    int bestNumberOfAreas;
    double* bestWeights = new double[numberGene];
    bool cont = true;
    while(cont)
    {
        numCasesTotalPre = numCasesTotalCur;
        numCasesAveragePre = numCasesAverageCur;
        bestNumberOfAreas = numberOfAreas;

        if(!bestAreas.empty())
        {
            for (int i=0; i<bestAreas.size(); i++)
                delete[] bestAreas[i].caseIDs;
        }

        bestAreas = listAreas;
        listAreas = getClusters(minpts ,weightsCalculated, epsilonIsMinRadius);

        numberOfAreas = listAreas.size();
        numCasesTotalCur = 0;

        for (int i=0; i<numberOfAreas; i++)
        {
            numCasesTotalCur+= listAreas[i].numCases;
            printf("\n\nCases Covered %d:\n", listAreas[i].numCases);
            for (int j = 0; j < listAreas[i].numCases; ++j)
            {
                printf("Samples Being Covered %i %s\n", countX, listSamples[listAreas[i].caseIDs[j]].sampleName);
            }

        }

        numCasesAverageCur=numCasesTotalCur/numberOfAreas;


        int numberOfMaximumDimensions;
        if(dimensionParameter==0)
            numberOfMaximumDimensions = (int)(2*pow(numCasesTotalCur,0.5));
        else
            numberOfMaximumDimensions = 5*(dimensionParameter+1);
        cont = continueLoop(countX, bestAreas, listAreas, bestWeights, numCasesTotalCur>numCasesTotalPre);

        if(cont)
        {
            int** caseIDsInAreas = new int*[numberOfAreas];
            int* caseCounts = new int[numberOfAreas];
            int* centerIds = new int[numberOfAreas];
            for (int i=0; i<numberOfAreas; i++)
            {
                centerIds[i] = listAreas[i].centerId;
                caseIDsInAreas[i] = listAreas[i].caseIDs;
                caseCounts[i] = listAreas[i].numCases;
            }
            
            for(int i = 0; i < numberGene; i++)
                bestWeights[i] = weightsCalculated[i];

            calculateWeights(numberOfAreas, caseIDsInAreas, caseCounts,numberOfMaximumDimensions, centerIds, averagedCenters, lpType, noConstraints);
            delete[] centerIds;
            delete[] caseIDsInAreas;
            delete[] caseCounts;

            for (int count2=0; count2<numberGene; count2++)
            {
                tempWeights[count2]=weightsCalculated[count2];
                if (tempWeights[count2]>0)
                    printf("Temp Weights %lf\n", tempWeights[count2]);
            }

            qsort(tempWeights, numberGene, sizeof(double), compare);
            previousWeight = currentWeight;
            currentWeight = tempWeights[numberGene-numberOfMaximumDimensions];
            for (int count2=0; count2<numberGene; count2++)
            {
                
                if (weightsCalculated[count2]<=tempWeights[numberGene- numberOfMaximumDimensions ])
                {
                    weightsCalculated[count2]=0;
                }
                else
                    printf("Non Zero Weights %lf %lf\n", weightsCalculated[count2], tempWeights[numberGene- numberOfMaximumDimensions]);
            }
        }
        countX++;
        printf("\nIteration Count: %d\t", countX);
        if(cont)
            printf("continuing\n");
        else printf("finished\n");
    }
    fprintf(fpOut, "%i\n",bestNumberOfAreas);

    for (int count=0; count<numberGene; count++)
    {
        fprintf(fpOut,"%lf ", bestWeights[count]);
    }

    for(int numCenters = 0; numCenters<bestNumberOfAreas; numCenters++)
    {
        printf("\nArea %d Size: %d\n", numCenters, bestAreas[numCenters].numCases);
        
        for (int count2=0; count2<bestAreas[numCenters].numCases; count2++)
            printf("Samples Being Covered %s\n", listSamples[bestAreas[numCenters].caseIDs[count2]].sampleName);

        fprintf(fpOut,"\n%d\n", bestAreas[numCenters].numCases);
        for (int count2=0; count2<bestAreas[numCenters].numCases; count2++)
            fprintf(fpOut,"%s ", listSamples[bestAreas[numCenters].caseIDs[count2]].sampleName);
    }

    updateDistanceMatrix(bestWeights);
    double* oneLastFinal = getRadiuses();
    double* radisusControls = new double[totalSamples];


    int totalPredictions = 0, trueCase = 0, trueControl = 0, falseCase = 0, falseControl = 0; 
    std::vector<int>* neighborPts;

    int* nbCounts = getCaseNeighbors(oneLastFinal, neighborPts, epsilonIsMinRadius);

    for(int i = 0; i < totalSamples; i++)
        if(listSamples[i].trainOrTest==1 && listSamples[i].filterSample==0){
            totalPredictions++;
            double radius = 0;
            if(epsilonIsMinRadius){
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
                    
                    if(epsilonIsMinRadius)
                        gotcha = std::min(oneLastFinal[bestAreas[j].caseIDs[k]], radius) > 
                            distanceMatrix[i][bestAreas[j].caseIDs[k]];
                    else gotcha = oneLastFinal[bestAreas[j].caseIDs[k]] > distanceMatrix[i][bestAreas[j].caseIDs[k]];

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
        delete[] neighborPts;
        delete[] nbCounts;
        delete[] radisusControls;
        delete[] oneLastFinal;

        printf("Out of total %d predictions:\n", totalPredictions);
        printf("%d is guessed to be case, %d was control\n",  falseCase+trueCase, falseCase);
        printf("%d is guessed to be control, %d was case\n", falseControl+trueControl, falseControl);
        printf("false prediction: %f", falseCase/((double)falseCase+trueCase));

    fclose(fpOut);
}




