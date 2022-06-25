#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>

#include "util.h"

using namespace std;

float calcMean(float *in, int size) {
    float mean = 0.0;
    float sum = 0.0;

    for (int i = 0; i < size; ++i)
    {
        sum += in[i];
    }
    mean = sum / (float)size;

    return mean;
}

float calcRMS(float** x,float** y, int numRow, int numCol) {
    float res = 0;
    float denom;
    float numerator;

    // sum of squres
    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            res += pow((x[i][j] - y[i][j]), 2.0);
        }
    }
	
    numerator = sqrt(res);
    denom = sqrt(numRow * numCol);
    res = numerator / denom;
    return res;
}

float* findMaxOfEachColumn(float **in, int numRow, int numCol) {
    float* maxVals = new float[numCol];
    int currMax;

    for (int j = 0; j < numCol; j++) {
        currMax = 0;
        for (int i = 0; i < numRow; i++) {
            if (in[i][j] > currMax) {
                currMax = in[i][j];
            }
        }
        maxVals[j] = currMax;
    }
    return maxVals;
}


float** loadData(string file, int numRow, int numCol) {

    float** out_data;
    out_data = new float*[numRow];
    for (int i = 0; i < numRow; i++) {
        out_data[i] = new float[numCol];
    }

    ofstream fileErrors;
    fileErrors.open("./results/data_file_errors.txt", ofstream::out | ofstream::app);
    
    ifstream ip;
    ip.open(file.c_str());
    
    
    if (ip.fail()) {
        fileErrors << file << " failed" << endl;
    }

    string line, field;
    
    while (!ip.eof()) {
        for (int i = 0; i < numRow; i++) {
            // out_data[i] = new float[numCol];

            getline(ip, line);
            istringstream s(line);

            for (int j = 0; j < numCol; j++) {
                getline(s, field, ',');
                istringstream str(field);
                str >> out_data[i][j];
            }
        }
    }

    ip.close();
    fileErrors.close();
    return out_data;
}

int randGenerator (int i) {return rand() % i;}

void shuffleTrials(int numTrials, int numTrain, int numTest, int *trainTrials, int *testTrials) {
    vector<int> allTrials;
    for (unsigned int i = 0; i < numTrials; i++) {
        allTrials.push_back(i);
    }

    random_shuffle(allTrials.begin(), allTrials.end(), randGenerator);
    
    for (unsigned int i = 0; i < numTrain; i++) {
        trainTrials[i] = allTrials[i];
    }
    for (unsigned int i = 0; i < numTest; i++) {
        testTrials[i] = allTrials[i+numTrain];
    }
}

void calcCorrCoef(float** X, float** Y, int numRow, int numCol, float* corrCoef) {
    float sub=-1.0e-20;// number to avoide division with zero
    float meanX, meanY, sumXX, sumYY, sumXY;

    for (unsigned int i = 0; i < numCol; i ++) {
        meanX = meanY = 0.0;
        for (unsigned int j = 0; j < numRow; j ++) {
            meanX += X[j][i];
            meanY += Y[j][i];
        }
        meanX /= numRow;
        meanY /= numRow;

        sumXX = sumYY = sumXY = 0.0;
        for (unsigned int j = 0; j < numRow; j ++) {
            sumXX += (X[j][i] - meanX) * (X[j][i] - meanX);
            sumYY += (Y[j][i] - meanY) * (Y[j][i] - meanY);
            sumXY += (X[j][i] - meanX) * (Y[j][i] - meanY);
        }
        corrCoef[i] = sumXY / (sqrt(sumXX * sumYY) + sub);
    }
}

float calcCorr(float** X, float** Y, int numRow, int numCol) {
    /*  correlation of two 2D matrices */
    float corrCoef;
    int length = numRow * numCol;
    float* X_flat = new float [length];
    float* Y_flat = new float [length];

    // convert input matrices X and Y to one dimensional arrays
    int count = 0;
    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            X_flat[count] = X[i][j];
            Y_flat[count] = Y[i][j];
            count ++;
        }
    }
    // compute correlation between two arrays
    corrCoef = calcPopCorrCoef(X_flat, Y_flat, length);

    return corrCoef;
}


float calcPopCorrCoef(float* x,float* y, int length)
{
    const float TINY=1.0e-20;//tiny number to avoide division with zero
    float res = 0;
    int j;
    float syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

    for (j=0;j<length;j++) 
        {
        ax += x[j];
        ay += y[j];
    }
    ax /= length;
    ay /= length;
    for (j=0;j<length;j++) 
        {     
        sxx += (x[j]-ax)*(x[j]-ax);
        syy += (y[j]-ay)*(y[j]-ay);
        sxy += (x[j]-ax)*(y[j]-ay);
    }
    res=sxy/(sqrt(sxx*syy)+TINY);
        return res;
}

