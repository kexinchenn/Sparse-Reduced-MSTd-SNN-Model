#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>

#include "util.h"

using namespace std;

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

