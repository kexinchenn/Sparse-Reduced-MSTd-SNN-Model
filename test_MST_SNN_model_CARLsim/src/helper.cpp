#include <assert.h>
#include <cmath>
#include <vector>
#include <stdio.h>

#include "helper.h"
#include "util.h"

using namespace std;

void reconstructMT(vector<vector<float> > weights, vector<float> FRs, float* recMT) {
    int wNumRow = weights.size();
    int wNumCol = weights[0].size();
    int numNeur = FRs.size();

    assert (wNumCol == numNeur);

	float productSum;

    for (unsigned int i = 0; i < wNumRow; i ++) {
    	productSum = 0.0;
    	for (unsigned int j = 0; j < wNumCol; j ++) {
    		if (isnan(weights[i][j]) || weights[i][j] < 0) {
    			weights[i][j] = 0;
    		}
    		productSum += weights[i][j] * FRs[j];
	    }
	    recMT[i] = productSum;
    }
}  

