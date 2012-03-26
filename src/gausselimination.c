/*subroutines for Gaussian Elimination solver */

#include "ddamemory.h"

void forwardSubstitution(double ** a, int n) {
	int i, j, k, max;
	int notZerosCount;
	double t;
	for (i = 0; i < n; ++i) {
		notZerosCount = 0;
		/*deal with first row of all zeros, if present*/
		for (j = 0; j<= n; j++){
			if(a[i][j] != 0)
			notZerosCount++;
		}

		/* if all zeros, skip that row */
		if (notZerosCount == 0)
			continue;


		max = i;
		for (j = i + 1; j < n; ++j)
			if (a[j][i] > a[max][i])
				max = j;
		
		for (j = 0; j < n + 1; ++j) {
			t = a[max][j];
			a[max][j] = a[i][j];
			a[i][j] = t;
		}
		
		for (j = n; j >= i; --j)
			for (k = i + 1; k < n; ++k)
				a[k][j] -= a[k][i]/a[i][i] * a[i][j];

/*		for (k = 0; k < n; ++k) {
			for (j = 0; j < n + 1; ++j)
				printf("%.2f\t", a[k][j]);
			printf("\n");
		}*/
	}
}

void reverseElimination(double ** a, int n) {
	int i, j, notZerosCount;
	for (i = n - 1; i >= 0; --i) {
		notZerosCount = 0;
		/*deal with first row of all zeros, if present*/
		for (j = 0; j<= n; j++){
			if(a[i][j] != 0)
			notZerosCount++;
		}

		/* if all zeros, skip that row */
		if (notZerosCount == 0)
			continue;

		a[i][n] = a[i][n] / a[i][i];
		a[i][i] = 1;
		for (j = i - 1; j >= 0; --j) {
			a[j][n] -= a[j][i] * a[i][n];
			a[j][i] = 0;
		}
	}
}

void gauss(double ** A, double * b, int sizeAMatrix, double * returnVec) {

	double ** combinedMatrix;
	double * resultVec;
	int i,j;

	combinedMatrix = DoubMat2DGetMem(sizeAMatrix, sizeAMatrix+1);
	
	/*form a combined matrix from A and B*/
	for(j = 0; j < (sizeAMatrix+1); j++){
		for(i = 0; i < (sizeAMatrix); i++){
			if (j!= sizeAMatrix)
				combinedMatrix[i][j] = A[i][j];
			else
				combinedMatrix[i][j] = b[i];
		}/*i*/
	}/*j*/	

	forwardSubstitution(combinedMatrix, sizeAMatrix);
	reverseElimination(combinedMatrix, sizeAMatrix);

	//resultVec = (double *)malloc(sizeAMatrix*sizeof(double));
	//for (i=0; i<sizeAMatrix; i++)
	//	resultVec[i] = combinedMatrix[i][sizeAMatrix];
	
	for (i=0; i < sizeAMatrix; i++){
	//	returnVec[i] = resultVec[i-1];
		returnVec[i] = combinedMatrix[i][sizeAMatrix];
	}

	//for (i = 0; i < sizeAMatrix; i++)
	//	free(combinedMatrix[i]);
	//free(combinedMatrix);

	free2DMat(combinedMatrix,sizeAMatrix);

	//free(resultVec);

}