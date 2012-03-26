/*File contains all of the subroutines used to calculate flow through the pipe network
Actual pipe and node creation occurs through functions in pipes.c*/

/*Single Phase Network flow, follows the algorithm developed in Brebbia and Ferrante 1983 */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

#include <math.h>

#include "dda.h"
#include "ddamemory.h"

void calckCoeff(Fluidsdata *fd, Geometrydata *gd, Analysisdata *ad);
void generateKMatrix(Fluidsdata *fd);
void assignBoundaryConditionsToNodes(Fluidsdata *fd);
void generateCVector(Fluidsdata *fd);
//int Gaussian_Elimination(double **A, int n, double *B);  
void gauss(double ** A, double * b, int sizeAMatrix, double * returnVec);
void calcHeads(Fluidsdata * fd, Analysisdata *ad);



void calcPressures(Fluidsdata * fd, Geometrydata * gd, Analysisdata *ad){

	/*follows the algorithm in Brebbia and Ferrante 1983*/

	calckCoeff(fd,gd,ad);

	generateKMatrix(fd);

	assignBoundaryConditionsToNodes(fd);

	generateCVector(fd);

	calcHeads(fd, ad);


}/*calcPressures*/

void calckCoeff(Fluidsdata * fd, Geometrydata * gd, Analysisdata *ad){
	
	int pipe_i;
	double k, width, length, g;
	double x1, y1, x2, y2;
	double friction_factor, rho, mu;

	int ** connectivityTable;
	double ** nodes;

	connectivityTable = fd->connectivityTable;
	nodes = fd->nodes;
	

	/*for now, using the book method, for testing */
	/*k coefficients are stored in the first column (index 0) of the connectivity matrix
	Currently, k is equal to the constant developed from using Poiseuille's Law*/

	/*both Jing 2001 and Ben 2011 just solve for linear flow*/

	for (pipe_i = 1; pipe_i <= fd->nPSegments; pipe_i++){
		
		x1 = nodes[connectivityTable[pipe_i][1]][1];
		y1 = nodes[connectivityTable[pipe_i][1]][2];
		x2 = nodes[connectivityTable[pipe_i][2]][1];
		y2 = nodes[connectivityTable[pipe_i][2]][2];

		length = sqrt(pow((x2-x1),2)+pow((y2-y1),2));
		width = fd->pipeSegments[pipe_i][5];
		friction_factor = fd->pipeSegments[pipe_i][6];
		g = ad->gravaccel;
		rho = fd->fluidDensity;
		mu = fd->fluidViscosity;
		
		if(fd->ksolverflag == 1){
			/*laminar flow, pg 72 Brebbia and Ferrante*/
			k = PI*rho*g*pow(width,4)/(128*mu*length);
			fd->pipeSegments[pipe_i][7] = k;	
		}

		if(fd->ksolverflag == 2){
			/*cubic flow, Jing 2001 */
			k = rho*g*pow(width,3)/(12*mu*length);
			fd->pipeSegments[pipe_i][7] = k;	
		}
	
	}/*pipe_i*/

}/* calckCoeff */

void generateKMatrix(Fluidsdata *fd){
	
	int i, j, pipe_i;
	double k;

	/*
	fd->kmatsize1 = fd->kmatsize2 = fd->nNodes+1;
	/*generate matrix of ks*/
    /*
	fd->kMatrix = DoubMat2DGetMem(fd->kmatsize1, fd->kmatsize2);
	*/

	/*fd->kMatrix has to be reset to zeros*/
	for(i = 0; i<fd->kmatsize1; i++){
		for(j = 0; j<fd->kmatsize2; j++){
			fd->kMatrix[i][j] = 0.00;
		}
	}	

	/*for each pipe, add k to positions ii and jj, and subtract it from ji and ij,
	where i and j are the assumed start and end nodes, respectively*/

	for (pipe_i = 1; pipe_i <= fd->nPSegments; pipe_i++){

		i = fd->connectivityTable[pipe_i][1];
		j = fd->connectivityTable[pipe_i][2];
		k = fd->pipeSegments[pipe_i][7];

		fd->kMatrix[i][i] += k;
		fd->kMatrix[j][j] += k;
		fd->kMatrix[j][i] -= k;
		fd->kMatrix[i][j] -= k;

	}/* pipe_i */
} /*generateKMatrix */

void assignBoundaryConditionsToNodes (Fluidsdata *fd){
	
	int bound_i, node_i;
	double xb, yb, xn, yn;

	/*assigns a node index to each boundary condition based on its xy coord*/
	for(bound_i = 1; bound_i <= fd->nBoundaryCond; bound_i++){
		
		xb = fd->boundaryConditions[bound_i][1];
		yb = fd->boundaryConditions[bound_i][2];

		for(node_i = 1; node_i <= fd->nNodes; node_i++){
			
			xn = fd->nodes[node_i][1];
			yn = fd->nodes[node_i][2];

			/*compare xs and ys */
			if(xb == xn && yb == yn){
				fd->boundaryConditions[bound_i][0] = node_i;	
				break;
			}/*if*/

		}/*node_i*/
	}/*bound_i*/

}/*assignBoundaryConditionsToNodes*/
		
void generateCVector(Fluidsdata *fd){

	int cvecsize1, bound_i, node_ind, k_i,k_j;
	double k_orig, bound_head;


	cvecsize1 = fd->cvectorsize1;

	/*sets all cs equal to zero*/
	for (node_ind = 1; node_ind<=cvecsize1; node_ind++){
		fd->cVector[node_ind] = 0.0;
	}			

	/*currently, only either flow or head (not both) boundaries 
	may be specified, based on the allowed inputs in the fluids dtd*/

	for (bound_i = 1; bound_i <= fd->nBoundaryCond; bound_i++){

		/*for flow, add the flowrates at node i to Ci for each flow boundary condition*/
		if(fd->boundaryConditions[bound_i][4] == 1){
			
			node_ind = (int)fd->boundaryConditions[bound_i][0];
			fd->cVector[node_ind] += fd->boundaryConditions[bound_i][3];

		}/*if*/
		else if(fd->boundaryConditions[bound_i][4] == 2){
			
			node_ind = (int)fd->boundaryConditions[bound_i][0];
			k_orig = fd->kMatrix[node_ind][node_ind];
			bound_head = fd->boundaryConditions[bound_i][3];

			/*for head, first subtract k[k_i][node_ind]*bound_head from cvec[k_i]
			then set all k's related to that node to 0,
			and k[node_ind][node_ind] to 1*/
			for (k_i = 1; k_i <= fd->nNodes; k_i++){
				fd->cVector[k_i] -= fd->kMatrix[k_i][node_ind]*bound_head;
				fd->kMatrix[k_i][node_ind] = 0;
			}/*k_i*/
			for (k_j = 1; k_j <= fd->nNodes; k_j++){
				fd->kMatrix[node_ind][k_j] = 0;
			}/*k_j*/
			fd->kMatrix[node_ind][node_ind] = 1;

			/*finally, fd->cVector[node_ind] is set to the boundary head, 
			and it should not be changed again over the course of the 
			algorithm*/
			fd->cVector[node_ind] = bound_head;
			
		}/*if*/

	}/*bound_i*/

}/*generateCVector*/

void calcHeads(Fluidsdata *fd, Analysisdata *ad){

	/*Note, at least one head needs to be specified, or all of the heads
	will be off by some random constant of the same magnitude */

	int node_i;
	double gamma;

	gamma = fd->fluidDensity*ad->gravaccel; /*specific gravity*/

	/*write new c's to Heads*/

	/*currently using basic Gaussian Elimination */
	gauss(fd->kMatrix, fd->cVector, fd->kmatsize1, fd->heads);

	if (fd->gravityflag == 1){
		/*calculate pressures from heads (subtracting the elevation head)*/
		for(node_i = 1; node_i<=fd->nNodes; node_i++){
			fd->pressures[node_i] = (fd->heads[node_i]-fd->nodes[node_i][2])*gamma;	
		}
	} 
	else if (fd->gravityflag == 0){
		/*calculate pressures from heads (ignoring the elevation head)*/
		for(node_i = 1; node_i<=fd->nNodes; node_i++){
			fd->pressures[node_i] = (fd->heads[node_i])*gamma;	
		}
	}

}/* calcHeads */




