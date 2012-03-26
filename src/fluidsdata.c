/*Initital template for this file copied over from analysisdata.c*/

#ifdef WIN32
#pragma warning( disable : 4115 )        
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>


#include "dda.h"
#include "ddamemory.h"
#include "ddafile.h"
#include "ddaml.h"

#define I1 "   "
#define I2 "      "



void initializeFData(DDA * dda, char * ffilename){
			
	   Fluidsdata * fd = fdata_new();
	   
	   ddaml_read_fluids_file(fd,ffilename);

	   dda_set_fluidsdata(dda,fd);
};

void initFluidMatrices(Fluidsdata *fd, Geometrydata *gd){

	/*all of these were pulled out of fluidsdriver.c*/
	fd->kmatsize1 = fd->kmatsize2 = fd->nNodes+1;
	fd->kMatrix = DoubMat2DGetMem(fd->kmatsize1, fd->kmatsize2);
	
	fd->cvectorsize1 = fd->nNodes+1;
	fd->cVector = (double *)calloc(fd->cvectorsize1, sizeof(double));

};

void 
fdata_delete(Fluidsdata * fd) {

   free2DMat((void **)fd->jointFluidProp, fd->jointFluidPropsize1);
   free2DMat((void **)fd->boundaryConditions, fd->boundsize1);
   free2DMat((void **)fd->connectivityTable, fd->connsize1);
   free2DMat((void **)fd->kMatrix, fd->kmatsize1);
   free2DMat((void **)fd->nodes, fd->nodesize1);
   free2DMat((void **)fd->pipeSegments, fd->pipesize1);
   
   free(fd->heads);
   free(fd->pressures);
   free(fd->cVector);

   free(fd);

}

Fluidsdata *
fdata_new() {
 
   Fluidsdata * fdo;

  /* FIXME:  change to malloc. */
   fdo = (Fluidsdata *)calloc(1,sizeof(Fluidsdata));
   
/**************************/

  /* FIXME: Set up the private functions.  These can be moved in the
   * future to the analysis initialization function.
   */

   fdo->free = fdata_delete; //freeFluidsData;

   return fdo;
}  

void       
boundarycond_set_xcoord (BoundaryCond *bc, double xcoord) {
   bc->xCoord = xcoord;
}
void       
boundarycond_set_ycoord (BoundaryCond *bc, double ycoord) {
   bc->yCoord = ycoord;
}
void       
boundarycond_set_condition (BoundaryCond *bc, double condition) {
   bc->condition = condition;
}
void       
boundarycond_set_type (BoundaryCond *bc, int type) {
   bc->type = type;
}

double     
boundarycond_get_xcoord(BoundaryCond * bc) {
   return bc->xCoord;
}
double     
boundarycond_get_ycoord(BoundaryCond * bc) {
   return bc->yCoord;
}
double     
boundarycond_get_condition(BoundaryCond * bc) {
   return bc->condition;
}
int     
boundarycond_get_type(BoundaryCond * bc) {
   return bc->type;
}

BoundaryCond * 
BoundaryCond_new (void) {

   BoundaryCond * bc = (BoundaryCond*)malloc(sizeof(BoundaryCond));
   memset(bc,0x0,sizeof(BoundaryCond));
   return bc;
}

