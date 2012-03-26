
/* There should be no windows dependencies in this file. */

#ifndef __FLUIDSDATA_H__
#define __FLUIDSDATA_H__

#ifdef WIN32
#pragma warning( disable : 4115 )        
#endif

#include "ddatypes.h"
#include "geometrydata.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif


/*--------------Fluidsdata--------------*/

typedef struct _fluidsdata_tag Fluidsdata;

struct _fluidsdata_tag {

	int nFluidJoints;
	int jointFluidPropsize1;
	int jointFluidPropsize2;
	double ** jointFluidProp;
	/* width friction_factor*/

	int nPSegments;
	int pipesize1;
	int pipesize2;
	double ** pipeSegments;
	/* 0 v1a v1b v2a v2b meanHydraulicAperture friction_factor resistance(k) node1 node2 */

	int connsize1;
	int connsize2;
	int ** connectivityTable;
	/* 0 startnode endnode*/
	
	int nNodes;
	int nodesize1;
	int nodesize2;
	double ** nodes;
	/* # of pipes connected to that node, x coord, y coord, index of up to 10 pipes*/

	int kmatsize1;
	int kmatsize2;
	double ** kMatrix;
	/* i,j contains the k values (ie, element resitivity) between node i and node j */

	int nBoundaryCond;
	int boundsize1;
	int boundsize2;
	double ** boundaryConditions;
	/* Node, XCoord, YCoord, Condition, Flag 1 = Flow or Flag 2 = Head */

	int cvectorsize1;
	double * cVector;

	int headssize1;
	double * heads;

	int pressuressize1;
	double * pressures;

	int mpNodeIndexsize1;
	int * mpNodeIndex;

	int vertNodeIndexsize1;
	int * vertNodeIndex;

	double fluidDensity;
	double fluidViscosity;
	double fluidMu;

	int ksolverflag;
	/* 1 for laminar, 2 for cubic*/

	int gravityflag;
	/* 0 if a horizontal problem, 1 if vertical */

	void (*free)(Fluidsdata *);
 
};


void           fdata_read_input_file     (Fluidsdata *, 
                                          char * filename);


/* void           fdata_write_ddaml          (Fluidsdata * fd, 
                                           PrintFunc printer, 
                                           void * stream); */

Fluidsdata * fdata_new(void);

void fdata_delete(Fluidsdata *);

void initFluidMatrices(Fluidsdata *fd, Geometrydata *gd);


/*--------------Boundarycond----------------*/

typedef struct _boundarycond BoundaryCond;

struct _boundarycond {
  /* Need a union of structs and an 
   * enum in here to handle the various types
   * of friction laws.
   */
   double xCoord;
   double yCoord;
   double condition;
   int type;
};

BoundaryCond * BoundaryCond_new (void);

void boundarycond_set_xcoord (BoundaryCond *bc, double xcoord);
void boundarycond_set_ycoord (BoundaryCond *bc, double ycoord); 
void boundarycond_set_condition (BoundaryCond *bc, double condition);
void boundarycond_set_type (BoundaryCond *bc, int type);

double boundarycond_get_xcoord (BoundaryCond *bc);
double boundarycond_get_ycoord (BoundaryCond *bc); 
double boundarycond_get_condition (BoundaryCond *bc);
int boundarycond_get_type (BoundaryCond *bc);

/*Other functions */


#ifdef __cplusplus
}
#endif


#endif /* __FLUIDSDATA_H__ */


