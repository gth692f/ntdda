
/*
 *  Add some comments here!
 */

#include <math.h>
#include <assert.h>

#include "analysis.h"
#include "ddamemory.h"
#include "ddadlist.h"
#include "datalog.h"
#include "bolt.h"
#include "timehistory.h"

#include "stress.h"



/* Helper var for examining variables */
static char mess[80];

/** Move this into its own file or the timehistory file. */
static void seismicload(DList * seispoints,TimeHistory * timehistory, int timestep,
                        int *k1, double **F, double ** moments, 
                        double ** matprops, TransMap transmap);









/**************************************************/
/* df10: initiate coefficient and load matrix     */
/**************************************************/
/*  I think this means initialize rather than initiate.
 *  Either way this routine needs to be abstracted into 
 *  "set_matrix_to_zero" function, and the time step 
 *  stuff put somewhere else.
 */
void df10(Geometrydata *bd, Analysisdata *ad, int **locks, 
          double **F) {
   int i;
   int nBlocks = bd->nBlocks;
   int nContacts = bd->nContacts;  /* Not changed in df10(). */

   double ** globalTime = ad->globalTime;
   double ** K = ad->K;

  /* initial a[][]=0 f[][]=0  */
  /* These two calls get moved out of here eventually
   * and the matrices initialized directly before being
   * used in the analysis.  This function needs to be 
   * turned into a macro.  FIXME: These need to get turned in 
   * 1D arrays, then a simple call to memset.
   */
   setMatrixToZero(K, ad->ksize1, ad->ksize2);
   setMatrixToZero(F, nBlocks+1, 7);

  /* (GHS: for reduce time interval     m0[][1]=0 in df07) */
  /* This will have to stay here as the whole matrix is
   * not initialized, only one index. 
   * FIXME: GHS indicates this may be needed only if the 
   * time interval is reduced.
   */
   for (i=1; i<= nContacts; i++)
   {
      locks[i][1]=0;  // was m0
   }  /*  i */
 
   globalTime[ad->cts][0] = globalTime[ad->cts-1][0] + ad->delta_t;
   globalTime[ad->cts][2] = ad->delta_t;

}  




/**************************************************/
/* df12: submatrix of fixed points                */
/**************************************************/
void df12(Geometrydata *gd, Analysisdata *ad, int *k1, 
          double **F, double **blockArea, int **n,
          TransMap transmap)
{
   int i;
  /* i0 is the block number that the point is
   * associated with.
   */
   int blocknumber;  // was i0
   //int i1;  // redundant block number
  /* Memory locations */
   int i2, i3;
   int j;
   int j1;
   int ell;
   int nFPoints = gd->nFPoints;
   double ** points = gd->points;
   double e[7][7];
   double s[31];
   double ** K = ad->K;
  /* After the first time step, c seems to carry point displacement
   * from the previous time step.  This is used to compute a 
   * ``restoring force'' opposite to the direction that the 
   * fixed point moved during the last time step.
   */
   double ** c = ad->c;
   double x;
   double y;
   double T[7][7];

  /* g4 is a spring value set at the beginning from the 
   * initial contact spring stiffness.  The initial contact
   * spring stiffness may be changed adaptively during run
   * time, so we can't use that.  g4 has been changed to 
   * FPointSpring.  See GHS 1988, Chapter 2, p. 93.
   * FIXME:  Initialize FPointSpring in a more appropriate
   * place.
   */
   double p = (100.0)*ad->FPointSpring;


  /* i0: block number  */
   for (i=1; i<= nFPoints; i++) {


     /* i0 is the block number that the point is 
      * associated with.  Notice that the fixed points
      * are read in as a line so that the 3,4 slots
      * are an endpoint.  This is overwritten in the 
      * geometry pre-processor.
      */
      blocknumber = (int)points[i][3];
      x=points[i][1];
      y=points[i][2];

      transmap(blockArea[blocknumber],T,x,y);
      
      /* (GHS: 6*6 submatrix of fixed point for a[][]) */
      for (j=1; j<= 6; j++) {
         for (ell=1; ell<= 6; ell++) {

            e[j][ell] = p*(T[1][j]*T[1][ell] + T[2][j]*T[2][ell]);
         }
      } 

	    /* K[][] here is only for fixed points: Eq. 2.78,
      * p. 94, Chapter 2, Shi 1988.
      */  
      //i1=(int)points[i][3];
		//i2=k1[i1];
      i2 = k1[blocknumber];
		i3 = n[i2][1]+n[i2][2]-1;

      for (j=1; j<= 6; j++) {
         for (ell=1; ell<= 6; ell++) {

            j1=6*(j-1) + ell;
            K[i3][j1] += e[j][ell];
         }
      }
      
	  
     /* (GHS: 6*1 submatrix of fixed displacements for f[]) */
     /* points[][4] and points[][5] are most likely (x,y)
      * load values as set in df09() from the timeDeps array.
      * The matrix c stores displacement values of the fixed
      * points.  The penalty method does not eliminate displacement
      * of fixed points, so a ``restoring force'' oriented 
      * opposite to the displacement is added to the force vector
      * for the next (current) time step.  NOTE: For this to work 
      * correctly, the fixed points must at the beginning of the 
      * points array.
      */
      for (j=1; j<= 6; j++) {

         s[j] = (T[1][j]*c[i][1]) + (T[2][j]*c[i][2]);
      }
      
      
     /* (GHS: for constrained or fixed point)  */
     /* This next block of code is not derived for 
      * fixed points (and constrained derivation 
      * was noted as unclear).  The idea here is that 
      * motion of fixed points is compensated by pushing
      * in a direction opposite of that which the fixed
      * point displaced.  The fixed point displacement 
      * is computed and stored in df25().
      */
      for (j=1; j<= 6; j++) {

         F[i2][j] += p*s[j];
      }

   }

}  /* Close  df12().  */





/**************************************************/
/* df15: submatrix of point loading               */
/**************************************************/
void df15(Geometrydata *gd, Analysisdata *ad, int *k1, double **F, 
          double **blockArea, TransMap transmap) {


   int i, i0, i1, i2;
   int j;
   int nLPoints = gd->nLPoints;
   int nFPoints = gd->nFPoints;
   double ** points = gd->points;
   double s[7];
   double ** a = ad->K;
   double x, y;
   double T[7][7];
   double xload;
   double yload;

   for (i=nFPoints+1; i<= nFPoints+nLPoints; i++) {

     /* i0 is probably the number of the block 
      * associated with a particular load point.
      */
      i0=(int)points[i][3];
     /* x and y are the coordinates where the point 
      * load is applied.
      */
      x=points[i][1];
      y=points[i][2];

      xload = points[i][4];
	  yload = points[i][5];
      transmap(blockArea[0],T,x,y);

      y=points[i][2];

      xload = points[i][4];
      yload = points[i][5];

      transmap(blockArea[i0],T,x,y);

   
     /* points[i][4] points[i][5] 
      * is time dependent loads initialized 
      * df09()  
      */
      for (j=1; j<= 6; j++) {
         s[j] = T[1][j]*points[i][4] + T[2][j]*points[i][5];
      }  

      i1=(int)points[i][3];
      i2=k1[i1];
      
     /* BUG: There is a bug in here that results from 
      * having a point load located outside of the
      * rock mass.
      */
      for (j=1; j<= 6; j++) {
         F[i2][j] += s[j];
      }  
   }  

}  


void
seismicload(DList * seispoints, TimeHistory * timehistory, int timestep,
            int *k1, double ** F, double ** moments, 
            double ** e0, TransMap transmap)
{

   DList * ptr;
   DDAPoint * ptmp;
   int j;
   int i2; // memory location for block in forcing vector
   int i0;
   double x,y;
   double force, mass, accel;
   double T[7][7] = {0.0};
   double s[7] = {0.0};

   if (timestep >= th_get_number_of_datapoints(timehistory)) {
      return;
   }

   dlist_traverse(ptr, seispoints) {

      ptmp = (DDAPoint *)ptr->val;
      i0 = ptmp->blocknum;
      x = ptmp->x;
      y = ptmp->y;
      transmap(moments[i0],T,x,y);

      accel = th_get_data_value(timehistory, timestep);

      mass = moments[i0][1]*e0[i0][1];
      force = accel*mass;

      for (j=1; j<= 6; j++) {
         s[j] = T[1][j]*force + T[2][j]*0.0;
      }  

      i2 = k1[i0];

     /* BUG: There is a possible bug in here that results from 
      * having a point load located outside of the
      * rock mass.
      */
      for (j=1; j<= 6; j++) {
         F[i2][j] += s[j];
      }  
   }

}  



/**************************************************/
/* df16: submatrix of volume force                */
/**************************************************/
/*  k1 is ???
 *  F is the "force" vector.
 *  e0 is material props: ma we e0 u0 c11 c22 c12 t-weight
 */
void df16(Geometrydata *gd, int *k1, double **F, double **e0,
          double **moments, Fluidsdata* fd, Analysisdata * ad)
{
   int block;
   int i2;
   double o1;
   double o2;
   int nBlocks = gd->nBlocks;
   int gravityflag = 1;

   /*WEM: currently, option for turning off gravity only applicable
   with pipe flow.  Shouldn't be too hard to change, but this 
   is the way it was in the original program*/

   if (ad->pipeflowflag == 1)
	   gravityflag = fd->gravityflag;

   if (gravityflag == 1){
	   for (block=1; block<= nBlocks; block++) {

		  i2 = k1[block];
		  o1 = 0;
		 /* e0[1] needs to be unit weight. 
		  * FIXME: unit weight is used here and 
		  * density is used in inertia matrix.
		  * Problem is that gravity may not be consistent
		  * between the two.  The thing to do here is to
		  * use the density and the gravity acceleration 
		  * value given in the input file to turn gravity
		  * loading on and off.
		  * FIXME: Need to account for changing density here as 
		  * well, either by thickness or by density.
		  */
		  o2 = -e0[block][1];
		  F[i2][1] += o1*moments[block][1];
		  F[i2][2] += o2*moments[block][1];
	   }  
   }
} 

void 
fluid_pressure_load(Fluidsdata *fd, Geometrydata *gd, 
					Analysisdata * ad, double **blockArea,
					int * k1, TransMap transmap){
	
	/*Initial Algorithm for fluid flow taken from 
	Kim, Amadei and Pan (KAP) (1999),
	may need to be upgraded later on to the Jing 2001 paper. */

	/*Implementation of block loads begins on page 957 of KAP (1999)*/

    int v_b_start, v_b_end;
	int vp1, vp2; /*block and pipe vertex indices*/
	int block_i, pipe_i, force_i;
	int block_1or2, node1, node2;
	int i2;
	int vt1, vt2, vt3, vt4; /*temporary pipe vertex indices*/

	double L, d, p1, p2, F1, F2, F, Fx, Fy; 
	/*variable names are those assigned in KAP 1999*/
	double x1, y1, x2, y2, xf, yf, alpha;

	double n = 1; /*contact porosity, for now set to a unit value*/
	double T[7][7] = {0.0};


	/*algorithm works by applying the pressure from each pipe
	to the corresponding block.*/
	for (pipe_i = 1; pipe_i<=fd->nPSegments; pipe_i++){

		/*if width of the pipe equals zero, no water is in it,
		and thus no fluid pressure.  continue skips this pipe,
		assigning no pressure and moves on to the next one
		Shouldn't need a tolerance, since these values never change*/
		if (fd->pipeSegments[pipe_i][5] == 0)
			continue;

		/*Can't remember why I set it to if the width is zero, as opposed
		to if the flag in the zero index of pipes is set to one...
		Alternate between these for debugging... Reason for one or the other,
		apparent in the SquarePipe2 unit test*/

		//if (fd->pipeSegments[pipe_i][0] == 1)
		//	continue;
	
		/*each pipe applies pressure to two adjacent blocks
		(unless it's a boundary pipe), in which case pressure 
		applies to only one block*/

		for (block_1or2 = 1; block_1or2 <= 2; block_1or2++){
			/*using the first block's endpoints (index 1 and 3),
			then repeat using the second (index 2 and 4)*/

			vt1 = (int)fd->pipeSegments[pipe_i][1];
			vt2 = (int)fd->pipeSegments[pipe_i][2];	
			vt3 = (int)fd->pipeSegments[pipe_i][3];
			vt4 = (int)fd->pipeSegments[pipe_i][4];;

			/*if boundary node*/
			/*this should make the code only pass through
			the block_1or2 for loop once; Without this, it loops
			twice on the boundary nodes, no bueno*/

			if(vt1 == vt2 && vt3 == vt4){
				block_1or2 = 2;
			}

			if (block_1or2 == 1){
				vp1 = vt1;
				vp2 = vt3;
			} 
			else if (block_1or2 == 2) {
				/* indexing based on counterclockwise rotation
				on block vertices, allows easier formulation
				for calculating the force directions */
				vp1 = vt4;
				vp2 = vt2;
			}

			for (block_i = 1; block_i<=gd->nBlocks; block_i++){
				/*for each pipe, find the two blocks to 
				which it corresponds*/
				v_b_start = (int)gd->vindex[block_i][1];
				v_b_end = (int)gd->vindex[block_i][2];


				/*find the block range that the two endpoints are within.
				Should only be one block, so if the for loop is not 
				on that block, it will continue and test the next one*/
				/*need something different for measured points*/
				if (vp1 >= 0 && vp2 >= 0){
					if (v_b_start <= vp1 && vp1 <= v_b_end &&
						v_b_start <= vp2 && vp2 <= v_b_end){
						/*block is the correct one, assign x and y*/		

						/*boundaries, don't necessarily follow
						the correct order... with the blocks,
						the larger index should always go first*/
						if(vt1 == vt2 && vt3 == vt4){
							/*first, deal with block index overlap*/
							if ((vp1 == v_b_start) && (vp2 != (vp1+1)))
								vp1 = v_b_end+1;
							else if ((vp2 == v_b_start) && (vp1 != (vp2+1)))
								vp2 = v_b_end+1;

							/*now put points in appropriate order*/
							if (vp2 < vp1){
								vt1 = vp2;
								vp2 = vp1;
								vp1 = vt1;
							}
						}

						x1 = gd->vertices[vp1][1];
						y1 = gd->vertices[vp1][2];
						x2 = gd->vertices[vp2][1];
						y2 = gd->vertices[vp2][2];

						node1 = fd->vertNodeIndex[vp1];
						node2 = fd->vertNodeIndex[vp2];

						p1 = fd->pressures[node1];
						p2 = fd->pressures[node2];
					}
					else {
						continue;
					}
				}
				/*deal with the measured points*/
				else if (vp1 < 0 || vp2 < 0) {
					if (vp1 < 0 && vp2 >= 0){
						/*check if block_i matches the block for vp1*/
						if(gd->points[-vp1][3] != block_i)
							continue;
						
					    /*didn't continue, assign xy, node and pressure */
						x1 = gd->points[-vp1][1];
						y1 = gd->points[-vp1][2];
						x2 = gd->vertices[vp2][1];
						y2 = gd->vertices[vp2][2];

						node1 = fd->mpNodeIndex[-vp1];
						node2 = fd->vertNodeIndex[vp2];

						p1 = fd->pressures[node1];
						p2 = fd->pressures[node2];

					}
					else if (vp2 < 0 && vp2 >= 0){
						/*check if block_i matches the block for vp2*/
						if(gd->points[-vp2][3] != block_i)
							continue;

					    /*didn't continue, assign xy */
						x1 = gd->vertices[vp1][1];
						y1 = gd->vertices[vp1][2];
						x2 = gd->points[-vp2][1];
						y2 = gd->points[-vp2][2];

						node1 = fd->vertNodeIndex[vp1];
						node2 = fd->mpNodeIndex[-vp2];

						p1 = fd->pressures[node1];
						p2 = fd->pressures[node2];
					}
					else if (vp1 < 0 && vp2 < 0){
						/*blocks should be the same, check with vp1 */
						if(gd->points[-vp1][3] != block_i)
							continue;

						x1 = gd->points[-vp1][1];
						y1 = gd->points[-vp1][2];
						x2 = gd->points[-vp2][1];
						y2 = gd->points[-vp2][2];

						node1 = fd->mpNodeIndex[-vp1];
						node2 = fd->mpNodeIndex[-vp2];

						p1 = fd->pressures[node1];
						p2 = fd->pressures[node2];
					}
				}

				/*should be on the correct block at this point*/
				/*calculate the length of the block segment, L
				and the resultant force, F */


				L = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
				F = n*(p1+p2)*L/2; /*divide by 1000, since forces
									   go into the F matrix in MN*/
				
				if (F != 0){
					/* calculate centroid of the force load (KAP 1999 22-24) */
					if (p1 > p2){
						F1 = n*p2*L;
						F2 = n*1/2*(p1-p2)*L;
						d = (F1*(L/2)+F2*(L/3))/F;
					}
					else if (p2 > p1) {
						F1 = n*p1*L;
						F2 = n*1/2*(p2-p1)*L;
						d = (F1*(L/2)+F2*(2*L/3))/F;
					}
					else if (p2 == p1) {
						d = 0.5*L;
					}

					xf = x1 + d/L*(x2-x1);
					yf = y1 + d/L*(y2-y1);

					/* KAP 1999, eqns 25-28*/
					/* these should account for every case, since the indices
					more counter clockwise around the blocks*/
					if (x2 > x1){
						/*alpha = atan2((y2-y1),(x2-x1));*/
						alpha = atan((y2-y1)/(x2-x1));
						Fx = -F*sin(alpha);
						Fy = F*cos(alpha);
					}
					else if (x1 > x2){
						alpha = atan((y2-y1)/(x2-x1));
						Fx = F*sin(alpha);
						Fy = -F*cos(alpha);
					}
					else if (x1 == x2 && y2 > y1){
						Fx = -F;
						Fy = 0;
					}
					else if (x1 == x2 && y1 > y2){
						Fx = F;
						Fy = 0;
					}

					transmap(blockArea[block_i],T,xf,yf);

					/*indexing for blocks in the force vector */
					i2 = k1[block_i];

					/*KAP 1999 eqn 32*/
					for (force_i = 1; force_i <= 6; force_i++){
						ad->F[i2][force_i] += (Fx*T[1][force_i] + Fy*T[2][force_i]);

					}/*force_i*/

					/*if it reaches this point, it found the right block, 
					can break on the for loop for the rest of the blocks 
					for this pair of points*/
				}

				if(p1 > 0 || p2 > 0 )
				p1 = p1; /*stop for debugging*/


				break;
					
			}/*block_i*/
		}/*block_1or2*/
	}/*pipe_i*/
} /* fluid_pressure_load */

/* This was moved from ddanalysis to shorten the size
 * of the analysis loop.
 */
void
assemble(Geometrydata * gd, Analysisdata * ad,
         int ** locks, double ** e0,
         int * k1, int * kk, int ** n, double ** U,
         TransMap transmap, Fluidsdata * fd) {

   //double ** f = AData->F;
   double ** moments = gd->moments;

/** @todo Pull the globaltime array out of here
 * and just pass the current time into the 
 * interpolation function.
 */

   df09(ad->loadpoints,
          gd->points,
          ad->globalTime,
          ad->timeDeps,
          ad->tindex,
          gd->nFPoints,
          ad->nLPoints,
          ad->cts,
          ad->delta_t);

     /* Initialize submatrices back to zero, set
      * some other parameter related to the global
      * time step.
      * FIXME: see if locks can be moved out of here.
      */
	   df10(gd,ad,locks,ad->F);

	  /* submatrix of fixed points.  */
      df12(gd,ad,k1,ad->F,moments,n,transmap);

     /* df13 submatrix of stiffness */
	   stress_stiffness(gd->nBlocks,ad->K,k1,e0,moments,n,ad->planestrainflag);

     /* df14 submatrix of initial stress  */
	   stress_initial(gd->nBlocks,k1,ad->F,e0,moments);

     /* submatrix of point loading  */
	   df15(gd,ad, k1,ad->F,moments, transmap);

     /* seismic loading */
      if (ad->timehistory != NULL) {

         seismicload(gd->seispoints, ad->timehistory, ad->cts,
                     k1, ad->F, moments, e0, transmap);
      }

     /* submatrix of volume force */
	   df16(gd,k1,ad->F,e0,moments, fd, ad);
     
	 /* bolt stiffness */
      if (gd->nBolts > 0) {

         bolt_stiffness_a(gd->rockbolts, gd->nBolts, ad->K, 
            k1, kk, n, moments,ad->F,transmap);
      }

	  /*line loads from pipe network fluids */
	  if (ad->pipeflowflag == 1){
		  fluid_pressure_load(fd, gd, ad, moments, k1, transmap);		          
	  }
}  