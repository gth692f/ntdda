/*files for calculating pipe geometry, added by WEM */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "ddamemory.h"
#include "geometry.h"
#include "dda.h" // for dda_display_warning

int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }

void  
producePipeSegments(Geometrydata *gd, Fluidsdata *fd){
	
	int d1, d2, d3; /*dummy variables*/
	int block_i, vi, line_i, obsegment_i, obsegment_i2; /*index variables */
	int n0; /* number of blocks */
	int n3; /* number of points */
	int nlines, n_obsegments, n_psegments, n_new_MPoints;
	int vi_start, vi_end;
	int new_line_flag, new_seg_flag;
	int v1a, v1b, v2a, v2b; /* 1,2 - point #; a,b - block # */
	double mb_tolerance;

	int pipesize1;
    int pipesize2;
   
	unsigned long long inf = 0x7ff0000000000000; /*to define a positive infinity, which is giving trouble*/

	double m, b; /*slope and y intercept */
	double x1, y1, x2, y2;
	double x1a, x1b, y1a, y1b, x2a, x2b, y2a, y2b, ml1, ml2, bl1, bl2;

	int ** vindex;
	double ** vertices;
	double ** lines; /* for each line [# of blocks on it, slope, y-intercept, block #s from v index...] */
	double ** original_block_segments; /* [# of matches, x1, y1, x2, y2, v1, v2, m, b, block#]*/
		/*# of matches, will be zero for edges, 1 for regular pairs, 1 or more for irregular pairs..?*/
    double ** pipe_segments; /* [0 v1a v1b v2a v2b width]; */
	double ** points = gd->points;

	int inRangeExclusive(double x3, double y3, double x4a, double y4a, double x4b, double y4b);
	void newPipeSegment(double ** pipe_segments, int n_psegments, int v1a, int v1b, int v2a, int v2b);
	void assignMPoints (Geometrydata *gd, double ** pipe_segments, int n_psegments, double ** vertices,
							int ** vindex, double v1a_d, double v1b_d, double v2a_d, double v2b_d,
							int v1a_o, int v1b_o, int v2a_o, int v2b_o);


	
	n0 = gd->nBlocks;
	n3 = gd->n3;
	nlines = 1;
	n_obsegments = 0; 
	n_psegments = 0;
	n_new_MPoints = 0;
	vi_start = vi_end = 0;
	mb_tolerance = 0.000000001;
	
	/*information about which block was attached to the original
	line is not currently needed for anything, uncomment
	this line out if it is at any point*/
	/*lines = DoubMat2DGetMem((n3*(n3-1)/2+1),gd->nBlocks+3);*/
	lines = DoubMat2DGetMem((n3*(n3-1)/2+1),3);
	original_block_segments = DoubMat2DGetMem((n3*(n3-1)/2+1),10);

	pipesize1 = (n3*(n3-1)/2+1);
	pipesize2 = 8;

	pipe_segments = DoubMat2DGetMem(pipesize1,pipesize2);

	vindex = gd->vindex;
	vertices = gd->vertices;
	
	/* assign a line number to every line formed from the vertices of every block*/
	for (block_i = 1; block_i<=n0; block_i++){
		vi_start = vindex[block_i][1];
		vi_end = vindex[block_i][2];
		for (vi = vi_start; vi<=vi_end; vi++){
			
			new_line_flag = 1;

			x1 = vertices[vi][1];
			y1 = vertices[vi][2];
			x2 = vertices[vi+1][1];
			y2 = vertices[vi+1][2];

			m = (y2-y1)/(x2-x1);
			
			/*check for infinity case*/
			if (isinf(m)){
				m = *( double* )&inf ; /*forcing positive infinity */
				b = x2; /*x intercept*/
			}
			else
		    b = y2 - x2*m;

			/* checks if line is already present for line matching*/
			for (line_i = 1; line_i <= nlines; line_i++){
				/* if this doesn't work due to floating point errors, use abs > 0.0000001)*/
				if (lines[line_i][0] > 0 && lines[line_i][1] == m && lines[line_i][2] == b){
					lines[line_i][0]++;
					d1 = (int)lines[line_i][0];
					/*lines[line_i][d1+2] = block_i;*/
					/*information about which block was attached to the original
					line is not currently needed for anything, uncomment
					this line out if it is at any point*/

					new_line_flag = 0; /*since no new line was added*/
					break;
					
				} /*if*/	

			}/*line_i*/

			if (new_line_flag == 1){
				lines[nlines][0] = 1;
				lines[nlines][1] = m;
				lines[nlines][2] = b;
				/*lines[nlines][3] = block_i;*/
				/*information about which block was attached to the original
				line is not currently needed for anything, uncomment
				this line out if it is at any point*/
				

				nlines++;
			} /* if*/

			/*creates a new segment based on these values*/
			
			n_obsegments++;

			original_block_segments[n_obsegments][1] = x1;
			original_block_segments[n_obsegments][2] = y1;
			original_block_segments[n_obsegments][3] = x2;
			original_block_segments[n_obsegments][4] = y2;
			original_block_segments[n_obsegments][5] = vi;
			original_block_segments[n_obsegments][6] = vi+1;
			original_block_segments[n_obsegments][7] = m;
			original_block_segments[n_obsegments][8] = b;
			original_block_segments[n_obsegments][9] = block_i;
			
			/* resetting vi_end+1 to vi_start, since I think that is all that's tracked */
			if (vi == vi_end)
				original_block_segments[n_obsegments][6] = vi_start;



		}/* vi */
	
	}/* block_i*/

	/*Have to check a lot of conditions, based on the different possible geometries among the blocks.  
	Hence the large number of for loops and if statements*/
	for (obsegment_i = 1; obsegment_i < n_obsegments; obsegment_i++){
		x1a = original_block_segments[obsegment_i][1];
		y1a = original_block_segments[obsegment_i][2];
		x2a = original_block_segments[obsegment_i][3];
		y2a = original_block_segments[obsegment_i][4];

		v1a = (int)original_block_segments[obsegment_i][5];
		v2a = (int)original_block_segments[obsegment_i][6];

		ml1 = original_block_segments[obsegment_i][7];
		bl1 = original_block_segments[obsegment_i][8];

		new_seg_flag = 0;

		for (obsegment_i2 = obsegment_i+1; obsegment_i2 <= n_obsegments; obsegment_i2++){
			x1b = original_block_segments[obsegment_i2][3];
			y1b = original_block_segments[obsegment_i2][4];
			x2b = original_block_segments[obsegment_i2][1];
			y2b = original_block_segments[obsegment_i2][2];

			v1b = (int)original_block_segments[obsegment_i2][6];/*note 6 and 5 switched from v1a, b*/
			v2b = (int)original_block_segments[obsegment_i2][5];

			ml2 = original_block_segments[obsegment_i2][7];
			bl2 = original_block_segments[obsegment_i2][8];

			/*Find all perfect matches (fairly certain points will be on opposite sides..?)*/
			if (x1a == x1b && y1a == y1b && x2a == x2b && y2a == y2b){

				original_block_segments[obsegment_i][0]++;
				original_block_segments[obsegment_i2][0]++;
				n_psegments++;

				newPipeSegment(pipe_segments, n_psegments, v1a, v1b, v2a, v2b);



				continue;

			}/* if */

			/*Find all lines with one point matching*/
			
			/*first point matches, second doesn't*/
			if ((x1a == x1b && y1a == y1b) && (x2a != x2b || y2a != y2b)){
				/*First, check slope and intercept*/
				if (((fabs(ml1 - ml2) < mb_tolerance) || (isinf(ml1) && isinf(ml2))) && fabs(bl1-bl2) < mb_tolerance){
				/*Next, check if points are in the range defined by the other line.
				 Measure point assignment, if applicable, depends on which point is outside which line*/
					if (inRangeExclusive(x2a, y2a, x1b, y1b, x2b, y2b)){
						/*Assign measure point to x2b, y2b*/

						n_new_MPoints++;
						original_block_segments[obsegment_i][0]++;
						original_block_segments[obsegment_i2][0]++;
						n_psegments++;
						
						newPipeSegment(pipe_segments, n_psegments, v1a, v1b, v2a, -n_new_MPoints);
						assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
							pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
							pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
							v1a, v1b, v2a, v2b);



						continue;

					}/* if */
					else if (inRangeExclusive(x2b, y2b, x1a, y1a, x2a, y2a)){
						/*Assign measure point to x2a, y2a*/

						n_new_MPoints++;
						original_block_segments[obsegment_i][0]++;
						original_block_segments[obsegment_i2][0]++;
						n_psegments++;

						newPipeSegment(pipe_segments, n_psegments, v1a, v1b, -n_new_MPoints, v2b);
						assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
							pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
							pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
							v1a, v1b, v2a, v2b);



						continue;

					}/* else if*/
						
				/*Finally, find the distance and assign a measured point to the appropriate block*/

				} /* if */
			}/* if */

			/* second point matches, first doesn't */
			else if ((x1a != x1b || y1a != y1b) && (x2a == x2b && y2a == y2b)){

				if (((fabs(ml1 - ml2) < mb_tolerance) || (isinf(ml1) && isinf(ml2))) 
					&& fabs(bl1-bl2) < mb_tolerance){
				/*Next, check if points are in the range defined by the other line.
				 Measure point assignment, if applicable, depends on which point is outside which line*/
					if (inRangeExclusive(x1a, y1a, x1b, y1b, x2b, y2b)){
						/*Assign measure point to x1b, y1b*/
						/*maybe make these negative, to distinguish them from regular vertex indices*/

						n_new_MPoints++;
						original_block_segments[obsegment_i][0]++;
						original_block_segments[obsegment_i2][0]++;
						n_psegments++;

						newPipeSegment(pipe_segments, n_psegments, v1a, -n_new_MPoints, v2a, v2b);
						assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
							pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
							pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
							v1a, v1b, v2a, v2b);

						continue;

					}/* if */
					else if (inRangeExclusive(x1b, y1b, x1a, y1a, x2a, y2a)){
						/*Assign measure point to x1a, y1a*/

						n_new_MPoints++;
						original_block_segments[obsegment_i][0]++;
						original_block_segments[obsegment_i2][0]++;
						n_psegments++;

						newPipeSegment(pipe_segments, n_psegments, -n_new_MPoints, v1b, v2a, v2b );
						assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
							pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
							pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
							v1a, v1b, v2a, v2b);

						continue;

					}/* else if*/
						
				/*Finally, find the distance and assign a measured point to the appropriate block*/

				} /* if */
			}/* else if */

			/*Finally, check if no matches*/
			else if ((x1a != x1b || y1a != y1b) && (x2a != x2b || y2a != y2b)){
				/*First, check slope and intercept*/
				if (((fabs(ml1 - ml2) < mb_tolerance) || (isinf(ml1) && isinf(ml2)))
					&& fabs(bl1-bl2) < mb_tolerance){
					/*means we need to add two measured points, 6 scenarios*/
					/*checks if points are oriented correctly on the line
					First four, blocks are staggered, edge of one block hangs 
					over the edge of the other, such that no points would be registered
					but a pipe exists between.  Last two, one block is completely on top 
					of the other, same problem*/
					if (inRangeExclusive(x1a, y1a, x1b, y1b, x2b, y2b) &&
						inRangeExclusive(x1b, y1b, x1a, y1a, x2a, y2a) &&
						inRangeExclusive(x1a, y1a, x2a, y2a, x2b, y2b) &&
						inRangeExclusive(x1b, y1b, x2a, y2a, x2b, y2b)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;

							newPipeSegment(pipe_segments, n_psegments, v1a, v1b, d2, d3);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* if */
					else if (inRangeExclusive(x1b, y1b, x1a, y1a, x2a, y2a) &&
						inRangeExclusive(x2a, y2a, x1b, y1b, x2b, y2b) &&
						inRangeExclusive(x1b, y1b, x1a, y1a, x2b, y2b) &&
						inRangeExclusive(x2a, y2a, x1a, y1a, x2b, y2b)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;

							newPipeSegment(pipe_segments, n_psegments, d2, v1b, v2a, d3);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* else if */
					else if (inRangeExclusive(x2a, y2a, x1b, y1b, x2b, y2b) &&
						inRangeExclusive(x2b, y2b, x1a, y1a, x2a, y2a) &&
						inRangeExclusive(x2a, y2a, x1a, y1a, x1b, y1b) &&
						inRangeExclusive(x2b, y2b, x1a, y1a, x1b, y1b)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;

							newPipeSegment(pipe_segments, n_psegments, d2, d3, v2a, v2b);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* else if */
					else if (inRangeExclusive(x2b, y2b, x1a, y1a, x2a, y2a) &&
						inRangeExclusive(x1a, y1a, x1b, y1b, x2b, y2b) &&
						inRangeExclusive(x2b, y2b, x1b, y1b, x2a, y2a) &&
						inRangeExclusive(x1a, y1a, x1b, y1b, x2a, y2a)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;

							newPipeSegment(pipe_segments, n_psegments, v1a, d2, d3, v2b);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* else if */
					else if (inRangeExclusive(x1a, y1a, x1b, y1b, x2b, y2b) &&
						inRangeExclusive(x2a, y2a, x1b, y1b, x2b, y2b)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;

							newPipeSegment(pipe_segments, n_psegments, v1a, d2, v2a, d3);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* else if */
					else if (inRangeExclusive(x1b, y1b, x1a, y1a, x2a, y2a) &&
						inRangeExclusive(x2b, y2b, x1a, y1a, x2a, y2a)){

							d2 = -n_new_MPoints++;
							d3 = -n_new_MPoints++;

							original_block_segments[obsegment_i][0]++;
							original_block_segments[obsegment_i2][0]++;
							n_psegments++;


							newPipeSegment(pipe_segments, n_psegments, d2, v1b, d3, v2b);
							assignMPoints (gd, pipe_segments, n_psegments, vertices, vindex, 
								pipe_segments[n_psegments][1], pipe_segments[n_psegments][2], 
								pipe_segments[n_psegments][3], pipe_segments[n_psegments][4],
								v1a, v1b, v2a, v2b);

							continue;

					}/* else if */

				}/* if */

			
			}/* else if*/

			/*In case there's some state not accounted for above*/
			else{
				dda_display_warning("Unaccounted for case in pipe formation from geometry file");
					/*perhaps could include something that would allow the geometry to be produced, 
					just without the fluid flow module loaded? */
			}/*else*/

		}/*obsegment_i2*/
	}/*obsegment_i*/

	/*now, create pipes for the boundary pipes*/
	for (obsegment_i = 1; obsegment_i < n_obsegments; obsegment_i++){
		
		/* if true, segment was a boundary (unused) */
		if (original_block_segments[obsegment_i][0] == 0){
			v1a = (int)original_block_segments[obsegment_i][5];
			v2a = (int)original_block_segments[obsegment_i][6];
			/* label a boundary by reduplicating the indices
		    for the new pipe formed */
			n_psegments++;
			newPipeSegment(pipe_segments, n_psegments, v1a, v1a, v2a, v2a);

		}/* if */
	}/*obsegment_i*/

	/*Assign pipe data to the fluids struct*/
	fd->pipesize1 = n_psegments+1;
	fd->pipesize2 = pipesize2;
	fd->nPSegments = n_psegments; 
	fd->pipeSegments = DoubMat2DGetMem(pipesize1,pipesize2);

	for (d1 = 0; d1 < fd->pipesize1; d1++){
		for (d2 = 0; d2 < fd->pipesize2; d2++){
			fd->pipeSegments[d1][d2] = pipe_segments[d1][d2];		
		}/*d2*/
	}/*d1*/

	free2DMat((void **)lines, (n3*(n3-1)/2+1)); 
	free2DMat((void **)original_block_segments, (n3*(n3-1)/2+1));
	free2DMat((void **)pipe_segments, pipesize1);

}/* close producePipeSegments */

int inRangeExclusive(double xa, double ya, double x3b, double y3b, double x4b, double y4b){
	
	/* Determines if point 3 lies on the line defined by 4a and 4b.
	Only use if it is already established that the slopes of the lines are the same
	function does not return true if x3 is the same point as either x4*/
	
	double dist_a3b, dist_a4b, dist_34;
	
	/*case, is x3, y3 contained within [x4a, y4a, x4b, y4b]*/
	dist_a3b = sqrt(pow((x3b-xa),2)+pow((y3b-ya),2));
	dist_a4b = sqrt(pow((x4b-xa),2)+pow((y4b-ya),2));
	dist_34 = sqrt(pow((x4b-x3b),2)+pow((y4b-y3b),2));

	if (dist_a3b < dist_34 && dist_a4b < dist_34){
		/* point is within the line*/
		return 1;}
	else
		return 0;
	
} /* close inRangeExclusive */


int inRangeInclusive(double xa, double ya, double x3b, double y3b, double x4b, double y4b){
	
	/* Determines if point 3 lies on the line defined by 4a and 4b.
	Only use if it is already established that the slopes of the lines are the same
	function does return true if x3 is the same point as either x4*/
	
	double dist_a3b, dist_a4b, dist_34;
	
	/*case, is x3, y3 contained within [x4a, y4a, x4b, y4b]*/
	dist_a3b = sqrt(pow((x3b-xa),2)+pow((y3b-ya),2));
	dist_a4b = sqrt(pow((x4b-xa),2)+pow((y4b-ya),2));
	dist_34 = sqrt(pow((x4b-x3b),2)+pow((y4b-y3b),2));

	if (dist_a3b <= dist_34 && dist_a4b <= dist_34){
		/* point is within the line*/
		return 1;}
	else
		return 0;
	
} /* close inRangeExclusive */


void newPipeSegment(double ** pipe_segments, int n_psegments, int v1a, int v1b, int v2a, int v2b){
	/* adds a new pipe segment with indices v1a, v1b, v2a, and v2b*/
		
	pipe_segments[n_psegments][1] = v1a;
	pipe_segments[n_psegments][2] = v1b;
	pipe_segments[n_psegments][3] = v2a;
	pipe_segments[n_psegments][4] = v2b;
	
	} /* newPipeSegment */


void assignMPoints (Geometrydata *gd, double ** pipe_segments, int n_psegments, double ** vertices, 
					int ** vindex, double v1a_d, double v1b_d, double v2a_d, double v2b_d, int v1a_o, int v1b_o, int v2a_o, int v2b_o){
	/* Currently, there are a lot of extra measured points.  May want to leave them in for now, though,
	so that they can be changed to adequately reflect the joint positions later. */

	int pointCount, nPoints, nMPoints; /*number of total points, with and without hole points respectively*/
	int block_i;
	double ** points;
	double x, y;
	int v1a, v1b, v2a, v2b;

	nMPoints = gd->nMPoints;
	nPoints = gd->nPoints;
	pointCount = gd->pointCount;
	points = gd->points;

	v1a = (int)v1a_d;
	v1b = (int)v1b_d;
	v2a = (int)v2a_d;
	v2b = (int)v2b_d;

	/*add measure points at intersections that need them. Measure points are at the end of the points list,
	so each time they are added in at pointCount + 1*/
	if (v1a < 0){
		pointCount++;
		nPoints++;
		nMPoints++;
		
		x = vertices[v1b][1];
		y = vertices[v1b][2];

		pipe_segments[n_psegments][1] = -pointCount; /*again, negative will serve in the future to say which track measured points*/

		points[pointCount][0] = 1;
		points[pointCount][1] = x;
		points[pointCount][2] = y;
		
		/* find the block that holds v1a */
		for(block_i = 1; block_i <= gd->nBlocks; block_i++){
			if(vindex[block_i][1] <= v1a_o && v1a_o <= vindex[block_i][2]){
				points[pointCount][3] = block_i;
				break;
			}/*if*/
		} /*block_i*/
	}/* if */
	if (v1b < 0){
		pointCount++;
		nPoints++;
		nMPoints++;

		x = vertices[v1a][1];
		y = vertices[v1a][2];
		pipe_segments[n_psegments][2] = -pointCount; /*again, negative will serve in the future to say which track measured points*/

		points[pointCount][0] = 1;
		points[pointCount][1] = x;
		points[pointCount][2] = y;
		
		/* find the block that holds v1b */
		for(block_i = 1; block_i <= gd->nBlocks; block_i++){
			if(vindex[block_i][1] <= v1b_o && v1b_o <= vindex[block_i][2]){
				points[pointCount][3] = block_i;
				break;
			}/*if*/
		} /*block_i*/
	}/* if */
	if (v2a < 0){
		pointCount++;
		nPoints++;
		nMPoints++;

		x = vertices[v2b][1];
		y = vertices[v2b][2];
		pipe_segments[n_psegments][3] = -pointCount; /*again, negative will serve in the future to say which track measured points*/

		points[pointCount][0] = 1;
		points[pointCount][1] = x;
		points[pointCount][2] = y;
		
		/* find the block that holds v2a */
		for(block_i = 1; block_i <= gd->nBlocks; block_i++){
			if(vindex[block_i][1] <= v2a_o && v2a_o <= vindex[block_i][2]){
				points[pointCount][3] = block_i;
				break;
			}/*if*/
		} /*block_i*/
	}/* if */
	if (v2b < 0){
		pointCount++;
		nPoints++;
		nMPoints++;

		x = vertices[v2a][1];
		y = vertices[v2a][2];
		pipe_segments[n_psegments][4] = -pointCount; /*again, negative will serve in the future to say which track measured points*/

		points[pointCount][0] = 1;
		points[pointCount][1] = x;
		points[pointCount][2] = y;
		
		/* find the block that holds v2b */
		for(block_i = 1; block_i <= gd->nBlocks; block_i++){
			if(vindex[block_i][1] <= v2b_o && v2b_o <= vindex[block_i][2]){
				points[pointCount][3] = block_i;
				break;
			}/*if*/
		} /*block_i*/
	}/* if */

	gd->nPoints = nPoints;
	gd->pointCount = pointCount;
	gd->points = points;
	gd->nMPoints = nMPoints;

} /* assignMPoints */


void
calcPipeWidths(Geometrydata * gd, Fluidsdata *fd){

	/*calculates the widths between each pipe, for use in the flow calculations
	Width is calculated as the average across a pipe based on the location of 
	its endpoints. eg, for the pipe between rocks with edge indices 1-3 and 2-4:
	 1         3
	*-----------*
	 2         4
	 the width is the average of the distances between (1 and 2) and (3 and 4)*/

	/* should be three different cases:
	a) at least one index = 0: the pipe is a boundary edge, 
	for which it keeps its original width (boundflag = 1)
	b) all indices > 0: the width is computed normally
	c) at least one index < 0: the edge has one or more measured points, 
	   in which case something special 
     	might need to happen for the geometry to be correct. (measflag += 1) */
	
	int v1a, v1b, v2a, v2b;
	int n_psegments, pipe_i;
	double x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b;
	double dist1, dist2, dist_avg;
	double meanHydraulicAperture, r;

	double ** pipe_segments;
	double ** points;
	double ** vertices;
	
		
	pipe_segments = fd->pipeSegments;
	points = gd->points;
	n_psegments = fd->nPSegments;
	vertices = gd->vertices;

	for (pipe_i = 1; pipe_i <= n_psegments; pipe_i++){

		/* if one, width is preassigned a constant value,
		go on to the next pipe and leaves the hydraulic
		aperture the same as it was*/
		if (pipe_segments[pipe_i][0] == 1)
			continue;

		v1a = (int)pipe_segments[pipe_i][1];
		v1b = (int)pipe_segments[pipe_i][2];
		v2a = (int)pipe_segments[pipe_i][3];
		v2b = (int)pipe_segments[pipe_i][4];

		/* get appropriate x-y coordinates for each vertex or measured point
		measured points created for the pipe segments use negative indices, 
		hence the if statments */
        
		if (v1a != v1b && v2a != v2b){
			if (v1a > 0){
				x1a = vertices[v1a][1];
				y1a = vertices[v1a][2];
			}/* if */
			else if (v1a < 0){
				x1a = points[-v1a][1];
				y1a = points[-v1a][2];
			}/* else if */

			if (v1b > 0){
				x1b = vertices[v1b][1];
				y1b = vertices[v1b][2];
			}/* if */
			else if (v1b < 0){
				x1b = points[-v1b][1];
				y1b = points[-v1b][2];
			}/* else if */

			if (v2a > 0){
				x2a = vertices[v2a][1];
				y2a = vertices[v2a][2];
			}/* if */
			else if (v2a < 0){
				x2a = points[-v2a][1];
				y2a = points[-v2a][2];
			}/* else if */

			if (v2b > 0){
				x2b = vertices[v2b][1];
				y2b = vertices[v2b][2];
			}/* if */
			else if (v2b < 0){
				x2b = points[-v2b][1];
				y2b = points[-v2b][2];
			}/* else if */

			dist1 = sqrt(pow((x1a-x1b),2)+pow((y1a-y1b),2));
			dist2 = sqrt(pow((x2a-x2b),2)+pow((y2a-y2b),2));
			dist_avg = (dist1+dist2)/2;

			/*Jing 2001 eqn 7*/
			if (dist1 == 0 && dist2 == 0){
				r = 1e-20;
			}
			else {
				if (dist1 <= dist2)
					r = dist1/dist2;
				else
					r = dist2/dist1;
			}

			meanHydraulicAperture = dist_avg*pow(16*pow(r,2)/pow(1+r,4),1/3);
			
			pipe_segments[pipe_i][5] = meanHydraulicAperture; 
		}/* if not a boundary pipe*/
		else{
			pipe_segments[pipe_i][5] = pipe_segments[pipe_i][5];
			/*keeps the original value, don't really need this line,
			just included in case we need to do something else with 
			the boundary lines later */
		}

	}/*pipe_i*/

}/* end calcPipeWidths */


void
initializePipeProperties(Geometrydata *gd, Fluidsdata *fd){
	/*find a matching pair of two positive indices, use those to label the joint type*/
	/*this code takes the ddaml fluid property inputs for each joint type and matches them 
	with the newly developed pipes. Assigns the initial width, friction factor and a constant
	saying whether or not the widths change*/ 

	double ** pipe_segments;
	
	int  pipe_i, v_i;
	int n_psegments, n_negative;
	int v4indexes[4];
	int joint_type;

	double findJointTypeZeroToOneNeg(Geometrydata *gd, int * v4ind);
	double findJointTypeTwoNeg(Geometrydata *gd, int * v4ind);

	pipe_segments = fd->pipeSegments;
	n_psegments = fd->nPSegments;

	for (pipe_i = 1; pipe_i <= n_psegments; pipe_i++){

		v4indexes[0] = (int)pipe_segments[pipe_i][1];
		v4indexes[1] = (int)pipe_segments[pipe_i][2];
		v4indexes[2] = (int)pipe_segments[pipe_i][3];
		v4indexes[3] = (int)pipe_segments[pipe_i][4];

		n_negative = 0;

		/* count # of negative values */
		for (v_i = 0; v_i <= 3; v_i++){
			if (v4indexes[v_i] < 0)
				n_negative++;
		}/*v_i*/
		
		/*pipe type found based on properties of the vertex indices that 
		make up its endpoints.  Algorithm for assignment depends on if 
		zero, one or two endpoints have negative indices (implying measured
		points)*/
		switch (n_negative) {
			case 0: 
				joint_type = (int)findJointTypeZeroToOneNeg(gd, v4indexes);
				break;
			case 1:
				joint_type = (int)findJointTypeZeroToOneNeg(gd, v4indexes);
				break;
			case 2:
				joint_type = (int)findJointTypeTwoNeg(gd, v4indexes);
				break;
			default:
				joint_type = 0;
				break;
				   
				
		}/*n_negative*/
		
		/* if 0 or 1 negative */

		/* assign a width based on this joint_type*/

	    /* WARNING: This assumes that joint types are listed in order 
	     * of occurrence in xml file.  The type attribute is ignored
	     * for now.
	     */

		/*if the joint_type is less than the number of defined types in the fluid file*/
		if (joint_type <= fd->nFluidJoints){
			pipe_segments[pipe_i][5] = fd->jointFluidProp[joint_type][0]; /*width*/
			pipe_segments[pipe_i][6] = fd->jointFluidProp[joint_type][1]; /*friction_factor*/
			pipe_segments[pipe_i][0] = fd->jointFluidProp[joint_type][2]; 
				/*boolean for whether or not width remains constant */
		}
		else{
			pipe_segments[pipe_i][5] = 0.00;
			pipe_segments[pipe_i][6] = 1.00;
		}


	}/*pipe_i*/

}/* initializePipeProperties */




double findJointTypeZeroToOneNeg(Geometrydata *gd, int * v4ind){

	int isPositiveConsecutive(Geometrydata *gd, int va, int vb);
		/* for each pair, check if one follows the other.
		Inclusion of the negative measured points in each case means that 
		it's not as simple as looking at the first two; need to go through 
		and check all combinations, front and back take the first one that works. 
		Only points on the same block will have consecutive numbers*/

	int v_io, v_ii, v1, v2;
	double joint_type;

	double ** vertices;

	vertices = gd->vertices;

	/*default, in case nothing happens */
	joint_type = 0;

	for (v_io = 0; v_io <= 2; v_io++){
		for (v_ii = 1; v_ii <=3; v_ii++){
			
			v1 = v4ind[v_io];
			v2 = v4ind[v_ii];
			
			if (v1 > 0 && v2 > 0){
				if (isPositiveConsecutive(gd,v1,v2))
					joint_type = vertices[v2][0];
				else if (isPositiveConsecutive(gd,v2,v1))
					joint_type = vertices[v1][0];
			}/*if*/

		}/*v_ii*/
	}/*v_io*/
				
	return joint_type;
}/*findJointTypeZeroToOneNeg*/




double findJointTypeTwoNeg(Geometrydata *gd, int * v4ind){

	int isPositiveConsecutive(Geometrydata *gd, int va, int vb);
		/*Two possibilities here that I see.  Either two points
		exist that are positive consecutive, or two that aren't.
		for those that aren't, it looks like the joint type
		is stored in either one of the vertex indices+1.*/

	int v_i, v_io, v_ii, v1, v2;
	double joint_type;
	
	double x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b, ma, mb, ba, bb;
	double mb_tolerance;


	double ** vertices;
	unsigned long long inf = 0x7ff0000000000000; /*to define a positive infinity*/

	vertices = gd->vertices;
	mb_tolerance = 0.00000001;

	/*default, in case nothing happens */
	joint_type = 0;

	/*first case, two consecutive points exist*/
	for (v_io = 0; v_io <= 2; v_io++){
		for (v_ii = 1; v_ii <=3; v_ii++){
			
			v1 = v4ind[v_io];
			v2 = v4ind[v_ii];

			if (v1 > 0 && v2 > 0){
				if (isPositiveConsecutive(gd,v1,v2)){
					joint_type = vertices[v2][0];
					return joint_type;
				}/* if */
				else if (isPositiveConsecutive(gd,v2,v1)){
					joint_type = vertices[v1][0];
					return joint_type;
				}/* else if */
			
			}/*if*/
		}/*v_ii*/
	}/*v_io*/
				
	/*second case, no consecutive points exist*/
	/*here, find the first two positive numbers, and 
	use slope matching and intercepts to determine
	whether to get the joint_type from v_i or v_i+1*/

	/*assign positive numbers to v1 and v2*/
	v1 = 0; v2 = 0;
	for (v_i = 0; v_i <= 3; v_i++){
		if (v1 == 0 && v4ind[v_i] > 0){
			v1 = v4ind[v_i];
			continue;
		}/*if*/
		if (v4ind[v_i] > 0)
			v2 = v4ind[v_i];
	}/*v_i*/

	/*case 1, lines are shared by v_i and v_i+1
	Since only one other case is possible (v_i and v_i-1), 
	and the index for that would be controlled by v_i,
	checking case 1 first, and if not case 1 the case 2*/
	
	x1a = vertices[v1][1];
	y1a = vertices[v1][2];
	x2a = vertices[v1+1][1];
	y2a = vertices[v1+1][2];

	ma = (y2a-y1a)/(x2a-x1a);
			
	/*check for infinity case*/
	if (isinf(ma)){
		ma = *( double* )&inf ; /*forcing positive infinity */
		ba = x2a; /*x intercept*/
	}
	else
		ba = y2a - x2a*ma;

	x1b = vertices[v2][1];
	y1b = vertices[v2][2];
	x2b = vertices[v2+1][1];
	y2b = vertices[v2+1][2];

	mb = (y2b-y1b)/(x2b-x1b);
		
	/*check for infinity case*/
	if (isinf(mb)){
		mb = *( double* )&inf ; /*forcing positive infinity */
		bb = x2b; /*x intercept*/
	}
	else
		bb = y2b - x2b*mb;

	if (((fabs(mb-ma) < mb_tolerance) || (isinf(mb) && isinf(ma)))
		&& fabs(bb-ba) < mb_tolerance){
			joint_type = vertices[v1+1][0]; /* or v2+1, shouldn't matter */
			return joint_type;
	
	}/* if */
	else{
			joint_type = vertices[v1][0];
			return joint_type;
	}/* else*/




	return joint_type;
}/*findJointTypeTwoNeg*/




int isPositiveConsecutive(Geometrydata *gd, int va, int vb){
	
	int block_i;

	if(va > 0 && vb > 0){
		
		if (vb == (va+1))
			return 1;
		else{
			for (block_i = 1; block_i < gd->nBlocks; block_i++){
				/* want it to return 1 so that the joint number of the
				initial index gets assigned, not that of the latter */
				if (vb == gd->vindex[block_i][1] && va == gd->vindex[block_i][2])
					return 1;
			}/*block_i*/
		} /*else*/;
	}/*if*/
	
	/*if nothing was true */
	return 0;
}/* isPositiveConsecutive */



void produceNodes(Geometrydata *gd, Fluidsdata *fd){

	double ** vertices;
	double ** pipe_segments;
	double ** nodes;
	int ** connectivityTable;
	double xnode, ynode, xpipe, ypipe;

	int nodesize1, nodesize2, nNodes, newPipetoNodeInd;
	int connsize1, connsize2;
	int pipe_i, node_i, end_i, v_ind, d1, d2, block_i;
	int nodeAlreadyRegisteredFlag;
	int va, vb;

	vertices = gd->vertices;
	pipe_segments = fd->pipeSegments;
	
	/*Shouldn't have any more nodes than intersection points,
	and I think elsewhere in the code it says that no more than 
	10 lines may intersect at a point */

	nodesize1 = gd->nIntersectionPoints+1;
	nodesize2 = 3;
	nodes = DoubMat2DGetMem(nodesize1, nodesize2); 
	/*# of pipes at that node, node x coord, node y coord, 
	block or measured point indices that intersect there (up to 10)*/

	connsize1 = fd->nPSegments+1;
	connsize2 = 3;
	connectivityTable = IntMat2DGetMem(connsize1, connsize2);

	fd->vertNodeIndexsize1 = gd->vertexsize1+1;
    fd->vertNodeIndex = (int *)calloc(fd->vertNodeIndexsize1,sizeof(int));

	fd->mpNodeIndexsize1 = gd->nPoints+1;
	fd->mpNodeIndex = (int *)calloc(fd->mpNodeIndexsize1,sizeof(int));

	nNodes = 0;
    
	/*find the nodes connected to the two endpoints of each pipe*/
	for (pipe_i = 1; pipe_i <= fd->nPSegments; pipe_i++){
		
		/*nodes v1 and v2, should be at the same point*/
		for(end_i = 1; end_i <= 2; end_i++){

			nodeAlreadyRegisteredFlag = 0;

			/*assign v_ind, based on whichever index is not negative*/
			if (end_i == 1){
				if (pipe_segments[pipe_i][1] > 0)
					v_ind = 1;
				else
					v_ind = 2;
			}/* if */
			else{
				if (pipe_segments[pipe_i][3] > 0)
					v_ind = 3;
				else
					v_ind = 4;
			}/* if */		
		
			/*pipe_segments contains the index for the 
			x and y coordinates in vertices*/
			xpipe = vertices[(int)pipe_segments[pipe_i][v_ind]][1];
			ypipe = vertices[(int)pipe_segments[pipe_i][v_ind]][2];

			for (node_i = 1; node_i <= nNodes; node_i++){

				xnode = nodes[node_i][1];
				ynode = nodes[node_i][2];
				
				/*case, node already is registered, 
				just assigning another pipe to it*/
				if(xpipe == xnode && ypipe == ynode){
					
					nodes[node_i][0]++;

					nodeAlreadyRegisteredFlag = 1;

					connectivityTable[pipe_i][end_i] = node_i;

					if (end_i == 1){
						va = pipe_segments[pipe_i][1];
						vb = pipe_segments[pipe_i][2];
					}
					else {
						va = pipe_segments[pipe_i][3];
						vb = pipe_segments[pipe_i][4];
					}

					if (va > 0)
						fd->vertNodeIndex[va] = node_i;
					else
						fd->mpNodeIndex[-va] = node_i;

					if (vb > 0)
						fd->vertNodeIndex[vb] = node_i;
					else
						fd->mpNodeIndex[-vb] = node_i;

					fd->pipeSegments[pipe_i][end_i+7] = node_i;

					break;

				}/* if */
					
			}/* node_i*/

			if (nodeAlreadyRegisteredFlag == 0){
				nNodes++;

				nodes[nNodes][0]++;
				nodes[nNodes][1] = xpipe;
				nodes[nNodes][2] = ypipe;

				connectivityTable[pipe_i][end_i] = nNodes;

				/*Assign an original node to each of the 
				block indices (and deal with measured points*/

				if (end_i == 1){
					va = pipe_segments[pipe_i][1];
					vb = pipe_segments[pipe_i][2];
				}
				else {
					va = pipe_segments[pipe_i][3];
					vb = pipe_segments[pipe_i][4];
				}

				if (va > 0)
					fd->vertNodeIndex[va] = nNodes;
				else
					fd->mpNodeIndex[-va] = nNodes;

				if (vb > 0)
					fd->vertNodeIndex[vb] = nNodes;
				else
					fd->mpNodeIndex[-vb] = nNodes;

				fd->pipeSegments[pipe_i][end_i+7] = nNodes;
			}/* if */		


		}/* end_i */
	}/* pipe_i*/

	/*assign pressures to the indices that are doubled in the block indexing.
	e.g., a 4-sided block will have corner indices 1 2 3 4 5, with 1 and 5 
	being the same point.  This gives them the same pressure */

	for(block_i = 1; block_i<=gd->nBlocks; block_i++){
		va = gd->vindex[block_i][1];
		vb = gd->vindex[block_i][2];
		fd->vertNodeIndex[vb+1] = fd->vertNodeIndex[va];
	}/*block_i*/


	/* copy nodes and connectivity table into the fluidsdata struct */
	fd->nNodes = nNodes;
	fd->nodesize1 = nNodes+1;
	fd->nodesize2 = nodesize2;
	fd->nodes = DoubMat2DGetMem(fd->nodesize1, fd->nodesize2);

	for (d1 = 0; d1 < fd->nodesize1; d1++){
		for (d2 = 0; d2 < fd->nodesize2; d2++){
			fd->nodes[d1][d2] = nodes[d1][d2];		
		}/*d2*/
	}/*d1*/

	fd->connsize1 = connsize1;
	fd->connsize2 = connsize2;
	fd->connectivityTable = IntMat2DGetMem(fd->connsize1, fd->connsize2);

	for (d1 = 0; d1 < fd->connsize1; d1++){
		for (d2 = 0; d2 < fd->connsize2; d2++){
			fd->connectivityTable[d1][d2] = connectivityTable[d1][d2];		
		}/*d2*/
	}/*d1*/

	/*initialize heads and pressures for each node*/
    fd->headssize1 = fd->nNodes+1;
	fd->heads = (double *)calloc(fd->headssize1,sizeof(double));

	fd->pressuressize1 = fd->nNodes+1;
	fd->pressures = (double *)calloc(fd->pressuressize1,sizeof(double));

	free2DMat(connectivityTable, connsize1);
	free2DMat(nodes, nodesize1);

	nNodes = 0;

}/* produceNodes */


/*not currently implemented (or working), kept here though for future use, 
in case someone comes back to it */
void 
adjustPipeMeasurePoints(Geometrydata *gd, Fluidsdata *fd){
	/*This function adjusts the locations of measured points generated
	to monitor T intersections. For instance:
		|                                     |
	   1|2                                   1|2
	----------    will be changed to:  ---------------
	  (3,4)                                  3 4
	  with 3 and 4 moved slightly to the left and right, respectively,
	  so that the geometry for the widths looks more like:

	  1--2                       1--2
	  |  |   and less like        \/
	  3--4                        34                     */

	int pipe_i, pipe_ii, v_i, v_ii;
	int v_ind1, v_ind2,v_pair_1,v_pair_2;
	double xmp1, ymp1, xmp2, ymp2;
	double xp1, yp1, xp2, yp2;
	double x_dummy, y_dummy, m;
	double ** A_proj;
	double * b_proj, * newmpXY;
	
	void gauss(double ** A, double * b, int sizeAMatrix, double * returnVec);

	A_proj = DoubMat2DGetMem(2,2);
	b_proj = (double *)malloc(2*sizeof(double));
	newmpXY = (double *)malloc(2*sizeof(double));
	
	for (pipe_i = 1; pipe_i <= fd->nPSegments; pipe_i++){

		/*looking at point 1 (indices 1 and 2) first */

		for (v_i = 1; v_i <= 4; v_i++){

			v_ind1 = (int)fd->pipeSegments[pipe_i][v_i];
		
 			/*first, find the new measured points
			shouldn't be any more than one per pair*/
			if (v_ind1 < 0){
				/*next, find the point that it matches.
				Do this by comparing x and y coordinates 
				with other measured points */

				xmp1 = gd->points[-v_ind1][1];
				ymp1 = gd->points[-v_ind1][2];

				/*now, comparison with all other negative 
				values in pipes.  Should be able to take the first match 
				Starting from the next pipe in the list after pipe_i*/

				for(pipe_ii = pipe_i+1; pipe_ii <= fd->nPSegments; pipe_ii++){
					for (v_ii = 1; v_ii <= 4; v_ii++){

						v_ind2 = (int)fd->pipeSegments[pipe_ii][v_ii];

						if(v_ind2 < 0){

							xmp2 = gd->points[-v_ind2][1];
							ymp2 = gd->points[-v_ind2][2];

							/*check if they are the same point*/

							if(xmp1 == xmp2 && ymp1 == ymp2){

								/*now, calculate the requisite geometries and offsets*/
								/*first, need indices of the point on the opposite
								side of the pipe from the measured point*/
								if (v_i == 1 || v_i == 3)
									v_pair_1 = (int)fd->pipeSegments[pipe_i][v_i+1];
								else
									v_pair_1 = (int)fd->pipeSegments[pipe_i][v_i-1];

								if (v_ii == 1 || v_ii == 3)
									v_pair_2 = (int)fd->pipeSegments[pipe_ii][v_ii+1];
								else
									v_pair_2 = (int)fd->pipeSegments[pipe_ii][v_ii-1];

								xp1 = gd->vertices[v_pair_1][1];
								yp1 = gd->vertices[v_pair_1][2];
								xp2 = gd->vertices[v_pair_2][1];
								yp2 = gd->vertices[v_pair_2][2];

								/*The goal here is to find the projection of each of the
								two original points onto the line formed by the measured point
								and the slope of the line between the two original points.
								Special cases are the horizontal and vertical lines,
								where point slope form doesn't work.  Should just be 
								straightforward linear algebra. */

								/*first, find the slope of the original lines*/

								if (xp1 == xp2){
									/*vertical line, so ymp = yp*/
									gd->points[-v_ind1][2] = yp1;
									gd->points[-v_ind2][2] = yp2;
								} else if (yp1 == yp2){
									/*horizontal line, so xmp = xp */
									gd->points[-v_ind1][1] = xp1;
									gd->points[-v_ind2][1] = xp2;
								} else {
									/*first, find the slope*/
									m = (yp2-yp1)/(xp2-xp1);
									/*second, use the slope to find another point
									on the line (y-ymp) = m*(x-xmp).  Use the point
									x = 0.0 (as long as xmp != 0.0) */
									if (xmp1 != 0.0)
										x_dummy = 0.0;	
									else
										x_dummy = 1.0;

									y_dummy = m*(x_dummy-xmp1)+ymp1;

									/*third, project the points onto the new line formed,
									and save the result to the original measured points.
									This process uses the Gaussian elimination function
									written earlier.*/
									/*theory for this is taken from
									http://cs.nyu.edu/~yap/classes/visual/03s/hw/h2/math.pdf*/

									/*A_proj will be the same for both point pairs */

									A_proj[0][0] = xmp1-x_dummy;
									A_proj[0][1] = ymp1-y_dummy;
									A_proj[1][0] = y_dummy-ymp1;
									A_proj[1][1] = xmp1-x_dummy;

									/*initialize the return vector*/
									newmpXY[0] = 0.0;
									newmpXY[1] = 0.0;

									/*set b_proj for the first original point */
									b_proj[0] = -(-xp1*(xmp1-x_dummy)-yp1*(ymp1-y_dummy));
									b_proj[1] = -(-y_dummy*(xmp1-x_dummy)+x_dummy*(ymp1-y_dummy));
																

									/*solve the system for the first point pair, and assign the 
									returned values to the xy of the first measured point*/					
									gauss(A_proj, b_proj, 2, newmpXY);
									gd->points[-v_ind1][1] = newmpXY[1];
									gd->points[-v_ind1][2] = newmpXY[2];

									/*now, repeat for the second point */

									/*initialize the return vector*/
									newmpXY[0] = 0.0;
									newmpXY[1] = 0.0;

									/*set b_proj for the first original point */
									b_proj[0] = -(-xp2*(xmp1-x_dummy)-yp2*(ymp1-y_dummy));
									b_proj[1] = -(-y_dummy*(xmp1-x_dummy)+x_dummy*(ymp1-y_dummy));

									gauss(A_proj, b_proj, 2, newmpXY);
									gd->points[-v_ind2][1] = newmpXY[1];
									gd->points[-v_ind2][2] = newmpXY[2];
								}
							}/* if xmp1 == xmp2, ymp1 == ymp2 */
						}/* v_ind2 */
					}/* v_ii */
				}/* pipe_ii */
			}/* if v_ind <0*/
		}/* v_i*/
	}/*pipe_i*/

	free(b_proj);
	/* free(newmpXY); Already freed in gauss */
	free(A_proj[1]);
	free(A_proj[0]);
	free(A_proj);

}/*adjustPipeMeasurePoints*/