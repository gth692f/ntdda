
/*
File with the definitions, reading, etc. for the fluidsddaml
Currently copied over from analysisddaml.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "ddaml_private.h"
#include "ddaml.h"
#include "ddamemory.h"
#include "fluidsdata.h"
#include "timehistory.h"

/** Get rid of this later. 
 * The way to get rid of it is to define 
 * ddaml_display_* functions, then set them
 * to the dda_display_* functions.
 */
#include "dda.h"
#include "joint.h"


/* Mostly for debugging windows. */
static char mess[128];


typedef DList JOINTPROPLIST;
typedef DList BOUNDARYCONDLIST;

static JOINTPROPLIST * jointproplist;
static BOUNDARYCONDLIST * boundarycondlist;


static xmlNsPtr nspace;

void transferFData(void);

void transferJointProplistToAStruct(Fluidsdata * fdata, JOINTPROPLIST * jointmatlist);
void transferBoundaryCondlistToAStruct(Fluidsdata * fdata, BOUNDARYCONDLIST * boundarycondlist); 

static double ** DoubMat2DGetMem(int n, int m);

static void initializeFLists(void);

Fluidsdata * fdata;



void 
initializeFLists()
{
   jointproplist = dlist_new();
   boundarycondlist = dlist_new();

}  /* close initializeLists() */


void 
transferFData(void)
{

	transferJointProplistToAStruct(fdata,jointproplist);
	transferBoundaryCondlistToAStruct(fdata,boundarycondlist);

}  /* close transferFData() */


void 
transferJointProplistToAStruct(Fluidsdata * fd, JOINTPROPLIST * jointproplist)
{
   int i = 0;
   int njprop;
   JOINTPROPLIST * ptr;
   Jointprop * jpropmp;

  /* WARNING: Type attribute of xml tag is ignored for now.
   * Joint properties should be listed in order of occurrence
   * in the xml file.
   */
   njprop = dlist_length(jointproplist);
   fd->nFluidJoints = njprop;
   fd->jointFluidPropsize1 = njprop + 1;
   fd->jointFluidPropsize2 = 3; //Width, friction factor, boolean for if width stays constant
   fd->jointFluidProp = DoubMat2DGetMem(fdata->jointFluidPropsize1,fdata->jointFluidPropsize2);

  /* WARNING: This assumes that joint types are listed in order 
   * of occurrence in xml file.  The type attribute is ignored
   * for now.
   */
   dlist_traverse(ptr, jointproplist)
   {
      jpropmp = ptr->val;
      fd->jointFluidProp[i+1][0] = jointprop_get_width(jpropmp);/*->width (originally jointmat_get_cohesion, find this to change it)*/
      fd->jointFluidProp[i+1][1] = jointprop_get_friction_factor(jpropmp);/*->friction factor; */
	  fd->jointFluidProp[i+1][2] = jointprop_get_constant_width(jpropmp);/* constant width boolean*/
      i++;
   }
    
}  /*  close transferJointPropListToAStruct() */

void
transferBoundaryCondlistToAStruct(Fluidsdata * fd, BOUNDARYCONDLIST *boundarycondlist){

   int i = 0;
   int nboundcond;
   BOUNDARYCONDLIST * ptr;
   BoundaryCond * boundcondmp;

  /* WARNING: Type attribute of xml tag is ignored for now.
   * Joint properties should be listed in order of occurrence
   * in the xml file.
   */
   nboundcond = dlist_length(boundarycondlist);
   fd->nBoundaryCond = nboundcond;
   fd->boundsize1 = nboundcond + 1;
   fd->boundsize2 = 5; 
   fd->boundaryConditions = DoubMat2DGetMem(fdata->boundsize1,fdata->boundsize2);

   dlist_traverse(ptr, boundarycondlist)
   {
      boundcondmp = ptr->val;
	  fd->boundaryConditions[i+1][1] = boundarycond_get_xcoord(boundcondmp); 
      fd->boundaryConditions[i+1][2] = boundarycond_get_ycoord(boundcondmp);  
	  fd->boundaryConditions[i+1][3] = boundarycond_get_condition(boundcondmp); 
	  fd->boundaryConditions[i+1][4] = boundarycond_get_type(boundcondmp);
      i++;
   }
    
}  /*  close transferBoundaryCondlistToAStruct() */



void 
parseJointFluidProperties(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{
   Jointprop * jprop;
   double temp;   
   jprop = jointprop_new(); // was jointmat_new();

   if (!strcmp("yes",xmlGetProp(cur,"ConstantWidth"))) {
	   /*1 for constant condition, 0 for variable*/
		jointprop_set_constant_width(jprop,1);
   } 
   else if (!strcmp("no",xmlGetProp(cur,"ConstantWidth"))) {
	   /*1 for constant condition, 0 for variable*/
	   jointprop_set_constant_width(jprop,0);
   }

   cur = cur->children;

   /*WEM note, joints must be input in order, type field not used
   (maybe starting from one, going up (hasn't been tested)*/

   while (cur != NULL) {

      if ((!strcmp(cur->name, "Width")) ) {
         temp = atof(xmlNodeListGetString(doc, cur->children, 1));
         jointprop_set_width(jprop,temp);
      }

      if ((!strcmp(cur->name, "FrictionFactor")) ) {
         temp = atof(xmlNodeListGetString(doc, cur->children, 1));
         jointprop_set_friction_factor(jprop,temp);
      }

      cur = cur->next;
   }

   dl_insert_b(jointproplist,(void*)jprop);


}  /*close parseJointFluidproperties() */

void 
parseBoundaryCondition(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{
   BoundaryCond * bcond;
   double temp;
   bcond = BoundaryCond_new(); 

   if (!strcmp("flow",xmlGetProp(cur,"type"))) {
	   /*1 for flow condition, 2 for head*/
		boundarycond_set_type(bcond,1);
   } 
   else if (!strcmp("head",xmlGetProp(cur,"type"))) {
	   /*1 for flow condition, 2 for head*/
	   boundarycond_set_type(bcond,2);
   }

   cur = cur->children;

   while (cur != NULL) {

      if ((!strcmp(cur->name, "XCoord")) ) {
         temp = atof(xmlNodeListGetString(doc, cur->children, 1));
         boundarycond_set_xcoord(bcond,temp); 
      }

      if ((!strcmp(cur->name, "YCoord")) ) {
         temp = atof(xmlNodeListGetString(doc, cur->children, 1));
         boundarycond_set_ycoord(bcond,temp); 
      }

	  if ((!strcmp(cur->name, "Condition")) ) {
         temp = atof(xmlNodeListGetString(doc, cur->children, 1));
         boundarycond_set_condition(bcond,temp);
      }

      cur = cur->next;
   }

   dl_insert_b(boundarycondlist,(void*)bcond);


}  /*close parseBoundaryCondition */

void 
parsePipeKSolver(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{

   if (!strcmp("laminar",xmlGetProp(cur,"type"))) {
	   /*1 for flow condition, 2 for head*/
	   fdata->ksolverflag = 1;
   } 
   else if (!strcmp("cubic",xmlGetProp(cur,"type"))) {
	   /*1 for flow condition, 2 for head*/
	   fdata->ksolverflag = 2;
   }
}/*parse pipe solver*/

void 
parseGravityFlag(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{

   if (!strcmp("on",xmlGetProp(cur,"flag"))) {
	   /*1 for flow condition, 2 for head*/
	   fdata->gravityflag = 1;
   } 
   else if (!strcmp("off",xmlGetProp(cur,"flag"))) {
	   /*1 for flow condition, 2 for head*/
	   fdata->gravityflag = 0;
   }
}/*parse pipe solver*/

void 
parseFluidViscosity(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{
   fdata->fluidViscosity = atof(xmlGetProp(cur,"viscosity"));

}  /*close parseFluidViscosity() */


void 
parseFluidDensity(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) 
{
   fdata->fluidDensity = atof(xmlGetProp(cur,"density"));

}  /*close parseFluidDensity() */


static double **
DoubMat2DGetMem(int n, int m)
{
   int i;
   double **x;

   assert ( (m!=0) && (n!=0) );

   //x = (double **)malloc(sizeof(double *)*n);
   x = (double **)calloc(n,sizeof(double *));
   if (x == NULL)
      return NULL;

   for ( i = 0; i < n; ++i)
   {
      //x[i] = (double *)malloc(sizeof(double)*m);
      x[i] = (double *)calloc(m,sizeof(double));
      if(x[i] == NULL)
         return NULL;
      //else 
         //memset(x[i], 0xDDA, m);
   }

   return x;

}  


/* this block, originally #if 0'd out*/
#if 0
static JointMat *
getNewJointMat(void)
{
   JointMat * jm;
   //fprintf(stdout,"Getting new joint\n");
   jm = (JointMat *)calloc(1,sizeof(JointMat));
   return jm;

}  /*  close getNewJointMat() */
#endif

/*this one, not originally, making it so now*/

#if 0
static BlockMat *
getNewBlockMat(void)
{
   BlockMat * bm;
   //fprintf(stdout,"Getting new block mat\n");
   bm = (BlockMat *)malloc(sizeof(BlockMat));
   memset(bm,0xDA,sizeof(BlockMat));
   return bm;

}  /*  close getNewJointMat() */
#endif

KWDTAB ftab_type[] = {
{"JointFluidProperties", 0, *parseJointFluidProperties },
{"BoundaryCondition",    0, *parseBoundaryCondition    },
{"FluidViscosity",       0, *parseFluidViscosity       },
{"FluidDensity",         0, *parseFluidDensity         },
{"PipeKSolver",          0, *parsePipeKSolver          },
{"GravityFlag",          0, *parseGravityFlag          },
{NULL ,                  0,  0                         }			
};

/*More or less stopped here
TODO:
Write functions for some of the things that I just arbitrarily changed up top
*/


void
parseFluids(Fluidsdata * fd, xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) {

   int i = 0;

  /* Easiest to parse into a list, then count 
   * and grab array memory, then copy data into 
   * arrays.  This will be good, because I will 
   * now have lists of elements that can eventually
   * be used in the main code.
   */
   initializeFLists(); 

   while (cur != NULL) 
   {
      i = 0;
      while (ftab_type[i].kwd)
      {
         if ((!strcmp(ftab_type[i].kwd,cur->name)) ) //&& (cur->ns == ns)) 
         {
            fprintf(stderr,"A Keyscan loop: Current name %s:\n", cur->name);
            switch(ftab_type[i].ktok)
            {
               case node:
	                 ftab_type[i].parsefn.nodeparse(doc,ns,cur);
               break;

               case string:
	                 ftab_type[i].parsefn.stringparse(doc,cur,1);
               break;

               case prop:
	                 ftab_type[i].parsefn.propparse(doc,cur,1);
               break;

               default:;
               break;
            }  /* end switch on node type */
         }
         i++;
      }   //close key scan loop 

      cur = cur->next;

   }  /* end loop over current pointer */

  /** Transfer data from list format to array format. */
   transferFData();

}  


void
ddaml_read_fluids_file(Fluidsdata * fd, char * filename) {
   
   xmlDocPtr doc;
   xmlNsPtr ns;
   xmlNodePtr cur;

   xmlNode *root_element = NULL;


   //xmlDoValidityCheckingDefaultValue = 1;

   fdata = fd;

   //ddaml_display_error   = fd->display_error;

  /**
   * build an XML tree from the file;
   */
   doc = xmlParseFile(filename);


   //ddaml_check_document((void*)doc,"http://www.tsoft.com/~bdoolin/dda","DDA");
   //cur = doc->root;

   root_element = xmlDocGetRootElement(doc);
   cur = root_element;
   ns = nspace;

  /** Now, walk the tree.  All the following code will
   * be replaced with keyscanning.  Initially, do not
   * check for keyscanning failure.  Let the dtd validator
   * handle that.  Later, a few simple checks can be added,
   * just to handle exceptional cases.
   */
  /* FIXME:  This needs to be handled as an exception to
   * protect against dereferencing a null pointer.
   */

   if(cur == NULL){
	   dda_display_warning("Error in opening fluids file, using arbitrary values");
	   return ;
   };

   cur = cur->children;

   while (cur != NULL) {
         if ( !strcmp( cur->name, "Fluids") )  {
            cur = cur->children;
            parseFluids(fd,doc, ns, cur); 
            break;
         }
	     cur = cur->next;
   }

   assert(fd != NULL);
} 




#ifdef STANDALONE
int 
main(int argc, char **argv) 
{
   int i;
   fprintf(stderr,"Parsing Fluids data...\n");

   for (i = 1; i < argc ; i++) {
	     fdata = XMLparseDDA_Fluids_File(argv[i]);
   }

   return(0);
}  

#endif /* STANDALONE */

