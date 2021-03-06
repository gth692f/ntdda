
#ifndef __DDAFILE_H__
#define __DDAFILE_H__

#include <stdio.h>


typedef struct _filepath_tag Filepaths;

typedef struct _ddafile DDAFile;

enum inputfiletype { original = 0, extended, ddaml } inputfiletype;


void initFilePaths                  (Filepaths * filepaths);

void openAnalysisFiles              (Filepaths * filepaths);

int  ddafile_get_type               (char * filename);


/* Maybe move this back into adata file later... */
void ddafile_read_original_analysis (void * userdata,
                                     char *, 
                                     int numfixedpoints,
                                     int pointcount,
                                     int numloadpoints);


#define FNAME_BUFSIZE  256


/** Yes yes, this is all very ugly,
 * but it turns out that carrying all 
 * this cruft ariound makes things 
 * easier to deal with for Win32.
 * So it has to stay in for now.
 */
struct _filepath_tag {

char gpath[FNAME_BUFSIZE];
char path[FNAME_BUFSIZE];
char apath[FNAME_BUFSIZE];
char vpath[FNAME_BUFSIZE];
char fpath[FNAME_BUFSIZE];
char oldpath[FNAME_BUFSIZE];
char gfile[FNAME_BUFSIZE];
char afile[FNAME_BUFSIZE];
char ffile[FNAME_BUFSIZE];
char vfile[FNAME_BUFSIZE];
char oldfile[FNAME_BUFSIZE];
char resfile[FNAME_BUFSIZE];
char replayfile[FNAME_BUFSIZE];
char logfile[FNAME_BUFSIZE];
char errorfile[FNAME_BUFSIZE];
char blockfile[FNAME_BUFSIZE];
char measfile[FNAME_BUFSIZE];
char porefile[FNAME_BUFSIZE];
char timefile[FNAME_BUFSIZE];
char parfile[FNAME_BUFSIZE];
char momentfile[FNAME_BUFSIZE];
char massfile[FNAME_BUFSIZE];
char gravfile[FNAME_BUFSIZE];
char htmlfile[FNAME_BUFSIZE];
char datafile[FNAME_BUFSIZE];
char gnuplotfile[FNAME_BUFSIZE];
char spyfile1[FNAME_BUFSIZE];
char spyfile2[FNAME_BUFSIZE];
char mfile[FNAME_BUFSIZE];
char dfile[FNAME_BUFSIZE]; // deformation (solution) vector
char eqfile[FNAME_BUFSIZE];
char fpointfile[FNAME_BUFSIZE];
char headfile[FNAME_BUFSIZE];
char pipefile[FNAME_BUFSIZE];
char cforce[FNAME_BUFSIZE];
char fforce[FNAME_BUFSIZE];
char stressfile[FNAME_BUFSIZE];
char boltfile[FNAME_BUFSIZE];
char boltlogfile[FNAME_BUFSIZE];
char vertexfile[FNAME_BUFSIZE];
char vertexlogfile[FNAME_BUFSIZE];
/* Note that the rootname contains the entire path.  This 
 * is probably not real good, but is messy to handle 
 * otherwise on win32.
 */
char rootname[FNAME_BUFSIZE];
};



typedef struct file_pointers {
   FILE * logfile;
   FILE * replayfile;
   FILE * errorfile;
   FILE * measfile;
   FILE * porefile;
   FILE * timefile;
   FILE * parfile; /* parameters */
   FILE * momentfile; /* block moments */
   FILE * massfile;
   FILE * gravfile;
   FILE * htmlfile;
   FILE * datafile;
   FILE * gnuplotfile;
   FILE * spyfile1;
   FILE * spyfile2;
   FILE * mfile;
   FILE * dfile;
   FILE * eqfile;
   FILE * fpointfile;
   FILE * headfile;
   FILE * pipefile;
   FILE * cforce;
   FILE * fforce;
   FILE * stressfile;
   FILE * boltfile;   // matlab format
   FILE * boltlogfile;  // log format
   FILE * vertexfile;  // matlab format
   FILE * vertexlogfile;  // log format
} FILEPOINTERS;

#endif  /* __DDA_FILE_H__ */

