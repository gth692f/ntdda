
#ifndef __PIPES_H__
#define __PIPES_H__


#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include "dda.h"

void producePipeSegments(Geometrydata *gd, Fluidsdata *fd);
void produceNodes(Geometrydata *gd, Fluidsdata *fd);
void initializePipeProperties(Geometrydata *gd, Fluidsdata *fd);
void calcPipeWidths(Geometrydata *gd, Fluidsdata *fd);
void adjustPipeMeasurePoints(Geometrydata *gd, Fluidsdata *fd);

#ifdef __cplusplus
}
#endif

#endif /* __PIPES_H__ */
