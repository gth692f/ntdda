
#ifndef __JOINT_H__
#define __JOINT_H__

#include <ddatypes.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct _joint Joint;
typedef struct _jointmat Jointmat;
typedef struct _jointprop Jointprop;


/** @todo Change to incomplete type. */
struct _joint {
   int	type;
   double epx1, epy1, epx2, epy2;
   double length;
   int is_selected;
};




struct _jointmat {
  /* Need a union of structs and an 
   * enum in here to handle the various types
   * of friction laws.
   */
   double friction;
   double cohesion;
	 double tension;
   int type;
};



struct _jointprop {
  /* Need a union of structs and an 
   * enum in here to handle the various types
   * of friction laws.
   */
   double width;
   double friction_factor;
   int constant_width;
};



Jointmat * jointmat_new           (void);
void       jointmat_delete        (Jointmat *);

Jointmat * jointmat_array_new     (int numjoints);

double     jointmat_get_friction  (Jointmat *);
void       jointmat_set_friction  (Jointmat *,
                                   double friction);


double     jointmat_get_cohesion  (Jointmat *);
void       jointmat_set_cohesion  (Jointmat *,
                                   double cohesion);


double     jointmat_get_tension   (Jointmat *);
void       jointmat_set_tension   (Jointmat *,
                                   double tension);


int        jointmat_get_type      (Jointmat *);
void       jointmat_set_type      (Jointmat *,
                                   int type);

int        joint_get_type         (Joint * j);
void       joint_set_type         (Joint * j,
                                   int type);


void       joint_set_endpoints    (Joint * j,
                                   double epx1,
                                   double epy1,
                                   double epx2,
                                   double epy2);

Jointprop * jointprop_new           (void);
void       jointprop_delete        (Jointprop *);

Jointprop * jointprop_array_new     (int numjoints);

double     jointprop_get_width  (Jointprop *);
void       jointprop_set_width  (Jointprop *,
                                   double width);

double     jointprop_get_friction_factor  (Jointprop *);
void       jointprop_set_friction_factor  (Jointprop *,
                                   double friction_factor);

int     jointprop_get_constant_width  (Jointprop *);
void       jointprop_set_constant_width  (Jointprop *,
                                   int constant_width);


/**
 * The scale factor is new_scale/old_scale.
 */
void       joint_rescale          (Joint * j,
                                   double scale_factor);

void       joint_print_xml        (Joint * j,
                                   PrintFunc printer,
                                   void * stream);




#ifdef __cplusplus
}
#endif

#endif  /* __JOINT_H__ */