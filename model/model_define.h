#ifndef __MY_DEFINE_H__
#define __MY_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

/* MODEL START */

typedef enum _diffusible_elem_e {
	DIFFUSIBLE_ELEM0,
	NUM_DIFFUSIBLE_ELEMS
} diffusible_elem_e;

typedef enum _grid_model_real_e {
	GRID_MODEL_REAL_FV_RATIO,
	NUM_GRID_MODEL_REALS
} grid_model_real_e;

const S32 NUM_AMR_LEVELS = 2;

/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

