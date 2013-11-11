#ifndef __MY_DEFINE_H__
#define __MY_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

/* MODEL START */



/* ---Agents--- */

typedef enum _agent_type_e {
	AGENT_YEAST_CELL,
	NUM_AGENT_TYPES
} agent_type_e;

/* ---Difusable elements--- */

typedef enum _diffusible_elem_e {
	DIFFUSIBLE_ELEM0,
	NUM_DIFFUSIBLE_ELEMS
} diffusible_elem_e;

/* ---Grid Properties--- */

typedef enum _grid_model_real_e {
	GRID_MODEL_REAL_FV_RATIO,
	NUM_GRID_MODEL_REALS
} grid_model_real_e;

/* ---YEAST_CELL properties--- */

typedef enum _yeast_cell_model_real_e {
	YEAST_CELL_MODEL_REAL_BUD_DIR_X, // = 1, x dir for bud
	YEAST_CELL_MODEL_REAL_BUD_DIR_Y, // = 0, y dir for bud
	YEAST_CELL_MODEL_REAL_CC_CLOCK, // = 0, current cell cycle pos
	NUM_YEAST_CELL_MODEL_REALS
} yeast_cell_model_real_e;


/* ---General--- */
// pi
const REAL GEN_PI = 3.14159265359;
// 4/3 pi
const REAL GEN_PI43 = (4.0/3.0)*GEN_PI;

/* ---Grid Properties--- */
const S32 NUM_AMR_LEVELS = 2;

/* ---Cell Properties--- */
/* -physical properties- */

/* -Growth and Division- */
// dv/dt = GROWTH_RATE * V, vol double at 84min
const REAL GROWTH_RATE = 0.0.0001375; // [=] sec^-1, cite Charvin2009 
// linear rate of cell cycle clock, 
const REAL CC_CLOCK_RATE = 0.0002347; // [=] sec^-1, cite Charvin2009 
/* Critical volume maitained by cell, also strts CC clock
Based here on radius of 2 micrometers, roughly average size of yeast
*/
const REAL R_CRITICAL = 2.0; // [=] micro meters
const REAL VOL_CRITICAL = R_CRITICAL * R_CRITICAL * R_CRITICAL * GEN_PI43; // [=] um^3
// Reset limit on clock 
const REAL CC_CLOCK_CRITICAL = 1.0; // cite Charvin2009



/* ---Domain--- */
const REAL IF_GRID_SPACING = R_CRITICAL * 2.0;/* this should be equal to or larger than MAX_CELL_RADIUS * 2.0, domain size in the xml file = 128 * 128 * 4928 */
const REAL BASELINE_TIME_STEP_DURATION = 100; // sec
const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 10;
const REAL STATE_AND_GRID_TIME_STEP = BASELINE_TIME_STEP_DURATION / NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;



/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

