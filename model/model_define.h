#ifndef __MY_DEFINE_H__
#define __MY_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

/* MODEL START */

/* --- UNITS:
time = seconds (sec);
length = micrometers ($\mu m$);
volume = cubic micrometers ($\mu m^3$);
mass = pico grams (pg);
concentration = $frac\{ng}{\mu m^3}$;
*/

/* ---Agents--- */

typedef enum _agent_type_e {
	AGENT_YEAST_CELL,
	NUM_AGENT_TYPES
} agent_type_e;


/* ---YEAST_CELL properties--- */


typedef enum _yeast_cell_model_real_e {
	YEAST_CELL_MODEL_REAL_BUD_DIR_X, // = 1, x dir for bud
	YEAST_CELL_MODEL_REAL_BUD_DIR_Y, // = 0, y dir for bud
	YEAST_CELL_MODEL_REAL_CC_CLOCK, // = 0, current cell cycle pos
	YEAST_CELL_MODEL_REAL_MOTHER_VOL, // size of the mother at bud formation, 0 if not a bud.
	NUM_YEAST_CELL_MODEL_REALS
} yeast_cell_model_real_e;

typedef enum _yeast_cell_model_int_e {
	YEAST_CELL_MODEL_INT_MOTHER, // = 0, tracks if cell is a mother ie value of 1
	NUM_YEAST_CELL_MODEL_INTS
} yeast_cell_model_int_e;


/* Integer ExtraMech properties for all cell types */
typedef enum _extra_mech_yeast_cell_model_real_e {
	EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_X,  // addditonal force in x direction
	EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_Y,  // addditonal force in x direction
	EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X,  // updated direction of bud
	EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_Y,  // updated direction of bud
	NUM_EXTRA_MECH_YEAST_CELL_MODEL_REALS
} extra_mech_yeast_cell_model_real_e;

/* ---Aditional Output---
the size must agree with run_param.xml
*/
typedef enum _output_model_real_e {
	OUTPUT_MODEL_REAL_BUD_DIR_X,
	OUTPUT_MODEL_REAL_BUD_DIR_Y,
	OUTPUT_MODEL_REAL_ID,
	NUM_OUTPUT_MODEL_REALS
} output_model_real_e;


/* ---Junction End--- */
typedef enum _junction_end_type_e {
	JUNCTION_END_BUD,
	NUM_JUNCTION_END_TYPES
} junction_end_type_e;


typedef enum _junction_end_model_int_e {
	NUM_JUNCTION_END_MODEL_INTS
} junction_end_model_int_e;

/* ---Difusable elements--- */

typedef enum _diffusible_elem_e {
	DIFFUSIBLE_ELEM_GLUCOSE,
	NUM_DIFFUSIBLE_ELEMS
} diffusible_elem_e;

/* ---Grid Properties--- */

typedef enum _grid_model_real_e {
	GRID_MODEL_REAL_GLUCOSE_DELTA, // amount of glucose change in grid
	GRID_MODEL_REAL_GLUCOSE_FRAC_AVAILABLE, // amount of glucose avalible / needed by agents
	GRID_MODEL_REAL_AGENT_VOL,  // grid volume taken up by agents.
	GRID_MODEL_REAL_AGENT_VOL_test, // volume of of all agents in this grid.
	NUM_GRID_MODEL_REALS
} grid_model_real_e;

typedef enum _grid_model_int_e {
	NUM_GRID_MODEL_INTS
} grid_model_int_e;

typedef enum _model_rng_type_e {
	MODEL_RNG_UNIFORM,
	NUM_MODEL_RNGS
} model_rng_type_e;



/* ---General--- */
// pi
const REAL GEN_PI = 3.14159265359;
// 4/3 pi
const REAL GEN_PI43 = (4.0/3.0)*GEN_PI;
const REAL GEN_SMALL = 1.0E-10;
const REAL GEN_EPS = 1E-52;
const BOOL WRITE_WARNING = true; // set to one to write modeling warnings, set to 0 to ignore.



/* ---Domain--- */
const REAL IF_GRID_SPACING = 6;/* this should be equal to or larger than MAX_CELL_RADIUS * 2.0, domain size in the xml file = 128 * 128 * 4928 */
const REAL BASELINE_TIME_STEP_DURATION = 1; // sec
/* grid steps for balancing diffusion and uptake,
was estimated by how much glucose needed for 4 cells in full packed box, kappa = 0.1
and the effective concentration equal to BC, ie not diffusion limited.
*/
const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1;
const REAL STATE_AND_GRID_TIME_STEP = BASELINE_TIME_STEP_DURATION / ( REAL ) NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;
// maximum displacement per step
const REAL MAX_DISP = IF_GRID_SPACING; 
/* -Grid Properties- */
const S32 NUM_AMR_LEVELS = 1;
const S32 NUM_PDE_TIME_STEPS_PER_STATE_AND_GRID_STEP = 1;
const S32 OVLP_MAX_LEVEL = 4; //number of max leves for determining agent overlap



/* ---Diffusion--- */
// concentration of elements in the bulk fluid in the flow channel
const REAL ELEM_BULK_CONCENTRATION[NUM_DIFFUSIBLE_ELEMS] = {2.0e-2}; // pg/um^3
const REAL ELEM_BETA[NUM_DIFFUSIBLE_ELEMS] = {600}; //um^2/sec
const REAL KAPPA_MIN = 0.1; // minimum kappa
const REAL BETA_MIN_SCALE = 0.05;
const REAL MAX_UPTAKE_FRAC = 0.9; //maximum uptake in UB = total*MAX_UPTAKE_FRAC

/* ---Cell Properties--- */

/* -Growth and Division- */
// dv/dt = GROWTH_RATE * V, vol double at 84min
const REAL GROWTH_RATE = 0.0001375; // [=] sec^-1, cite Charvin2009 
// linear rate of cell cycle clock, 
const REAL CC_CLOCK_RATE = 0.0002347; // [=] sec^-1, cite Charvin2009
// ammount of mass given from mother to daughter, should be less then min mass gained in g1
const REAL VOL_MOTHER_BUD = 2; // [=] um^3

const REAL BUD_DIR_SCALE = 0.025; // scale factor for random perturbation of next bud location
const REAL BUD_OVERLAP = 1; // um, length of overlap allowed for 2 agents that are connected by bud
/* Critical volume maitained by cell, also strts CC clock
Based here on radius of 2 micrometers, roughly average size of yeast
*/
const REAL R_CRITICAL = 2.0; // [=] micro meters
const REAL VOL_CRITICAL = R_CRITICAL * R_CRITICAL * R_CRITICAL * GEN_PI43; // [=] um^3
// Reset limit on clock 
const REAL CC_CLOCK_G1 = .25; // cite Charvin2009, end of g1 and start of bud formaiton 
const REAL CC_CLOCK_CRITICAL = 1.0; // cite Charvin2009
/* -physical properties- */
// maximum radius we can have for yeast cell, any bigger the sim crashes:
const REAL CELL_R_MAX = IF_GRID_SPACING/2.; // maximum allowed for IF dist     dep:2.4 // 1.1*R_CRITICAL; //1.2599*R_CRITICAL;
const REAL CELL_VOL_MAX = CELL_R_MAX * CELL_R_MAX * CELL_R_MAX * GEN_PI43;
// maximum radius of yeast mother at budding, prevents yeast cell from growing past max radius in next G1 at cost of increased bud size
const REAL CELL_R_BUD_MAX = CELL_R_MAX - 0.5 ; // maximum allowed for IF dist     dep:2.4 // 1.1*R_CRITICAL; //1.2599*R_CRITICAL;
const REAL CELL_VOL_BUD_MAX = CELL_R_BUD_MAX * CELL_R_MAX * CELL_R_MAX * GEN_PI43;
// Maximum interaction distance for 2 cells
const REAL CELL_INTRCT_DIST_MAX = IF_GRID_SPACING;// min geo size  dep:2.0*CELL_R_MAX;
// standard uptake of difusable elements 

/* -links to enviornment- */
const REAL CELL_ELEM_CONSTANT_UPTAKE[NUM_DIFFUSIBLE_ELEMS] = {4.17e-2}; // pg/(sec*cell)


/* ---Cell Properties 2--- */
/* -mechanics- */
// friction coefficent 
const REAL CELL_DAMP_COEF = 1.0;
/* cell stiffness for shoving currently set to disp cells 
by 1/2 of total disp per step based on single shoving interaction
*/
const REAL CELL_STIFF = 0.5 * ( CELL_DAMP_COEF / BASELINE_TIME_STEP_DURATION );
/* cell wall combine stiffness for shoving when yeast hits wall
assuming cell is 2 orders more stiff then wall, 
geometric mixture -> one order less
*/
const REAL CELL_WALL_STIFF = CELL_STIFF;// / 2.0;
const REAL BUD_STIFF = CELL_STIFF;




/* ---Chip Design---
A large matrix identifying each UB 
in the simulation as habitable or not
assuing flow is at the low and high
x-axis points.
Trapping Chip
11 in x, 15 in Y
5 by 5 block uninhabitable in each corrner.
*/
const REAL ADD_WALL = 0; // thinkness set by IF  dep:( IF_GRID_SPACING / 2.0 ) - R_CRITICAL ;//0.5; // additional thickenss to wall in um, make trap size of typical cell 
const S32 UB_NUM[2] = {16,24};
const S32 CHIP_DESIGN_MATRIX[16][24] =		       {{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
							{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}};



/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

