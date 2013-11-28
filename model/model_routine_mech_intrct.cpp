/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <float.h>
#include <math.h>

#include <iostream>

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::initJunctionSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const VIdx& vIdx1, const SpAgent& spAgent1, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */
) {
	/* MODEL START */

	link = false;

	/* MODEL END */

	return;
}

void ModelRoutine::computeForceSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const VIdx& vIdx1, const SpAgent& spAgent1, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, VReal& force/* force on spAgent0 due to interaction with spAgent1 (force on spAgent1 due to interaction with spAgent0 has same magnitude but the opposite direction), if force has the same direction with dir, two cells push each other, if has the opposite direction, two cells pull each other. */ ) {
	/* MODEL START */

	
	REAL R = spAgent0.state.getRadius() + spAgent1.state.getRadius();
	REAL mag;/* + for repulsive force, - for adhesive force */

	if( dist <= R ) {/* shoving to remove the overlap */
		mag = CELL_STIFF * ( R - dist );
	}
	else {/* adhesion */
		mag = 0.0;/* no adhesion */
	}

	for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
		force[dim] = mag * dir[dim];
	}


	/* MODEL END */

	return;
}

void ModelRoutine::computeExtraMechIntrctSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const VIdx& vIdx1, const SpAgent& spAgent1, const VReal& dir/* unit direction vector from spAgent1 to spAgent0 */, const REAL& dist, ExtraMechIntrctData& extraMechIntrctData0, ExtraMechIntrctData& extraMechIntrctData1, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */, BOOL& unlink ) {
	/* MODEL START */

	/* ---Bud--- 
	check to see if these two agents are bound by a bud:
	*/
	if ( spAgent0.junctionInfo.isLinked( spAgent1.junctionInfo ) == true ) {
		//cout <<"in junction" <<endl;
		/* check if junction should be broken,
		mother variable will be reset when appropriate,
		and here is the only chance we have to disolve the bond
		*/
		if ( spAgent0.state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 0 && spAgent1.state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 0 ) {
			unlink = true;
		}
		else {
			// caclulate harmonic potential
			REAL R = spAgent0.state.getRadius() + spAgent1.state.getRadius();
			REAL mag = 0;/* + for repulsive force, - for adhesive force */
			// remove original force for the interaction between these
			if( dist <= R ) {/* shoving to remove the overlap */
				mag += -1.0 * CELL_STIFF * ( R - dist );
			}
			// add new force allowing for some overlap with this bud

			mag += BUD_STIFF * ( ( R - BUD_OVERLAP ) - dist );

			extraMechIntrctData0.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_X, mag * dir[0] );
			extraMechIntrctData0.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_Y, mag * dir[1] );
			extraMechIntrctData1.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_X, -1 * mag * dir[0] );
			extraMechIntrctData1.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_Y, -1 * mag * dir[1] );

			// do the rotation, no actual forces involed at this time, simply tracking of direction
			extraMechIntrctData0.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X, -1 * dir[0] );
			extraMechIntrctData0.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_Y, -1 * dir[1] );
			extraMechIntrctData1.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X, dir[0] );
			extraMechIntrctData1.setModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_Y, dir[1] );

		}



	}


	/* MODEL END */

	return;
}
#endif

