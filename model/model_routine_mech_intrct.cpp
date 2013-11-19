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


	/* MODEL END */

	return;
}
#endif

