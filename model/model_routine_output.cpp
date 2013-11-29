/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <math.h>

#include "biocellion.h"

#include "model_routine.h"

/* UESR START */

#include "model_define.h"

/* UESR END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentOutput( const VIdx& vIdx, const SpAgent& spAgent, REAL& color, Vector<REAL>& v_extra ) {
	/* MODEL START */
	
	color = spAgent.state.getType();
	// show agents involved in buding  
	if (  spAgent.junctionInfo.getNumJunctions()==1 ) {
		color += 1;
		//distinguish mother ends of bud
		if ( spAgent.state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 1 ) {
			color += 1;
		}
	}

	CHECK( v_extra.size() == 0 );

	/* MODEL END */

	return;
}
#endif

