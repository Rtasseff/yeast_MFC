/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <math.h>

#include <iostream>

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox, const Vector<double>& v_y, Vector<double>& v_f ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state ) {
	/* MODEL START */

	REAL R0 = state.getRadius();
	REAL V0 = GEN_PI43 * R0 * R0 * R0;




	/* ---Growth--- 
	Assuming exponential growth
	*/
	
	// FD to calc step size
	REAL V1 = V0 + STATE_AND_GRID_TIME_STEP * GROWTH_RATE * V0;
	state.setRadius( CBRT( V1 / GEN_PI43 ) );
	
	


	/* MODEL END */

	return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, BOOL& divide, BOOL& disappear ) {
	/* MODEL START */

	/* ---Divide Decision--- 
	Assuming critical, prefect timer,
	Assuming an all at once divide,
	should improve to account for bud
	*/
	if ( spAgent.state.getModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK ) > CC_CLOCK_CRITICAL ){
		divide = true;
	}
	

	/* MODEL END */

	return;
}

void ModelRoutine::adjustSpAgentState( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state, VReal& disp ) {
	/* MODEL START */
	REAL R0 = state.getRadius();
	REAL V0 = GEN_PI43 * R0 * R0 * R0;

	/* ---Cell Cycle Progression--- 
	Assuming critical size, const theta model
	Charvin2009
	*/
	// check critical V:
	if ( V0 > VOL_CRITICAL ) {
		// move internal clock:
		state.incModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, (CC_CLOCK_RATE * BASELINE_TIME_STEP_DURATION) );
	

	/* MODEL END */

	return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& motherState, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

	REAL R0 = motherState.getRadius();
	REAL V0 = GEN_PI43 * R0 * R0 * R0;

	/* ---Divide--- 
	assuming mother critical volume
	assuming an all at once divide,
	should improve to account for budding
	*/
	daughterState.setType( motherState.getType() );
	CHECK( V0 > VOL_CRITICAL); 
	// set radius assuming all V > VC goes to d
	REAL Vd = V0 - VOL_CRITICAL;
	motherState.setRadius( CBRT( VOL_CRITICAL / GEN_PI43 ) );
	daughterState.setRadius( CBRT( Vd / GEN_PI43 ) );
	// displacement determined by bud dir and radius
	motherDisp[0] = -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) * motherState.getRadius();
	motherDisp[1] = -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) * motherState.getRadius();
	motherDisp[2] = 0.0;
	daughterDisp[0] = motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) * daughterState.getRadius();
	daughterDisp[1] = motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) * daughterState.getRadius();
	daughterDisp[2] = 0.0;

	// set up state
	CHECK( NUM_YEAST_CELL_MODEL_REALS==3 );
	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, -1.0 *  motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ));
	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y, -1.0 *  motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ));
	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK );

	// reset mothers cell cycle clock
	motherState.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, 0.0 )

	/* MODEL END */

	return;
}

void ModelRoutine::disappearSpAgent( const VIdx& vIdx, const SpAgent& spAgent, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<IfGridUpdate>& v_gridDelta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}
#endif

