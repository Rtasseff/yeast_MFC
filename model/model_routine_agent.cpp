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

	if( init == true ) {
		if ( startVIdx[2] == 0 && startVIdx[1] <= 7 && (startVIdx[1]+ regionSize[1]) >= 7 && startVIdx[0] <= 5 && (startVIdx[0]+ regionSize[0]) >= 5 ) {
			VIdx posVIdx;
			VReal posOffset;
			SpAgentState state;
			REAL radius = R_CRITICAL;

			posVIdx[0] = 5;
			posVIdx[1] = 7;
			posVIdx[2] = 0;
			posOffset[0] = 0.0;	
			posOffset[1] = 0.0;	
			posOffset[2] = 0.0;	

			// set up the state
			state.setType( AGENT_YEAST_CELL );
			state.setRadius( radius );
			for ( S32 i = 0; i < NUM_YEAST_CELL_MODEL_REALS; i++ ) {
				state.setModelReal( i, 0 );
			}
			state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, 1 );
			CHECK( ifGridHabitableBoxData.get( posVIdx ) == true );
			v_spAgentVIdx.push_back( posVIdx );
			v_spAgentState.push_back( state );
			v_spAgentOffset.push_back( posOffset );

		}
	}


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

	if (state.getType() == AGENT_YEAST_CELL ) { // yeast cell type -->

		REAL R0 = state.getRadius();
		REAL V0 = GEN_PI43 * R0 * R0 * R0;
		BOOL hasGlucose = false;


		/* ---Nutriant Uptake--- 
		Uptake calculated at grid,
		check if enough was avalible:
		*/
		if ( v_gridModelIntNbrBox[GRID_MODEL_INT_GLUCOSE_AVAILABLE].getVal( 0, 0, 0) == 1 ){
			//cell uptakes the stuff
			hasGlucose = true;
		}

		/* ---Growth--- 
		Assuming exponential growth
		only if glucose is avalible
		*/
		if (hasGlucose == true) {
			// FD to calc step size
			REAL V1 = V0 + STATE_AND_GRID_TIME_STEP * GROWTH_RATE * V0;
			state.setRadius( CBRT( V1 / GEN_PI43 ) );
		}


		/* ---Cell Cycle Progression--- 
		Assuming critical size, const theta model
		Charvin2009
		Only do if we have glucose
		*/
		if (hasGlucose == true) {
			// check critical V:
			if ( V0 > VOL_CRITICAL ) {
				// move internal clock:
				state.incModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, (CC_CLOCK_RATE * STATE_AND_GRID_TIME_STEP) );
			}
		}

	} // <-- Yeast cell type

	


	/* MODEL END */

	return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */


	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, BOOL& divide, BOOL& disappear ) {
	/* MODEL START */
	disappear = false;
	divide = false;
	if (spAgent.state.getType() == AGENT_YEAST_CELL ) { // yeast cell type -->
		/* Dissapper Decision
		assume that if cell touches boundry on 
		x limits, it gets carried away by flow
		*/
		REAL r = spAgent.state.getRadius();
		REAL dist = IF_GRID_SPACING*10.0;
		// low wall
		if ( vIdx[0]==0 ){
			dist = spAgent.vOffset[0] + 0.5 * IF_GRID_SPACING;
		}
		if ( vIdx[0] == (UB_NUM[0] - 1) ){
			dist = 0.5 * IF_GRID_SPACING - spAgent.vOffset[0];
		}

		if ( r > dist ){
			disappear = true;
		
		/* ---Divide Decision--- 
		Assuming critical, prefect timer,
		Assuming an all at once divide,
		should improve to account for bud
		dont divide if you disappear
		*/
		}
		else if ( spAgent.state.getModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK ) > CC_CLOCK_CRITICAL ){
			divide = true;
		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& vOffset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& state/* INOUT */, VReal& disp ) {
	/* MODEL START */
	REAL R0 = state.getRadius();
	REAL V0 = GEN_PI43 * R0 * R0 * R0;
	VReal force = mechIntrctData.force;


	
	/* ---Cell Boundry Collision---
	The cell boundry collision will consist of shoving.
	We must specificly look for boundry with chip desing (dummy spaces).
	We will not look at boundry of x high or low domain.
	*/

	REAL dist;
	VIdx vIdxTmp;


	// impact with y bound low
	if ( vIdx[1]==0 ){
		dist = vOffset[1] + 0.5 * IF_GRID_SPACING;
		if ( dist < R0 ) {
			// hitting low bound, push up
			force[1] += CELL_WALL_STIFF * (R0 - dist);
		}

	}
	// impact with y bound high
	if ( vIdx[1] == (UB_NUM[1] - 1) ){
		dist = 0.5 * IF_GRID_SPACING - vOffset[1];
		if (dist < R0) {
			//hitting low bound, push down 
			force[1] += CELL_WALL_STIFF * (dist - R0);
		}
	}
	// impact with chip wall
	for( S32 dim = 0 ; dim < 2 ; dim++ ) {
		if ( vIdx[dim]>0 ){
			// check below
			dist = vOffset[dim] + 0.5 * IF_GRID_SPACING;
			vIdxTmp = vIdx;
			vIdxTmp[dim] = vIdxTmp[dim]-1;
			// if within range and its not habitable...
			if ( ( CHIP_DESIGN_MATRIX[vIdxTmp[0]][vIdxTmp[1]] == 0 ) && ( dist < R0 ) ){
				// hitting low bound, push up
				force[dim] += CELL_WALL_STIFF * (R0 - dist);
			}
		}
		if ( vIdx[dim]<(UB_NUM[dim] - 1) ){
			// check above
			dist = 0.5 * IF_GRID_SPACING - vOffset[dim];
			vIdxTmp = vIdx;
			vIdxTmp[dim] = vIdxTmp[dim]+1;
			// if within range and its not habitable...
			if ( ( CHIP_DESIGN_MATRIX[vIdxTmp[0]][vIdxTmp[1]] == 0 ) && ( dist < R0 ) ){
				// hitting upper bound, push down
				force[dim] += CELL_WALL_STIFF * (dist - R0);
			}
		}
	}











	/* ---DISP--- 
	Caclulating displacement based on 
	baseline time step and
	overdamped kinetics 
	*/
	for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
		// dx/dt = 1/c * sum(Force); delta_x = dx/dt*delta_t
		disp[dim] = ( 1.0 / (CELL_DAMP_COEF) ) * force[dim] * BASELINE_TIME_STEP_DURATION  ;
		if ( disp[dim] > MAX_DISP ) {
			disp[dim] = MAX_DISP;
			cout <<"WARNING:MRA_ASAS_0001 - Cell moving too much" <<endl;
		}
		else if ( disp[dim] < (-1.0 * MAX_DISP) ) {
			disp[dim] = -1.0 * MAX_DISP;
			cout <<"WARNING:MRA_ASAS_0001 - Cell moving too much" <<endl;
		}
	}




	

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
	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, 0.0 );

	/* Random budding dir in the x,y plane */
	REAL scale;
	VReal dir;
	do {
		scale = 0.0;
		for( S32 dim = 0 ; dim < 2 ; dim++ ) {
			dir[dim] = Util::getModelRand( MODEL_RNG_UNIFORM ) - 0.5;
			scale += dir[dim] * dir[dim];
		}
		scale = SQRT( scale );
	} while( scale < GEN_EPS );

	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, dir[0]/scale );
	daughterState.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y, dir[1]/scale );
	




	// reset mothers cell cycle clock
	motherState.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, 0.0 );

	/* MODEL END */

	return;
}

void ModelRoutine::disappearSpAgent( const VIdx& vIdx, const SpAgent& spAgent, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<IfGridUpdate>& v_gridDelta ) {
	/* MODEL START */

	// do nothing	

	/* MODEL END */

	return;
}
#endif

