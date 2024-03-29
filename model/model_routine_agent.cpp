/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <math.h>

#include <iostream>

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

/* Extra Utility Functions */

static void getIsHabitableArray( const VIdx vIdx, BOOL aa_isHabitable[3][3] );

static void getIsHabitable( const VIdx vIdx, BOOL& isHabitable );

static REAL computeGridDelta( const REAL dir[2], const VReal vOffset, const REAL r, const REAL IF );


#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {
	/* MODEL START */

	if( init == true ) {
		if ( startVIdx[2] == 0 && startVIdx[1] <= 8 && (startVIdx[1]+ regionSize[1]) > 8 && startVIdx[0] <= 5 && (startVIdx[0]+ regionSize[0]) >= 5 ) {


			VIdx posVIdx;
			VReal posOffset;
			SpAgentState state;
			REAL radius = R_CRITICAL;

			posVIdx[0] = 5;
			posVIdx[1] = 8;
			posVIdx[2] = 0;
			posOffset[0] = 0.0;	
			posOffset[1] = 0.0;	
			posOffset[2] = 0.0;	

			// set up the state
			state.setType( AGENT_YEAST_CELL );
			state.setRadius( radius );
			CHECK ( NUM_YEAST_CELL_MODEL_REALS == 4 );
			for ( S32 i = 0; i < NUM_YEAST_CELL_MODEL_REALS; i++ ) {
				state.setModelReal( i, 0 );
			}
			// set a direction for budding 
			state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, -0.7071 );
			state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y, -0.7071 );

			CHECK ( NUM_YEAST_CELL_MODEL_INTS == 1 );
			for ( S32 i = 0; i < NUM_YEAST_CELL_MODEL_INTS; i++ ) {
				state.setModelInt( i, 0 );
			}

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
		REAL deltaV = 0.0;


		/* ---Nutriant Uptake--- 
		Uptake calculated at grid,
		check if enough or how much was avalible
		Since the cell uptakes from grids proportional 
		to its occupation fraction we need to use this 
		to do a weighted sum of avalibility.
		Takes from all grids its in.
		*/
		REAL fracGlucoseAvail = 0;
		REAL frac;
		REAL fracSum = 0;
		S32 xOffset;
		S32 yOffset;
		VIdx ubVIdxOffset;
		BOOL aa_isHabitable[3][3] = {{false,false,false},{false,false,false},{false,false,false}};
		getIsHabitableArray( vIdx, aa_isHabitable );


		/* ---Avalible Glucose--- 
		Here we calcualte the glucose available to this agent
		Since we have already calculated the fraction availible in each grid (based on all agents)
		and we have a fraction needed for this agent and all grids,
		we must calcualte a weighted fraction of glucose avalible.
		*/

		// loop though all the neighbor boxes
		for ( S32 x_i = 0; x_i < 3; x_i++ ) {
			for ( S32 y_i = 0; y_i < 3; y_i++ ) {
				xOffset = x_i - 1;
				yOffset = y_i - 1;
								
				if ( aa_isHabitable[x_i][y_i]==true ){
					ubVIdxOffset[0] = xOffset;
					ubVIdxOffset[1] = yOffset;
					ubVIdxOffset[2] = 0;

					// calculate fraction of overlap for this box 
					frac = Util::computeSphereUBVolOvlpRatio( OVLP_MAX_LEVEL, offset, R0, ubVIdxOffset );
					CHECK( frac <= 1 && frac >= 0 );
					fracSum += frac;
					fracGlucoseAvail += frac*v_gridModelRealNbrBox[GRID_MODEL_REAL_GLUCOSE_FRAC_AVAILABLE].getVal( xOffset, yOffset, 0);
				}
			}
		}
		
		
		// fraction can be off if overlapping uninhabitable regions, lets see by how much
		CHECK( fracSum <= ( 1 + GEN_SMALL ) );
		if ( fracSum < ( 1 - GEN_SMALL ) && WRITE_WARNING ) {
			cout << "WARNING:MRA_USAS_0001 - Not all of agent volume is accounted for, frac = "<<fracSum<<endl;
		}


		/* ---Growth--- 
		Assuming exponential growth
		limit growth by fraction of glucose available 
		*/
		if (fracGlucoseAvail > 0.0) {
			
			// FD to calc step size
			// allowing multiple situations:
			// 1 No Junction, no bud:
			if ( junctionInfo.getNumJunctions() == 0 ){
				deltaV = STATE_AND_GRID_TIME_STEP * GROWTH_RATE * V0;
			}
			else if ( junctionInfo.getNumJunctions() == 1 ){
				// has bud
				if ( state.getModelInt(YEAST_CELL_MODEL_INT_MOTHER) ==  0 ){
					/*2  -Daughter growth-
					grom at combined rate of this agent + mother
					where we assume mother vol set in model real at bud formation
					*/
					deltaV = STATE_AND_GRID_TIME_STEP * GROWTH_RATE * ( V0 + state.getModelReal(YEAST_CELL_MODEL_REAL_MOTHER_VOL) );
				}
				else if ( state.getModelInt(YEAST_CELL_MODEL_INT_MOTHER) ==  1 ){
					// 3, mothers do not grow
					deltaV = 0.0;
				}
				else {

					ERROR( "Mother int not set correctly." );
				}

			}
			else {
				ERROR( "Multiple buds, no logic to deal with this." );
			}


			// account for reduced glucose avalible
			deltaV *= fracGlucoseAvail;
			// To avoid error, only grow if you are not at maximum size.
			if ((V0 + deltaV)<(CELL_VOL_MAX-GEN_EPS)){
				state.setRadius( CBRT( ( V0 + deltaV ) / GEN_PI43 ) );
			}
			else if (WRITE_WARNING==true){
				cout <<"WARNING:MRA_USAS_0002 - cell near maximum size, growth prevented" <<endl;
			}
		}


		/* ---Cell Cycle Progression--- 
		Assuming critical size, const theta model
		Charvin2009
		Limit progression by fraciton of glucose available.
		*/
		if (fracGlucoseAvail > 0 ) {
			// check critical V:
			if ( V0 >= VOL_CRITICAL ) {
				// move internal clock:
				state.incModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, ( CC_CLOCK_RATE * STATE_AND_GRID_TIME_STEP * fracGlucoseAvail ) );
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
		if ( vIdx[0] == (UB_NUM[0] - 2) ){ // adjusted to remove one layer of IFs, which are a buffer for tech reasons
			dist = 0.5 * IF_GRID_SPACING - spAgent.vOffset[0];
		}

		if ( r > dist ){
			disappear = true;
		
		/* ---Divide Decision--- 
		Assuming critical, prefect timer,
		Division is equivlant to bud formation.
		If bud is already formed we do not form another,
		here we account for the by assuming a single bud;
		therefore, if a bud exists do not make another.
		*/
		}
		else if ( spAgent.state.getModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK ) >= CC_CLOCK_G1 && spAgent.junctionInfo.getNumJunctions()==0 && spAgent.state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 0 ){
			//cout <<spAgent.state.getModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK ) <<endl;
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
	const ExtraMechIntrctData extraData = mechIntrctData.extraMechIntrctData;


	// check bud conditions:
	if ( state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 1 ) {
		if ( junctionInfo.getNumJunctions()==0 ) {
			// this should not happen unless bud moved too far away, which is an error.
			state.setModelInt( YEAST_CELL_MODEL_INT_MOTHER, 0 );
			state.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, 0.0 );
			if ( WRITE_WARNING == true ) {
				cout <<"WARNING:MRA_ASAS_0005 - Mother found with no bud, set as not mother.  Reduce BL step size" <<endl;
			}
		}
		else if ( junctionInfo.getNumJunctions()>1 ) {
			ERROR("Mother found with multiple buds, no logic to handel this.");
		}
	}




	//cout << vOffset[2] <<endl;	
	/* ---Cell Boundry Collision---
	The cell boundry collision will consist of shoving.
	We must specificly look for boundry with chip desing (dummy spaces).
	We will not look at boundry of x high or low domain.
	Allowing walls to be a little closer to make up for large grid resolution.
	*/



	REAL dist;
	VIdx vIdxTmp;
	S32 xUBOffset;
	S32 yUBOffset;
	BOOL aa_isHabitable[3][3] = {{false,false,false},{false,false,false},{false,false,false}};
	// uninhabitable grids represent bounds
	getIsHabitableArray( vIdx, aa_isHabitable );
/*
	// loop though surounding grids
	for ( S32 x_i = 0; x_i < 3; x_i++ ) {
		for ( S32 y_i = 0; y_i < 3; y_i++ ) {
			xUBOffset = x_i - 1;
			yUBOffset = y_i - 1;
			// is this a boundry:
			if (aa_isHabitable[x_i][y_i]==false){
				computeGridDelta( {xUBOffset,yUBOffset}, vOffset, R0, const REAL IF )
			}

*/
	


	// impact with y bound low
	if ( vIdx[1]==0 ){
		dist = vOffset[1] + 0.5 * IF_GRID_SPACING - ADD_WALL;
		if ( dist < R0 ) {
			// hitting low bound, push up
			force[1] += CELL_WALL_STIFF * (R0 - dist);
		}

	}
	// impact with y bound high
	if ( vIdx[1] == (UB_NUM[1] - 1) ){
		dist = 0.5 * IF_GRID_SPACING - vOffset[1]  - ADD_WALL;
		if (dist < R0) {
			//hitting low bound, push down 
			force[1] += CELL_WALL_STIFF * (dist - R0);
		}
	}
	// impact with chip wall
	for( S32 dim = 0 ; dim < 2 ; dim++ ) {
		if ( vIdx[dim]>0 ){
			// check below
			dist = vOffset[dim] + 0.5 * IF_GRID_SPACING  - ADD_WALL;
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
			dist = 0.5 * IF_GRID_SPACING - vOffset[dim]  - ADD_WALL;
			vIdxTmp = vIdx;
			vIdxTmp[dim] = vIdxTmp[dim]+1;
			// if within range and its not habitable...
			if ( ( CHIP_DESIGN_MATRIX[vIdxTmp[0]][vIdxTmp[1]] == 0 ) && ( dist < R0 ) ){
				// hitting upper bound, push down
				force[dim] += CELL_WALL_STIFF * (dist - R0);
			}
		}
	}





	/* ---Extra Forces--- */



	force[0] += extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_X );
	force[1] += extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_ADD_FORCE_Y );

	/* ---Rotation--- 
	Currently no realisitc force mode, 
	we are just tracking the direction for future budding.
	*/
	if ( (FABS(extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X )) + FABS(extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X ))) > GEN_EPS ) {
		state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_X ) );
		state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y, extraData.getModelReal( EXTRA_MECH_YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) );
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
			if ( WRITE_WARNING == true ) {
				cout <<"WARNING:MRA_ASAS_0001 - Cell moving too much" <<endl;
			}
		}
		else if ( disp[dim] < (-1.0 * MAX_DISP) ) {
			disp[dim] = -1.0 * MAX_DISP;
			if ( WRITE_WARNING == true ) {
				cout <<"WARNING:MRA_ASAS_0001 - Cell moving too much" <<endl;
			}
		}
	}
	// no z-movment allowed unless its to move back to zero, warn if this happens:
	if ( vIdx[2] == 0 && FABS( vOffset[2] ) < GEN_SMALL ) {
		// in correct place, do not allow disp
		if ( FABS( disp[2] ) > GEN_SMALL ) {
			disp[2] = 0.0;
			if ( WRITE_WARNING == true ) {
				cout <<"WARNING:MRA_ASAS_0003 - Cell was set to move in z, corrected" <<endl;
			}
		}
	}
	else if ( vIdx[2] == 0 && FABS( vOffset[2] ) > GEN_SMALL ) {
		disp[2] = -1.0 * vOffset[2];
		if ( WRITE_WARNING == true ) {
			cout <<"WARNING:MRA_ASAS_0004 - Cell was found off center in z, corrected" <<endl;
		}

	}
	else {
		ERROR( "cell moved out of z=0 UB grid.");
	}




	/* ---Mitosis--- 
	In this case this is bud realeas at point in clock,
	here simply set back to not mother,
	junctionEnd will be removed else where.
	*/
	if ( state.getModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK ) >= CC_CLOCK_CRITICAL ) {
		CHECK( state.getModelInt( YEAST_CELL_MODEL_INT_MOTHER ) == 1 ); // no logic should prevent this
		state.setModelInt( YEAST_CELL_MODEL_INT_MOTHER, 0 );
		state.setModelReal( YEAST_CELL_MODEL_REAL_CC_CLOCK, 0.0 );



		// add a small random change for next division of mother, near, but not at scar
		REAL tmpDirX = state.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) + BUD_DIR_SCALE * 2.0 * (Util::getModelRand( MODEL_RNG_UNIFORM ) - 0.5 );
		REAL tmpDirY = state.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) + BUD_DIR_SCALE * 2.0 * (Util::getModelRand( MODEL_RNG_UNIFORM ) - 0.5 );
		REAL scale = SQRT( tmpDirX * tmpDirX + tmpDirY * tmpDirY );
		state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X, ( tmpDirX / scale ) );
		state.setModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y, ( tmpDirX / scale ) );


	}

	

	/* MODEL END */

	return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const AgentJunctionInfo& junctionInfo, const VReal& offset, const AgentMechIntrctData& mechIntrctData, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, SpAgentState& motherState, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	cout <<"in divide" <<endl;
	/* MODEL START */

	/* Only allow division one at a time,
	if we still have a junciton with other agent,
	ie another bud, then there is an error.
	In reality this maight actually occur,
	but current logic only handels one junciton.
	*/
	//cout <<junctionInfo.getNumJunctions() <<endl;

	if ( junctionInfo.getNumJunctions() > 0 ){
		ERROR( "New bud formed before existing daughter separated." );
	}
	//cout << "divide" <<endl;



	REAL R0 = motherState.getRadius();
	REAL V0 = GEN_PI43 * R0 * R0 * R0;
	REAL Vd = 0;

	/* ---Divide--- 
	assuming mother critical volume

	*/
	daughterState.setType( motherState.getType() );
	
	
	/* -Volume- 
	A portion of new mass goes to bud, pre set as constant
	However, if mother is too large and risks exceeding max size in next G1
	additioanl volume will be moved to bud to maintain mass conservation
	*/
	 
	if (V0 > CELL_VOL_BUD_MAX) {
		// set radius assuming all V > V_bud_max goes to d
		Vd = V0 - CELL_VOL_BUD_MAX;
		if (WRITE_WARNING==true){
		// will always happen if cell is in chamber long enough
		// not super important as mass is still conserved.
		//	cout <<"WARNING:MRA_DSA_0001 - mother over maximum size at budding, extra volume to daughter" <<endl;
		}
	}
	else {
		Vd = VOL_MOTHER_BUD;
	}
	// do not allow mother to go bellow v critical
	CHECK( (V0-Vd) > VOL_CRITICAL);	

	motherState.setRadius( CBRT( (V0-Vd) / GEN_PI43 ) );
	daughterState.setRadius( CBRT( Vd / GEN_PI43 ) );
	// displacement determined by bud dir and radius
	motherDisp[0] = -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) * ( motherState.getRadius() - BUD_OVERLAP );
	motherDisp[1] = -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) * ( motherState.getRadius() - BUD_OVERLAP );
	motherDisp[2] = 0.0;
	daughterDisp[0] = motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) * daughterState.getRadius();
	daughterDisp[1] = motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) * daughterState.getRadius();
	daughterDisp[2] = 0.0;

	// set up state for daughter	
	CHECK( NUM_YEAST_CELL_MODEL_REALS==4 );
	for ( S32 i = 0; i < NUM_YEAST_CELL_MODEL_REALS; i++ ) {
		daughterState.setModelReal( i, 0 );
	}

	daughterState.setModelReal(YEAST_CELL_MODEL_REAL_BUD_DIR_X, -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_X ) );
	daughterState.setModelReal(YEAST_CELL_MODEL_REAL_BUD_DIR_Y, -1.0 * motherState.getModelReal( YEAST_CELL_MODEL_REAL_BUD_DIR_Y ) );
	daughterState.setModelReal(YEAST_CELL_MODEL_REAL_MOTHER_VOL, (V0-Vd)  );

	// set agent to fully occupy the center box

	CHECK( NUM_YEAST_CELL_MODEL_INTS == 1 );
	daughterState.setModelInt( YEAST_CELL_MODEL_INT_MOTHER, 0 );

	// mark mother as mother for future use:
	motherState.setModelInt( YEAST_CELL_MODEL_INT_MOTHER, 1 );

	/* ---Bud Link--- */
	motherDaughterLinked = true;
	motherEnd.setType( JUNCTION_END_BUD );
	daughterEnd.setType( JUNCTION_END_BUD );

	cout <<"out of divide" <<endl;

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
/* Extra Utility Functions */


static void getIsHabitable( const VIdx vIdx, BOOL& isHabitable ){
	/* Custome for this model
	finds if the box indexed in vIdx is habitable
	accounts for boundries and chip design.
	*/
	isHabitable = false;
	if ( vIdx[0] >= 0 && vIdx[0] <= (UB_NUM[0] - 1) && vIdx[1] >= 0 && vIdx[1] <= (UB_NUM[1] - 1) && vIdx[2]==0 ){
		if ( CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]] == 1 ) {
			isHabitable = true;
		}
	}
	return;
}

static void getIsHabitableArray( const VIdx vIdx, BOOL aa_isHabitable[3][3] ){

	BOOL tmp = 0;
	VIdx vIdxTmp;

	for ( S32 x_i = 0; x_i < 3; x_i++ ){
		for ( S32 y_i = 0; y_i < 3; y_i++ ){
			vIdxTmp = vIdx;
			vIdxTmp[0] += x_i-1;
			vIdxTmp[1] += y_i-1;
			getIsHabitable( vIdxTmp, tmp );
			aa_isHabitable[x_i][y_i] = tmp;
		}
	}
}

static REAL computeGridDelta( const REAL dir[2], const VReal vOffset, const REAL r, const REAL IF ) {
	/* Compute specified grid delta, where delta is the overlap distance with the grid.
	Specifically delta is the difference between the radius and the shortest distance to the 
	specified grid boundary, edge for direct neighbors and corner for corner neighbors.
	If the agent is a sphere then delta is the maximum distance within the grid on a line
	that goes from the agent center to the agent boundary, or maximum overlap distance.
	Can be used to calculate overlap for mechanical force interactions, or overlap
	for use in estimating the occupation fraction.
	Use dir to specify neighbor grid as upper (1) or lower (-1) for each dimension (x,y) or zero to ignore that dimension.
	vOffset is agent offset at center box
	r is agent radius 
	IF is the interface Grid spacing 
	Will return zero for the center UB {0,0}.
	*/

	/* this could be done in shorter code by calculating the DHat and setting dist to sqrt(sum(DHat**2))
	but we want to avoid squaring and rooting when possible.
	*/
	REAL delta;

	if (dir[0] == 0 && dir[1] == 0 ) {
		delta = 0.0;
	}
	else {
		REAL DHat[2] = {0,0};
		REAL b = 0;
		REAL dist = 0;

		for( S32 dim = 0 ; dim < 2 ; dim++ ) {
			if ( dir[dim]== 1 || dir[dim] == -1 ) {
				// boundary location in this dim
				b = ( ( REAL ) dir[dim] ) * 0.5 * IF;
				// subtract locations in this dim
				DHat[dim] = vOffset[dim] - b;
			}
			else if ( dir[dim] != 0 ) {
				cout << dir[dim] <<endl;
				ERROR ( "Incorrect neighbor specification in dir in computeGridDelta" );
			}
		}
		if ( dir[0]==0 ) {
			//dim 1 must be not zero
			dist = FABS(DHat[1]);
		}
		else if ( dir[1]==0 ) {
			//dim 0 must be not zero
			dist = FABS(DHat[0]);
		}
		else {
			// corner
			dist = SQRT ( DHat[0]*DHat[0] + DHat[1]*DHat[1] );
		}

		if ( r > dist ) {
			delta = r - dist;
		}
		else {
			delta = 0.0;
		}


	}
	return( delta );

}

	

