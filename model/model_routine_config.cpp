/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
	/* MODEL START */

	ifGridSpacing = IF_GRID_SPACING;

	/* MODEL END */

	return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
	/* MODEL START */

	callInfo.numUpdateIfGridVarPreStateAndGridStepRounds = 1;
	callInfo.numUpdateIfGridVarPostStateAndGridStepRounds = 1;
//	callInfo.numUpdateIfGridVarPreStateAndGridSteps = 1;
//	callInfo.numUpdateIfGridVarPostStateAndGridSteps = 1;

	/* MODEL END */

	return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
	/* MODEL START */

	a_domainBdryType[0] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;
	a_domainBdryType[1] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;
	a_domainBdryType[2] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBdryType( pde_buffer_bdry_type_e& pdeBufferBdryType ) {
	/* MODEL START */

	pdeBufferBdryType = PDE_BUFFER_BDRY_TYPE_HARD_WALL;

	/* MODEL END */

	return;
}

void ModelRoutine::updateTimeStepInfo( TimeStepInfo& timeStepInfo ) {
	/* MODEL START */

	timeStepInfo.durationBaselineTimeStep = BASELINE_TIME_STEP_DURATION;
	timeStepInfo.numStateAndGridTimeStepsPerBaseline = NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

	/* MODEL END */

	return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& extraMechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
	/* MODEL START */

	extraMechIntrctSyncMethod = SYNC_METHOD_DELTA;
	updateIfGridVarSyncMethod = SYNC_METHOD_DELTA;

	/* MODEL END */

	return;
}

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {
	/* MODEL START */

	CHECK( NUM_AGENT_TYPES == 1 );

	SpAgentInfo info;

	v_spAgentInfo.resize( NUM_AGENT_TYPES );

	info.dMax = CELL_INTRCT_DIST_MAX;
	info.hasBool = false;
	info.numBoolVars = 0;
	info.numStateModelReals = NUM_YEAST_CELL_MODEL_REALS;
	info.numStateModelInts = 0;
	info.numExtraMechIntrctModelReals = 0;
	info.numExtraMechIntrctModelInts = 0;
	info.v_odeNetInfo.clear();

	v_spAgentInfo[AGENT_YEAST_CELL] = info;

	
	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateJunctionEndInfo( Vector<JunctionEndInfo>& v_junctionEndInfo ) {
	/* MODEL START */

	v_junctionEndInfo.clear();

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEInfo( Vector<PDEInfo>& v_pdeInfo ) {
	/* MODEL START */

	CHECK( NUM_DIFFUSIBLE_ELEMS == 1 );

	PDEInfo pdeInfo;
	GridPhiInfo gridPhiInfo;

	v_pdeInfo.resize( NUM_DIFFUSIBLE_ELEMS );

	/* stand in element 1 */

	pdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_TIME_DEPENDENT_LINEAR;
	pdeInfo.numLevels = NUM_AMR_LEVELS;
	pdeInfo.numTimeSteps = NUM_PDE_TIME_STEPS_PER_STATE_AND_GRID_STEP;
	pdeInfo.callAdjustRHSTimeDependentLinear = false;

	gridPhiInfo.elemIdx = DIFFUSIBLE_ELEM_GLUCOSE;
	gridPhiInfo.name = "glucose";
	gridPhiInfo.aa_bcType[0][0] = BC_TYPE_DIRICHLET_CONST;
	gridPhiInfo.aa_bcVal[0][0] = ELEM_BULK_CONCENTRATION[DIFFUSIBLE_ELEM_GLUCOSE];
	gridPhiInfo.aa_bcType[0][1] = BC_TYPE_DIRICHLET_CONST;
	gridPhiInfo.aa_bcVal[0][1] = ELEM_BULK_CONCENTRATION[DIFFUSIBLE_ELEM_GLUCOSE];
	gridPhiInfo.aa_bcType[1][0] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[1][0] = 0.0;
	gridPhiInfo.aa_bcType[1][1] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[1][1] = 0.0;
	gridPhiInfo.aa_bcType[2][0] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[2][0] = 0.0;
	gridPhiInfo.aa_bcType[2][1] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[2][1] = 0.0;

	gridPhiInfo.minVal = Info::getMGParabolicNormThreshold() * -1.0;
	gridPhiInfo.setNegToZero = true;
	gridPhiInfo.printSummary = true;
	gridPhiInfo.summaryType = SUMMARY_TYPE_SUM;
	gridPhiInfo.fileOutput = true;

	pdeInfo.v_gridPhiInfo.assign( 1, gridPhiInfo );

	v_pdeInfo[DIFFUSIBLE_ELEM_GLUCOSE] = pdeInfo;


	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {
	/* MODEL START */

	CHECK( NUM_GRID_MODEL_REALS == 2 );

	IfGridModelVarInfo info;

	v_ifGridModelRealInfo.resize( NUM_GRID_MODEL_REALS );

	info.name = "glucose_delta";
	info.printSummary = false;
	info.summaryType = SUMMARY_TYPE_SUM;/* dummy */

	v_ifGridModelRealInfo[GRID_MODEL_REAL_GLUCOSE_DELTA] = info;

	info.name = "agent_volume";
	info.printSummary = true;
	info.summaryType = SUMMARY_TYPE_SUM;

	v_ifGridModelRealInfo[GRID_MODEL_REAL_AGENT_VOL] = info;


	CHECK( NUM_GRID_MODEL_INTS == 1 );


	v_ifGridModelIntInfo.resize( NUM_GRID_MODEL_INTS );

	info.name = "glucose_avalible";
	info.printSummary = false;
	info.summaryType = SUMMARY_TYPE_SUM;/* dummy */

	v_ifGridModelIntInfo[GRID_MODEL_INT_GLUCOSE_AVAILABLE] = info;

	/* MODEL END */

	return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
	/* MODEL START */

	CHECK( NUM_MODEL_RNGS == 1 );

	v_rngInfo.resize( 1 );

	RNGInfo rngInfo;

	rngInfo.type = RNG_TYPE_UNIFORM;
	rngInfo.param0 = 0.0;
	rngInfo.param1 = 1.0;
	rngInfo.param2 = 0.0;/* dummy */

	v_rngInfo[MODEL_RNG_UNIFORM] = rngInfo;

	/* MODEL END */

	return;
}

void ModelRoutine::updateGlobalData( Vector<U8>& v_globalData ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::init( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::term( void ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::setPDEBuffer( const VIdx& startVIdx, const VIdx& regionSize, BOOL& isPDEBuffer ) {
	/* MODEL START */

	isPDEBuffer = false;
	/* MODEL END */

	return;
}

void ModelRoutine::setHabitable( const VIdx& vIdx, BOOL& isHabitable ) {
	/* MODEL START */


	CHECK( Info::getDomainSize( 0 ) == UB_NUM[0] );
	CHECK( Info::getDomainSize( 1 ) == UB_NUM[1] );

	/* us the predefined chip matrix to 
	test if this in habitable
	*/
	if ( vIdx[2] > 0 ) {
		 isHabitable = false;
	}
	else if ( CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]] == 0 ){
		isHabitable = false;
	}
	else if ( CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]] == 1 ) {
		isHabitable = true;
	}
	else {
		ERROR("unknown int in chip design matrix")
	}


	/* MODEL END */

	return;
}


