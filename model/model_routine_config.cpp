/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <math.h>

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

void ModelRoutine::updateGlobalInfo( Vector<U8>& v_globalData ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
	/* MODEL START */

	ifGridSpacing = 2.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
	/* MODEL START */

	callInfo.numUpdateIfGridVarPreStateAndGridSteps = 0;
	callInfo.numUpdateIfGridVarPostStateAndGridSteps = 0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
	/* MODEL START */

	CHECK( DIMENSION == 3 );

	for( S32 dim = 0 ; dim < 2 ; dim++ ) {
		a_domainBdryType[dim] = DOMAIN_BDRY_TYPE_PERIODIC;
	}
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

	timeStepInfo.durationBaselineTimeStep = 1.0;
	timeStepInfo.numStateAndGridTimeStepsPerBaseline = 1;

	/* MODEL END */

	return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& extraMechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
	/* MODEL START */

	extraMechIntrctSyncMethod = SYNC_METHOD_OVERWRITE;
	updateIfGridVarSyncMethod = SYNC_METHOD_OVERWRITE;

	/* MODEL END */

	return;
}

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {
	/* MODEL START */

	v_spAgentInfo.clear();

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

	pdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_TIME_DEPENDENT_LINEAR;
	pdeInfo.numLevels = NUM_AMR_LEVELS;
	pdeInfo.numTimeSteps = 2;
	pdeInfo.callAdjustRHSTimeDependentLinear = false;

	pdeInfo.advectionInfo.courantNumber = 0.5;/* dummy */

	pdeInfo.splittingInfo.v_diffusionTimeSteps.assign( 1, 1 );/* dummy */
	pdeInfo.splittingInfo.odeStiff = ODE_STIFF_NORMAL;/* dummy */
	pdeInfo.splittingInfo.odeH = 0.5;/* dummy */
	pdeInfo.splittingInfo.odeHm = 0.1;/* dummy */

	gridPhiInfo.elemIdx = DIFFUSIBLE_ELEM0;
	gridPhiInfo.name = "elem0";
	gridPhiInfo.aa_bcType[0][0] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[0][0] = 0.0;
	gridPhiInfo.aa_bcType[0][1] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[0][1] = 0.0;
	gridPhiInfo.aa_bcType[1][0] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[1][0] = 0.0;
	gridPhiInfo.aa_bcType[1][1] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[1][1] = 0.0;
	gridPhiInfo.aa_bcType[2][0] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[2][0] = 0.0;
	gridPhiInfo.aa_bcType[2][1] = BC_TYPE_NEUMANN_CONST;
	gridPhiInfo.aa_bcVal[2][1] = 0.0;

	gridPhiInfo.minVal = 0.0;
	gridPhiInfo.setNegToZero = false;
	gridPhiInfo.printSummary = true;
	gridPhiInfo.summaryType = SUMMARY_TYPE_SUM;
	gridPhiInfo.fileOutput = false;

	pdeInfo.v_gridPhiInfo.assign( 1, gridPhiInfo );
	v_pdeInfo.push_back( pdeInfo );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {
	/* MODEL START */

	IfGridModelVarInfo gridModelVarInfo;

	gridModelVarInfo.name = "free_vol_ratio";
	gridModelVarInfo.printSummary = false;
	gridModelVarInfo.summaryType = SUMMARY_TYPE_AVG;/* dummy */

	v_ifGridModelRealInfo.clear();
	v_ifGridModelRealInfo.push_back( gridModelVarInfo );

	v_ifGridModelIntInfo.clear();

	/* MODEL END */

	return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
	/* MODEL START */

	v_rngInfo.clear();

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

	isHabitable = true;

	/* MODEL END */

	return;
}

