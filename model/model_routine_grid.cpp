/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include <math.h>

#include <iostream>

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, Vector<REAL>& v_gridPhi, Vector<REAL>& v_gridModelReal, Vector<S32>& v_gridModelInt ) {
	/* MODEL START */

	CHECK( NUM_DIFFUSIBLE_ELEMS == 1 );
	CHECK( NUM_GRID_MODEL_REALS == 1 );

	v_gridPhi[DIFFUSIBLE_ELEM0] = ( REAL )( vIdx[0] + vIdx[1] + vIdx[2] ) + 1e-10;
	v_gridModelReal[GRID_MODEL_REAL_FV_RATIO] = 1.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 round, const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* INOUT, [elemIdx] */, Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* INOUT, [elemIdx] */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridKappa( const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	CHECK( v_gridModelReal[GRID_MODEL_REAL_FV_RATIO] > 0.0 );
	gridKappa = v_gridModelReal[GRID_MODEL_REAL_FV_RATIO];

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAlpha( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0, const Vector<REAL>& v_gridPhi1, const Vector<REAL>& v_gridModelReal0, const Vector<REAL>& v_gridModelReal1, const Vector<S32>& v_gridModelInt0, const Vector<S32>& v_gridModelInt1, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 1e-1;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 1e-1;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 1e-1;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0, const Vector<REAL>& v_gridPhi1, const Vector<REAL>& v_gridModelReal0, const Vector<REAL>& v_gridModelReal1, const Vector<S32>& v_gridModelInt0, const Vector<S32>& v_gridModelInt1, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	gridAdvectionVelocity = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	gridAdvectionVelocity = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	gridAdvectionVelocity = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	gridRHS = ( REAL )( vIdx[0] + vIdx[1] + vIdx[2] ) / ( 32.0 * 32.0 * 32.0 );

	/* MODEL END */

	return;
}

void ModelRoutine::adjustIfGridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const REAL gridPhi, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<double>& v_gridPhi/* [idx] */, const Vector<REAL>& v_gridModelReal, const Vector<S32>& v_gridModelInt, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
	/* MODEL START */

	for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
		v_finestLevel[pdeIdx] = 0;
		if( vIdx[2] == 16 ) {
			if( ( ( vIdx[0] * 3785 + vIdx[1] * 2931 + 7 ) % 40 ) == 0 ) {
				cout << "tagging." << endl;
				v_finestLevel[pdeIdx] = NUM_AMR_LEVELS - 1;
			}
		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridDirichletBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], const Vector<Vector<REAL> >& vv_gridModelReal/* [elemIdx] */, const Vector<Vector<S32> >& vv_gridModelInt/* [elemIdx] */, REAL& bcVal ) {
	/* note that a_gridPhi[2], vv_gridModelReal[elemIdx][2], and vv_gridModelInt[elemIdx][2] are valid only when the domain size in the direction of dim is not smaller than 3 */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridNeumannBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], const Vector<Vector<REAL> >& vv_gridModelReal/* [elemIdx] */, const Vector<Vector<S32> >& vv_gridModelInt/* [elemIdx] */, REAL& bcVal ) {
	/* note that a_gridPhi[2], vv_gridModelReal[elemIdx][2], and vv_gridModelInt[elemIdx][2] are valid only when the domain size in the direction of dim is not smaller than 3 */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::initPDEBufferPhi( const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, Vector<REAL>& v_gridPhi ) {
	/* MODEL START */

	v_gridPhi[DIFFUSIBLE_ELEM0] = 0.0;
	for( S32 i = 0 ; i < pdeBufferBoxSize[0] ; i++ ) {
		for( S32 j = 0 ; j < pdeBufferBoxSize[1] ; j++ ) {
			for( S32 k = 0 ; k < pdeBufferBoxSize[2] ; k++ ) {
				VIdx vIdx = startVIdx;
				vIdx[0] += i;
				vIdx[1] += j;
				vIdx[2] += k;
				v_gridPhi[DIFFUSIBLE_ELEM0] += ( REAL )( vIdx[0] + vIdx[1] + vIdx[2] ) + 1e-10;
			}
		}
	}
	v_gridPhi[DIFFUSIBLE_ELEM0] /= ( REAL )pdeBufferBoxSize[0] * ( REAL )pdeBufferBoxSize[1] * ( REAL )pdeBufferBoxSize[2];

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	gridKappa = 0.0;
	for( S32 i = 0 ; i < pdeBufferBoxSize[0] ; i++ ) {
		for( S32 j = 0 ; j < pdeBufferBoxSize[1] ; j++ ) {
			for( S32 k = 0 ; k < pdeBufferBoxSize[2] ; k++ ) {
				VIdx vIdx = startVIdx;
				vIdx[0] += i;
				vIdx[1] += j;
				vIdx[2] += k;
				gridKappa += 1.0 - ( REAL )( vIdx[0] + vIdx[1] + vIdx[2] ) / ( REAL )( vIdx[0] * 3 + vIdx[1] * 2 + vIdx[2] * 5 + 10 );
			}
		}
	}
	gridKappa /= ( REAL )pdeBufferBoxSize[0] * ( REAL )pdeBufferBoxSize[1] * ( REAL )pdeBufferBoxSize[2];

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAlpha( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

	gridAlpha = -0.01;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 1e-1;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 1e-1;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAdvectionVelocityInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	gridAdvectionVelocity = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	gridAdvectionVelocity = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	gridRHS = 0.0;
	for( S32 i = 0 ; i < pdeBufferBoxSize[0] ; i++ ) {
		for( S32 j = 0 ; j < pdeBufferBoxSize[1] ; j++ ) {
			for( S32 k = 0 ; k < pdeBufferBoxSize[2] ; k++ ) {
				VIdx vIdx = startVIdx;
				vIdx[0] += i;
				vIdx[1] += j;
				vIdx[2] += k;
				gridRHS = ( REAL )( vIdx[0] + vIdx[1] + vIdx[2] ) / ( 32.0 * 32.0 * 32.0 );
			}
		}
	}
	gridRHS /= ( REAL )pdeBufferBoxSize[0] * ( REAL )pdeBufferBoxSize[1] * ( REAL )pdeBufferBoxSize[2];

	/* MODEL END */

	return;
}

void ModelRoutine::adjustPDEBufferRHSTimeDependentLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferDirichletBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferBoxSize, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], REAL& bcVal ) {
	/* note that a_gridPhi[2], vv_gridModelReal[elemIdx][2], and vv_gridModelInt[elemIdx][2] are valid only when the domain size in the direction of dim is not smaller than 3 */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferNeumannBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferBoxSize, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], REAL& bcVal ) {
	/* note that a_gridPhi[2], vv_gridModelReal[elemIdx][2], and vv_gridModelInt[elemIdx][2] are valid only when the domain size in the direction of dim is not smaller than 3 */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

