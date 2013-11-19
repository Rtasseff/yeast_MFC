/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, Vector<REAL>& v_gridPhi/* [elemIdx] */, Vector<REAL>& v_gridModelReal/* [elemIdx] */, Vector<S32>& v_gridModelInt ) {
	/* MODEL START */

	CHECK( NUM_DIFFUSIBLE_ELEMS == 1 );
	CHECK( NUM_GRID_MODEL_REALS == 2 );

	v_gridModelReal[GRID_MODEL_REAL_GLUCOSE_DELTA] = 0.0;
	v_gridModelReal[GRID_MODEL_REAL_AGENT_VOL] = 0.0;
	
	if (CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]]==1 && vIdx[2]==0){
		v_gridPhi[DIFFUSIBLE_ELEM_GLUCOSE] = ELEM_BULK_CONCENTRATION[DIFFUSIBLE_ELEM_GLUCOSE];
	}
	else {
		v_gridPhi[DIFFUSIBLE_ELEM_GLUCOSE] = 0.0;
	}


	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 round, const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* INOUT, [elemIdx] */, Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* INOUT, [elemIdx] */ ) {
	/* MODEL START */

	CHECK( NUM_DIFFUSIBLE_ELEMS == 1 );
	// only have one round set
	CHECK( round == 0 );

	// onle z=0 is real, others are dummies
	if( vIdx[2]==0 ){
		if( pre == true ) {
			/* ---Agent Uptake--- 
			Only glucose, amount listed in model real,
			add up the total amount taken in this step
			*/

			REAL glucoseTotalDelta = 0;
			
			const UBAgentData& ubAgentData = *( ubAgentDataPtrNbrBox.getVal( 0, 0, 0 ) );
			for( S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
				const SpAgentState& state = ubAgentData.v_spAgent[i].state;
				if ( state.getType() == AGENT_YEAST_CELL  ){
					// note cells are taking out of grid, consuming 
					glucoseTotalDelta -= state.getModelReal( YEAST_CELL_MODEL_REAL_ELEM_GLUCOSE_UPTAKE );
				}
			}
			/* well we are going to set it, assuming nothing else will change it */
			v_gridModelRealNbrBox[GRID_MODEL_REAL_GLUCOSE_DELTA].setVal( 0, 0, 0, glucoseTotalDelta );
			
		}
		else { /* post */
			/* calculate cell volume */
			REAL agentVol = 0;
			const UBAgentData& ubAgentData = *( ubAgentDataPtrNbrBox.getVal( 0, 0, 0 ) );
			for( S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
				const SpAgentState& state = ubAgentData.v_spAgent[i].state;
				REAL r = state.getRadius();
				agentVol += r * r * r * GEN_PI43;
			}
			v_gridModelRealNbrBox[GRID_MODEL_REAL_AGENT_VOL].setVal( 0, 0, 0, agentVol );

		}
	}
	else { // looking at dummy space
		if (pre == true){
			if ( v_gridPhiNbrBox[DIFFUSIBLE_ELEM_GLUCOSE].getVal(0,0,0) > GEN_SMALL ){
				cout <<"WARNING:MRG_UIGV_0002 - glucose in dummy region!." <<endl;
			}
		}
	}




	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridKappa( const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	if ( CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]]==1 && vIdx[2] == 0){
	// legitimate region 

		REAL volUB = IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING;
		REAL tmpKappa = ( volUB - v_gridModelReal[GRID_MODEL_REAL_AGENT_VOL] ) / volUB ;
		
		if ( tmpKappa > GEN_SMALL ) {
			gridKappa = tmpKappa;
		}
		else {
			gridKappa = GEN_SMALL;
			cout <<"WARNING:MRG_UIGK_0001 - Kappa was calculated as negative." <<endl;
		}

	}
	else {
		gridKappa = 1;
	}
	

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAlpha( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */

	/* check my understanding */
	VIdx vIdx1Tmp = vIdx0;
	vIdx1Tmp[dim] = vIdx1Tmp[dim] + 1;
	CHECK( vIdx1Tmp == vIdx1 );

	/* set to defulat */
	gridBeta = ELEM_BETA[elemIdx];

	/* using the desing matrix to check if we are at a wall (no diffusion in z)*/
	if ( CHIP_DESIGN_MATRIX[vIdx0[0]][vIdx0[1]] == 0 || CHIP_DESIGN_MATRIX[vIdx1[0]][vIdx1[1]] == 0 || vIdx0[2] > 0 || vIdx1[2] > 0 ) {
		// no diffusion at chip wall
		gridBeta = 0;
	}



	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 0;

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridBeta ) {
	/* MODEL START */
	/* set defualt as zero */
	gridBeta = 0;

	// check if dummy region
	if ( vIdx[2] == 0 ) {
		/* look to see if we are at the flow channel */
		if ( dim == 0 && CHIP_DESIGN_MATRIX[vIdx[0]][vIdx[1]] == 1 ) {
			gridBeta = ELEM_BETA[elemIdx];
		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& vIdx1, const UBAgentData& ubAgentData0, const UBAgentData& ubAgentData1, const Vector<REAL>& v_gridPhi0/* [elemIdx] */, const Vector<REAL>& v_gridPhi1/* [elemIdx] */, const Vector<REAL>& v_gridModelReal0/* [elemIdx] */, const Vector<REAL>& v_gridModelReal1/* [elemIdx] */, const Vector<S32>& v_gridModelInt0/* [elemIdx] */, const Vector<S32>& v_gridModelInt1/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<REAL>& v_gridPhi/* [elemIdx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	/* Given total change over this step, calulate source term
	mass/(vol time)
	DO NOT ALLOW NEGATIVVE
	*/
	CHECK( NUM_DIFFUSIBLE_ELEMS == 1 );

	REAL total =  v_gridPhi[DIFFUSIBLE_ELEM_GLUCOSE] * IF_GRID_SPACING * IF_GRID_SPACING *IF_GRID_SPACING; // ng in UB
	REAL delta = v_gridModelReal[GRID_MODEL_REAL_GLUCOSE_DELTA];

	/* well we are going to set it, assuming nothing else will change it */
	if ( ( total + delta ) < 0 ) {
		// taking all of the glucose
		delta = -1.0 * total;
		cout <<"WARNING:MRG_UIGRHSL_0002 - Taking too much glucose, correcting." <<endl;
	}

	gridRHS = delta / ( IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING * STATE_AND_GRID_TIME_STEP );

	/* MODEL END */

	return;
}

void ModelRoutine::adjustIfGridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const REAL gridPhi, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, REAL& gridRHS/* INOUT, uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const UBAgentData& ubAgentData, const Vector<double>& v_gridPhi/* [idx] */, const Vector<REAL>& v_gridModelReal/* [elemIdx] */, const Vector<S32>& v_gridModelInt/* [elemIdx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrBox<const UBAgentData*>& ubAgentDataPtrNbrBox, const Vector<NbrBox<REAL> >& v_gridPhiNbrBox/* [elemIdx] */, const Vector<NbrBox<REAL> >& v_gridModelRealNbrBox/* [elemIdx] */, const Vector<NbrBox<S32> >& v_gridModelIntNbrBox/* [elemIdx] */, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
	/* MODEL START */

	for( S32 pdeIdx = 0 ; pdeIdx < NUM_DIFFUSIBLE_ELEMS ; pdeIdx++ ) {
		v_finestLevel[pdeIdx] = NUM_AMR_LEVELS - 1;

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

	CHECK( NUM_DIFFUSIBLE_ELEMS==1 );
	v_gridPhi[DIFFUSIBLE_ELEM_GLUCOSE] = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridKappa/* to consider cell volume exclusion in computing diffusion flux, set to 1.0 to ignore the volume occupied by cells */ ) {
	/* MODEL START */

	gridKappa = 1.0;
	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAlpha( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAlpha/* decay (-) */ ) {
	/* MODEL START */

	gridAlpha = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
	/* MODEL START */

	gridBeta = 0.0;

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAdvectionVelocityInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplemented." );
	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferAdvectionVelocityDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAdvectionVelocity ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	gridRHS = 0.0;
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

void ModelRoutine::updatePDEBufferDirichletBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferFaceSize, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], REAL& bcVal ) {
	/* note that a_gridPhi[2] is valid only when the domain size in the direction of dim is not smaller than 3 * (refinementment ratio between the coarsest AMR level and the finest AMR level) */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferNeumannBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferFaceSize, const S32 dim, const BOOL lowSide, const REAL a_gridPhi[3], REAL& bcVal ) {
	/* note that a_gridPhi[2] is valid only when the domain size in the direction of dim is not smaller than 3 * (refinementment ratio between the coarsest AMR level and the finest AMR level) */

	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

