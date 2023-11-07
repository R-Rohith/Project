#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct pid{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=(aod::track::pt>0.2f);//&&(aod::track::pt<1.3f);		//till 0.5 no delta peak....peak seems to begin from 0.25.
	//0.2 fot tpc, 0.35 for tof
	Filter globalfilter=requireGlobalTrackInFilter();
	HistogramRegistry histos{"hist",{},OutputObjHandlingPolicy::AnalysisObject};
	
	void init(InitContext const&)
	{
//		histos.add("rho1","rho1",kTH2D,{{24,-0.8,0.8,"eta"},{36,0,2*PI,"phi"}});
//		histos.add("rho2","rho2",kTH2D,{{864,0,864,"eta1phi1"},{864,0,864,"eta2phi2"}});
//		histos.add("mult","mult",kTH1I,{{200,0,200,"multiplicity"}});
//		histos.add("eta","eta",kTH1D,{{100,-1,1,"eta"}});
//		histos.add("phi","phi",kTH1D,{{100,-1,2*PI+1,"phi"}});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","bb_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","bb_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","bb_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","bb_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{400,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{400,-5,5,"nsigma"}});
		histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});
		histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{{24,-0.8,0.8,"eta"},{36,0,2.0*constants::math::PI,"phi"}});
		histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{{24,-0.8,0.8,"eta"},{36,0,2.0*constants::math::PI,"phi"}});
		}
	
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCFullEl,aod::pidTOFEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
//		float tpccut,tofcut;
		for(auto track1:tracks)
		{
/*		histos.fill(HIST("bb3"),track1.pt(),track1.tpcSignal());	
		if(fabs(track1.tpcNSigmaPi())>=3)continue;
		histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
		if(track1.tpcSignal()>80){ histos.fill(HIST("bb2_b"),track1.pt(),track1.tpcNSigmaPi());
		histos.fill(HIST("bb2_a"),track1.pt(),track1.tpcExpSignalDiffPi()/track1.tpcExpSigmaPi());}*/
			histos.fill(HIST("pt1"),track1.pt());
			histos.fill(HIST("bb1_b"),track1.pt(),track1.tpcNSigmaEl());
			histos.fill(HIST("bb2_b"),track1.pt(),track1.tofNSigmaEl());
			histos.fill(HIST("bb2"),track1.pt(),track1.beta());
			histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
			if(track1.hasTOF())
			histos.fill(HIST("pt2"),track1.pt());
			
/*//------------------------For Pion no tpc cut seems to required.....turns out its true 
			tpccut=2.5;
			tofcut=2.5;
			if(fabs(track1.tpcNSigmaPi())>=tpccut) continue;	//tpc doesn't seem to be working
			
//			{
				if(track1.hasTOF())
				{
//				if(track1.pt()<1.5)tofcut=3;
//				else tofcut=3;	
				if(fabs(track1.tofNSigmaPi())>=tofcut) continue;
				}
				else if(track1.pt()>=0.7)
					continue;
//			}
//			else			//including this creates a huge loss in low pt particle counts
//			{
//				if(track1.hasTOF())
//				{
//				tofcut=3;
//				if(fabs(track1.tofNSigmaPi())>=tofcut) continue;
//				}
				
				
//			}
//-------------------------------------------------------------------------------------needs work*/
/*//for Kaon there's a small region with e, since its in low momentum TOF is not very effective.
			tpccut=2;
			tofcut=2;
			if(track1.pt()<0.7)
			{
				if(track1.pt()<0.45)tpccut=2;
				else if(track1.pt()<0.55)tpccut=1;
				else tpccut=0.6;
				if(track1.hasTOF())
				{
					tpccut=2;
					if((fabs(track1.tofNSigmaKa())>tofcut)||(fabs(track1.tpcNSigmaKa())>tpccut))continue;
				}
				else
					if(fabs(track1.tpcNSigmaKa())>tpccut) continue;
			}
			else if(track1.hasTOF())
			{
//				if(track1.pt()<1)tofcut=3;
//				else if(track1.pt()<1.2)tofcut=2;
//				else if(track1.pt()<1.3)tofcut=1;
//				else tofcut=3;
				if((fabs(track1.tpcNSigmaKa())>tpccut)||(fabs(track1.tofNSigmaKa())>tofcut)) continue;
			}
			else
			continue;
//----------------------------------------------------------------------------------------------*/
/*//Proton----------------------------------------------------------------------------little loss
			tpccut=2.2;
			tofcut=2;
			if(track1.tpcNClsCrossedRows()<70) continue;
			if(track1.pt()<1.25)
			{
				if(track1.pt()<0.85)tpccut=2.2;
				else tpccut=1;
				if(track1.hasTOF())
				{
					tpccut=2.2;
					if((fabs(track1.tofNSigmaPr())>tofcut)||(fabs(track1.tpcNSigmaPr())>tpccut))continue;
				}
				else
					if(fabs(track1.tpcNSigmaPr())>tpccut) continue;
			}
			else if(track1.hasTOF())
			{
//				if(track1.pt()<1.6)tofcut=3;
//				else if(track1.pt()<1.9)tofcut=2;
//				else if(track1.pt()<2.4)tofcut=1;
//				else if(track1.pt()<3.5)tofcut=0;	
//				else tofcut=3;
				if((fabs(track1.tpcNSigmaPr())>tpccut)||(fabs(track1.tofNSigmaPr())>tofcut)) continue;
			}
			else
			continue;
//----------------------------------------------------------------------------------*/
//Electron---------------------------------------------------------------------
			if(fabs(track1.tpcNSigmaEl())>1.2) continue;
//----------------------------------------------------------*/			

			histos.fill(HIST("bb4"),track1.pt(),track1.beta());
			histos.fill(HIST("bb1_a"),track1.pt(),track1.tpcNSigmaEl());
			histos.fill(HIST("bb3"),track1.pt(),track1.tpcSignal());
			histos.fill(HIST("proj1"),track1.tpcNSigmaEl());
			histos.fill(HIST("bb2_a"),track1.pt(),track1.tofNSigmaEl());
			histos.fill(HIST("proj2"),track1.tofNSigmaEl());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP"),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptM"),track1.pt());
			}


		
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return {
	adaptAnalysisTask<pid>(cfgc)
	};
}
/*

o2physics_add_dpl_workflow(p2
                  SOURCES p2.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)
                  
o2physics_add_dpl_workflow(p2pre
                  SOURCES p2pre.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

o2physics_add_dpl_workflow(r2p2-4-id
                  SOURCES r2p2-4-id.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)
             
o2physics_add_dpl_workflow(mytask1
                  SOURCES mytask1.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

o2physics_add_dpl_workflow(r2
                  SOURCES r2.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

o2physics_add_dpl_workflow(slfr2
                  SOURCES slfr2.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

     
*/


