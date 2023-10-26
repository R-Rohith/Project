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
		histos.add("bb2","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","bb_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","bb_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","bb_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","bb_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{400,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{400,-5,5,"nsigma"}});
		}
	
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
//		float tpccut,tofcut;
		for(auto track1:tracks)
		{	
		if(fabs(track1.tpcNSigmaPi())>=3)continue;
		histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return {
	adaptAnalysisTask<pid>(cfgc)
	};
}



