#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct p2pre
{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
//	Filter ptfilter=(aod::track::pt>0.2f)&&(aod::track::pt<2.0f);		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	
	HistogramRegistry histos{"pT",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const&)
	{
		histos.add("pt","pt",kTH2D,{{100,-0.8,0.8,"eta"},{100,0,2*PI,"phi"}});
		histos.add("rho1","rho1",kTH2D,{{100,-0.8,0.8,"eta"},{100,0,2*PI,"phi"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
	{
		for(auto track:tracks)
		{
			histos.fill(HIST("pt"),track.eta(),track.phi(),track.pt());
			histos.fill(HIST("rho1"),track.eta(),track.phi());
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<p2pre>(cfgc)
    };
}
