#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include <TH2D.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct r2{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=(aod::track::pt>0.2f)&&(aod::track::pt<2.0f);		
	Filter globalfilter=requireGlobalTrackInFilter();

	
	HistogramRegistry histos{"hist",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{

		histos.add("BB","#frac{dE}{dx} vs p_T",kTH2D,{{300,0.2,2,"p_T"},{500,0,300,"#frac{dE}{dx}"}});
		histos.add("Beta","#beta vs p_T",kTH2D,{{300,0.2,2,"p_T"},{300,0.3,1.5,"#beta"}});
	}
void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::TracksExtra,aod::pidTOFbeta>> const& tracks)
	{
		for(auto track:tracks)
		{
			histos.fill(HIST("BB"),track.pt(),track.tpcSignal());
			histos.fill(HIST("Beta"),track.pt(),track.beta());
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<r2>(cfgc)
    };
}
