#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
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
		histos.add("rho1","rho1",kTH2D,{{24,-0.8,0.8,"eta"},{36,0,2*PI,"phi"}});
		histos.add("rho2","rho2",kTH2D,{{864,0,864,"eta1phi1"},{864,0,864,"eta2phi2"}});
		histos.add("mult","mult",kTH1I,{{200,0,200,"multiplicity"}});
		histos.add("eta","eta",kTH1D,{{100,-2,1,"eta"}});
		histos.add("2rho2","2rho2",kTH2D,{{864,0,864,"eta1phi1"},{864,0,864,"eta2phi2"}});
				
	}
void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& track)
//void process(aod::Collision const& coll,aod::TracksIU const& track)
	{
		int etabin1,phibin1,etabin2,phibin2;
		int mult=0;
		double rho1[24][36],rho2[864][864];
		for(int i=0;i<864;i++)
		{rho1[i/36][i%36]=0;
		for(int j=0;j<864;j++)
		rho2[i][j]=0;}
		for(auto track1:track)
		{
			mult++;
			if(fabs(track1.eta())>=0.8)continue;
			etabin1=(track1.eta()+0.8)*15;
			phibin1=(track1.phi()*18)/PI;
			histos.fill(HIST("rho1"),track1.eta(),track1.phi());
			rho1[etabin1][phibin1]++;
			histos.fill(HIST("eta"),track1.eta());
			for(auto track2:track)
			{
				if(track1.index()>=track2.index())continue;
				if(fabs(track2.eta())>=0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=(track2.phi()*18)/PI;
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				rho2[36*etabin1+phibin1][36*etabin2+phibin2]++;
//				histos.fill(HIST("rho2"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
			}
		}
		histos.fill(HIST("mult"),mult);
		for(int i=0;i<864;i++)
		{
			for(int j=0;j<864;j++)
			{
				histos.fill(HIST("2rho2"),i+1,j+1,rho2[i][j]/(rho1[i/36][i%36]*rho1[j/36][j%36]));
			}
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<r2>(cfgc)
    };
}
