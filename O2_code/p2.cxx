#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct p2
{
//		TFile* file=new TFile("p2pre.root","read");
//		TH2D *pt=(TH2D*)file->Get("pT/pt"),*pt_den=(TH2D*)file->Get("pT/rho1");

	HistogramRegistry histos{"P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const&)
	{
		histos.add("num","num",kTH2D,{{10000,0,10000,"eta1phi1"},{10000,0,10000,"eta2phi2"}});
		histos.add("den","den",kTH2D,{{10000,0,10000,"eta1phi1"},{10000,0,10000,"eta2phi2"}});
		histos.add("rho1","rho1",kTH2D,{{100,-0.8,0.8,"eta"},{100,0,2*PI,"phi"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
	{
//		double p1,p2;
		for(auto track:tracks)
		histos.fill(HIST("rho1"),track.eta(),track.phi());
		int etabin1,phibin1,etabin2,phibin2;
		for(auto track1:tracks)
		{
			if(fabs(track1.eta())>=0.8)continue;
			etabin1=(track1.eta()+0.8)*100/1.6;
			phibin1=track1.phi()*50/PI;
			for(auto track2:tracks)
			{
				if(fabs(track2.eta())>=0.8)continue;
				etabin2=(track2.eta()+0.8)*100/1.6;
				phibin2=track2.phi()*50/PI;
//				p1=pt->GetBinContent(etabin1+1,phibin1+1)/pt_den->GetBinContent(etabin1+1,phibin1+1);
//				p2=pt->GetBinContent(etabin2+1,phibin2+1)/pt_den->GetBinContent(etabin2+1,phibin2+1);
//				histos.fill(HIST("den"),etabin1*100+phibin1+0.5,etabin2*100+phibin2+0.5,histos.get<TH2D>(HIST("rho1"))->GetBinContent(etabin1,phibin1)*histos.get<TH2D>(HIST("rho1"))->GetBinContent(etabin2,phibin2));
//				histos.fill(HIST("num"),etabin1*100+phibin1+0.5,etabin2*100+phibin2+0.5,(track1.pt()-p1)*(track2.pt()-p2));
			}
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<p2>(cfgc)
    };
}
