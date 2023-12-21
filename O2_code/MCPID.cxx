#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Kaon{
//	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=(aod::track::pt>0.2f)&&(aod::track::pt<2.0f);		//till 0.5 no delta peak....peak seems to begin from 0.25.
	//0.2 fot tpc, 0.35 for tof
	Filter globalfilter=requireGlobalTrackInFilter();
	HistogramRegistry histos{"hist",{},OutputObjHandlingPolicy::AnalysisObject};
	
	void init(InitContext const&)
	{
		histos.add("p","momentum",kTH1D,{{100,0,10,"p"}});
		histos.add("eta","#eta",kTH1D,{{100,-1,1,"eta"}});
		histos.add("phi","#phi",kTH1D,{{200,0,2*constants::math::PI,"#phi"}});
		histos.add("p1","momentum1",kTH1D,{{100,0,10,"p"}});
		histos.add("eta1","#eta1",kTH1D,{{100,-1,1,"eta"}});
		histos.add("phi1","#phi1",kTH1D,{{200,0,2*constants::math::PI,"#phi"}});
		histos.add("p2","momentum2",kTH1D,{{100,0,10,"p"}});
		histos.add("eta2","#eta2",kTH1D,{{100,-1,1,"eta"}});
		histos.add("phi2","#phi2",kTH1D,{{200,0,2*constants::math::PI,"#phi"}});
		histos.add("tpc","TPC",kTH2D,{{100,0,3,"p"},{300,0,300,"${#frac{dE}{dx}}$"}});
		histos.add("tof","TOF",kTH2D,{{100,0,3,"p"},{300,-1,2,"#beta"}});
		histos.add("nsigmatpc","N#sigma_{TPC}",kTH2D,{{100,0,3,"p"},{100,-6.5,6.5,"N#sigma_{TPC}"}});
		histos.add("nsigmatof","N#sigma_{TOF}",kTH2D,{{100,0,3,"p"},{100,-6.5,6.5,"N#sigma_{TOF}"}});
		histos.add("tpc1","TPC",kTH2D,{{100,0,3,"p"},{300,0,300,"#frac{dE}{dx}"}});
		histos.add("tof1","TOF",kTH2D,{{100,0,3,"p"},{300,-1,2,"#beta"}});
		histos.add("nsigmatpc1","N#sigma_{TPC}",kTH2D,{{100,0,3,"p"},{100,-6.5,6.5,"N#sigma_{TPC}"}});
		histos.add("nsigmatof1","N#sigma_{TOF}",kTH2D,{{100,0,3,"p"},{100,-6.5,6.5,"N#sigma_{TOF}"}});
		histos.add("tpc2","TPC",kTH2D,{{100,0,3,"p"},{300,0,300,"#frac{dE}{dx}"}});
		histos.add("tof2","TOF",kTH2D,{{100,0,3,"p"},{300,-1,2,"#beta"}});
		histos.add("pdist","p distribution",kTH1D,{{100,0,3,"p"}});
		histos.add("pdist1","p distribution",kTH1D,{{100,0,3,"p"}});
		histos.add("cnt","counter",kTH1D,{{20,0,20,"category"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidBayesEl,aod::pidBayesPi,aod::pidBayesKa,aod::pidBayesPr,aod::pidTPCEl,aod::pidTPCKa,aod::pidTOFKa,aod::pidTOFbeta,aod::TracksExtra,aod::McTrackLabels>> const& tracks,aod::McParticles const&)
	{
		float tpccut,tofcut;
		for(auto track1: tracks)
		{
			if(!track1.has_mcParticle()){ histos.fill(HIST("cnt"),19.5);continue;}
			histos.fill(HIST("cnt"),0.5);
			if(fabs(track1.mcParticle().pdgCode())==11)histos.fill(HIST("cnt"),1.5);
			if(fabs(track1.mcParticle().pdgCode())==211)histos.fill(HIST("cnt"),2.5);
			if(fabs(track1.mcParticle().pdgCode())==321)histos.fill(HIST("cnt"),3.5);
			if(fabs(track1.mcParticle().pdgCode())==2212)histos.fill(HIST("cnt"),4.5);
			histos.fill(HIST("eta"),track1.eta());
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("p"),track1.p());
			tofcut=2;
			if(track1.p()<0.45)tpccut=2;
			else if(track1.p()<0.55)tpccut=1;
			else if(track1.p()<0.6)tpccut=0.6;
			else tpccut=3;
			histos.fill(HIST("tpc"),track1.p(),track1.tpcSignal());
			histos.fill(HIST("tof"),track1.p(),track1.beta());
			histos.fill(HIST("nsigmatpc"),track1.p(),track1.tpcNSigmaKa());
			histos.fill(HIST("nsigmatof"),track1.p(),track1.tpcNSigmaKa());
			if(track1.p()<0.6)
			{
				if(fabs(track1.tpcNSigmaKa())>tpccut)continue;
				if((fabs(track1.tpcNSigmaEl())<2)&&(fabs(track1.tofNSigmaKa())>tofcut))continue;
			}
			else
			{
				if((fabs(track1.tpcNSigmaKa())>tpccut)||(fabs(track1.tofNSigmaKa())>tofcut))continue;
			}
			histos.fill(HIST("pdist"),track1.p());
			histos.fill(HIST("tpc1"),track1.p(),track1.tpcSignal());
			histos.fill(HIST("tof1"),track1.p(),track1.beta());
			histos.fill(HIST("nsigmatpc1"),track1.p(),track1.tpcNSigmaKa());
			histos.fill(HIST("nsigmatof1"),track1.p(),track1.tofNSigmaKa());
			histos.fill(HIST("eta1"),track1.eta());
			histos.fill(HIST("phi1"),track1.phi());
			histos.fill(HIST("p1"),track1.p());
			histos.fill(HIST("cnt"),6.5);
			if(fabs(track1.mcParticle().pdgCode())==11)histos.fill(HIST("cnt"),7.5);
			if(fabs(track1.mcParticle().pdgCode())==211)histos.fill(HIST("cnt"),8.5);
			if(fabs(track1.mcParticle().pdgCode())==321)histos.fill(HIST("cnt"),9.5);
			if(fabs(track1.mcParticle().pdgCode())==2212)histos.fill(HIST("cnt"),10.5);
		}
		float tmp;
		int pid;
		for(auto track1: tracks)
		{
			if(!track1.has_mcParticle())continue;
			tmp=0;pid=0;
			if(track1.bayesEl()>=tmp)
			{
				if(track1.bayesEl()==tmp)pid=pid*10+1;
				else pid=1;
				tmp=track1.bayesEl();
			}
			if(track1.bayesPi()>=tmp)
			{
				if(track1.bayesPi()==tmp)pid=pid*10+2;
				else pid=2;
				tmp=track1.bayesPi();
			}
			if(track1.bayesKa()>=tmp)
			{
				if(track1.bayesKa()==tmp)pid=pid*10+3;
				else pid=3;
				tmp=track1.bayesKa();
			}
			if(track1.bayesPr()>=tmp)
			{
				if(track1.bayesPr()==tmp)pid=pid*10+4;
				else pid=4;
				tmp=track1.bayesPr();
			}
			if(pid!=3)continue;
			histos.fill(HIST("pdist1"),track1.p());
			histos.fill(HIST("tpc2"),track1.p(),track1.tpcSignal());
			histos.fill(HIST("tof2"),track1.p(),track1.beta());
			histos.fill(HIST("eta2"),track1.eta());
			histos.fill(HIST("phi2"),track1.phi());
			histos.fill(HIST("p2"),track1.p());
			histos.fill(HIST("cnt"),12.5);
			if(fabs(track1.mcParticle().pdgCode())==11)histos.fill(HIST("cnt"),13.5);
			if(fabs(track1.mcParticle().pdgCode())==211)histos.fill(HIST("cnt"),14.5);
			if(fabs(track1.mcParticle().pdgCode())==321)histos.fill(HIST("cnt"),15.5);
			if(fabs(track1.mcParticle().pdgCode())==2212)histos.fill(HIST("cnt"),16.5);
		}
	}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return {
	adaptAnalysisTask<Kaon>(cfgc)
	};
}

/*
proton=2212
kaon=321
pion=211
e=11
*/
