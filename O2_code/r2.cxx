#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct r2
{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>1.0f&&aod::track::pt<3.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();

	HistogramRegistry r2histos{"R2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const&)
	{ 
		const AxisSpec phi{36,0,2.0*PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		r2histos.add("phi","phi",kTH1D,{phi});
		r2histos.add("eta","eta",kTH1D,{eta});
		r2histos.add("test","tst",kTH1D,{{1000,-100,900,"test"}});
		r2histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{30,0,6,"pt"}});
		r2histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{30,0,6,"pt"}});
		r2histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
//		r2histos.add("h1i_n1_multPM2","multiplicity",kTH1I,{{200,0,200,"mult"}});
		r2histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{eta,phi});
		r2histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{eta,phi});
		r2histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_n2_eta1Phi1Eta2Phi2PM","rho_2 for +ve-ve particle",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_pt_etaPhiP","pt for +ve",kTH2D,{eta,phi});
		r2histos.add("h2d_pt_etaPhiM","pt for -ve",kTH2D,{eta,phi});
		r2histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM","ptpt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM","ptn for +ve-ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_npt_eta1Phi1Eta2Phi2PM","npt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		r2histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});

	}
//	void process(/*soa::Filtered<aod::Collisions> const& filteredCollisions*/aod::Collision const& collision,aod::TracksIU const& tracks)
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
	{
//		int i=0;
//		for(auto& collision:filteredCollisions){
		int mult=0,etabin1,etabin2,phibin1,phibin2;
//		double phi1,phi2,eta1,eta2;
//		r2histos.fill(HIST("h1i_n1_multPM"),filteredCollisions.numContrib());
		for(auto track1:tracks)
		{
			mult++;etabin1=track1.eta();
			//keep data within range
			r2histos.fill(HIST("phi"),track1.phi());
			r2histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				r2histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				r2histos.fill(HIST("h1d_n1_ptP"),track1.pt());
				r2histos.fill(HIST("h2d_pt_etaPhiP"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				r2histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				r2histos.fill(HIST("h2d_pt_etaPhiM"),track1.eta(),track1.phi(),track1.pt());
				r2histos.fill(HIST("h1d_n1_ptM"),track1.pt());
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*PI);

			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//				eta1=track1.eta();
//				eta2=track2.eta();
//				phi1=track1.phi();
//				phi2=track2.phi();
//				if(phi1<0)phi1+=2*PI;
//				if(phi2<0)phi2+=2*PI;
//				etabin1=(track1.eta()+0.8)*15;
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
//				phibin1=36*track1.phi()/(2*PI);
				phibin2=36*track2.phi()/(2*PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
//					r2histos.fill(HIST("test"),etabin1);
					if(track1.sign()*track2.sign()==-1)
					{
						r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						r2histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						r2histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						r2histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else if(track1.sign()==1)
					{
						r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						r2histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						r2histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						r2histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						r2histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						r2histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						r2histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		r2histos.fill(HIST("h1i_n1_multPM"),mult);
/*		for(auto[t1,t2]:combinations(o2::soa::CombinationsFullIndexPolicy(tracks,tracks)))
		{
			if(t1.index()==t2.index())continue;
			if(t1.eta()>=(0.8)||t1.eta()<=-0.8)continue;
			etabin1=(t1.eta()+0.8)*15;
			phibin1=36*t1.phi()/(2*PI);
			if(t2.eta()>=(0.8)||t2.eta()<=-0.8)continue;
			etabin2=(t2.eta()+0.8)*15;
			phibin2=36*t2.phi()/(2*PI);
			if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
				if(t1.sign()*t2.sign()==-1)
				r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
				else if(t1.sign()==1)
				r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
				else
				r2histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
			}
		}*/
//		if(i>0)break;
		}
//	}
};
struct phi{
	HistogramRegistry histos{"tst",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const&)
	{
		histos.add("dphi","dphi",kTH1D,{{100,-0.5*PI,1.5*PI,"dphi"}});
	}
	void process(aod::Collision const& collision,aod::TracksIU const& track)
	{
		double dphi;
		for(auto track1:track)
		{
			for(auto track2:track)
			{
				if(track1.index()>=track2.index())continue;
				dphi=track1.phi()-track2.phi();
				if(dphi<(-0.5*PI))dphi+=2*PI;
				if(dphi>(1.5*PI))dphi-=2*PI;
				histos.fill(HIST("dphi"),dphi);
			}
		}
	}
};
/*
struct tst{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt<2.f;
	Filter globalfilter=requireGlobalTrackInFilter();
	
	HistogramRegistry tsthistos{"Test",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const&)
	{
//		tsthistos.add("posZ","posZ",kTH1D,{{100,-100,100,"pos_z"}});
		tsthistos.add("pt","pt",kTH1D,{{100,0,10,"pt"}});
		tsthistos.add("phi","phi",kTH1D,{{24,0,2.0*PI,"phi"}});
		tsthistos.add("eta","eta",kTH1D,{{36,-0.8,0.8,"eta"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>> const& filteredCollisions, soa::Filtered<soa::Join<aod::TracksIU, aod::TrackSelection>> const& tracks)	{
//		if(collision.posZ()>4){
//		tsthistos.fill(HIST("posZ"),collision.posZ());
		for(auto track:tracks)
		{
			tsthistos.fill(HIST("pt"),track.pt());
//			if(track.pt()>1.0)continue;
			tsthistos.fill(HIST("phi"),track.phi());
			tsthistos.fill(HIST("eta"),track.eta());
		}
//		}
	}
};*/
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<r2>(cfgc),
    adaptAnalysisTask<phi>(cfgc),
//    adaptAnalysisTask<tst>(cfgc),
    };
}
