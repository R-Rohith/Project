#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct r2p24ch{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);//&&(nabs(aod::collision::posZ) > 5.f);
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;
	Filter globalfilter=requireGlobalTrackInFilter();

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM","rho_2 for +ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP","pt for +ve",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM","pt for -ve",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM","ptpt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM","ptn for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM","npt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
		for(auto track1:tracks)
		{
			histos.fill(HIST("pt1"),track1.pt());
			histos.fill(HIST("bb2"),track1.p(),track1.beta());
			histos.fill(HIST("bb1"),track1.p(),track1.tpcSignal());
			histos.fill(HIST("bb4"),track1.p(),track1.beta());
			histos.fill(HIST("bb3"),track1.p(),track1.tpcSignal());

			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};


struct r2p24pi{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM","rho_2 for +ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP","pt for +ve",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM","pt for -ve",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM","ptpt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM","ptn for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM","npt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID--------------------------------------------
			histos.fill(HIST("pt1"),track1.pt());
			histos.fill(HIST("bb1_b"),track1.p(),track1.tpcNSigmaPi());
			histos.fill(HIST("bb2_b"),track1.p(),track1.tofNSigmaPi());
			histos.fill(HIST("bb2"),track1.p(),track1.beta());
			histos.fill(HIST("bb1"),track1.p(),track1.tpcSignal());
			if(track1.hasTOF())
			histos.fill(HIST("pt2"),track1.pt());
//-----------------------------PION---------------------------------------------------
			tpccut=2.5;
			tofcut=2.5;
			if(fabs(track1.tpcNSigmaPi())>=tpccut) continue;	
			if(track1.pt()>=0.6)
			{
				if(track1.hasTOF())
				{
					if(fabs(track1.tofNSigmaPi())>=tofcut) continue;
				}
				else 
				continue;
			}
//--------------------------------------------------------------------*/
//			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.p(),track1.beta());
			histos.fill(HIST("bb1_a"),track1.p(),track1.tpcNSigmaPi());
			histos.fill(HIST("proj1"),track1.tpcNSigmaPi());
			histos.fill(HIST("bb2_a"),track1.p(),track1.tofNSigmaPi());
			histos.fill(HIST("proj2"),track1.tofNSigmaPi());
			histos.fill(HIST("bb3"),track1.p(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
//				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID---------------------------------------
//-----------------------------PION-----------------------------------------------------------
				tpccut=2.5;
				tofcut=2.5;
				if(fabs(track2.tpcNSigmaPi())>=tpccut) continue;	
				if(track2.pt()>=0.6)
				{
					if(track2.hasTOF())
					{
						if(fabs(track2.tofNSigmaPi())>=tofcut) continue;
					}
					else 
					continue;
				}		
//--------------------------------------------------------------------------*/
//				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------*/
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};

struct r2p24ka{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM","rho_2 for +ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP","pt for +ve",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM","pt for -ve",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM","ptpt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM","ptn for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM","npt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCKa,aod::pidTOFKa,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID--------------------------------------------
			histos.fill(HIST("pt1"),track1.pt());
			histos.fill(HIST("bb1_b"),track1.pt(),track1.tpcNSigmaKa());
			histos.fill(HIST("bb2_b"),track1.pt(),track1.tofNSigmaKa());
			histos.fill(HIST("bb2"),track1.p(),track1.beta());
			histos.fill(HIST("bb1"),track1.p(),track1.tpcSignal());
//--------------------KAON------------------------------------------------------
			tpccut=2;
			tofcut=2;
			if(track1.pt()<0.6)
			{
				if(track1.pt()<0.45)tpccut=2;
				else if(track1.pt()<0.55)tpccut=1;
				else tpccut=0.6;
				if(fabs(track1.tpcNSigmaKa())>tpccut) continue;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaKa())>tpccut)||(fabs(track1.tofNSigmaKa())>tofcut)) continue;
			}
			else
			continue;
//-----------------------------------------------------------------------------------------------------------*/				
			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.p(),track1.beta());
			histos.fill(HIST("bb1_a"),track1.pt(),track1.tpcNSigmaKa());
			histos.fill(HIST("proj1"),track1.tpcNSigmaKa());
			histos.fill(HIST("bb2_a"),track1.pt(),track1.tofNSigmaKa());
			histos.fill(HIST("proj2"),track1.tofNSigmaKa());
			histos.fill(HIST("bb3"),track1.p(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID---------------------------------------
//----------------------------Kaon------------------------------------------------------------
				tpccut=2;
				tofcut=2;
				if(track2.pt()<0.6)
				{
					if(track2.pt()<0.45)tpccut=2;
					else if(track2.pt()<0.55)tpccut=1;
					else tpccut=0.6;
					if(fabs(track2.tpcNSigmaKa())>tpccut) continue;
				}
				else if(track2.hasTOF())
				{
					if((fabs(track2.tpcNSigmaKa())>tpccut)||(fabs(track2.tofNSigmaKa())>tofcut)) continue;
				}
				else
				continue;
//-----------------------------------------------------------------------------------------------*/	
				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};

struct r2p24pr{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP","pt for +ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1d_n1_ptM","pt for -ve",kTH1D,{{100,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP","rho_1 for +ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM","rho_1 for -ve particle",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM","rho_2 for +ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP","pt for +ve",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM","pt for -ve",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM","ptpt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM","ptn for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM","npt for +ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"p"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"p"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPr,aod::pidTOFPr,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID--------------------------------------------
			histos.fill(HIST("pt1"),track1.pt());
			histos.fill(HIST("bb1_b"),track1.pt(),track1.tpcNSigmaPr());
			histos.fill(HIST("bb2_b"),track1.pt(),track1.tofNSigmaPr());
			histos.fill(HIST("bb2"),track1.p(),track1.beta());
			histos.fill(HIST("bb1"),track1.p(),track1.tpcSignal());
//-------------------Proton------------------------------------
			tpccut=2.2;
			tofcut=2;
//			if(track1.tpcNClsCrossedRows()<70) continue;
			if(track1.pt()<1.1)
			{
				if(track1.pt()<0.85)tpccut=2.2;
				else tpccut=1;
				if(fabs(track1.tpcNSigmaPr())>tpccut) continue;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaPr())>tpccut)||(fabs(track1.tofNSigmaPr())>tofcut)) continue;
			}
			else
			continue;

//--------------------------------------*/
//			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.p(),track1.beta());
			histos.fill(HIST("bb1_a"),track1.pt(),track1.tpcNSigmaPr());
			histos.fill(HIST("proj1"),track1.tpcNSigmaPr());
			histos.fill(HIST("bb2_a"),track1.pt(),track1.tofNSigmaPr());
			histos.fill(HIST("proj2"),track1.tofNSigmaPr());
			histos.fill(HIST("bb3"),track1.p(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID---------------------------------------
//-------------------Proton------------------------------------
				tpccut=2.2;
				tofcut=2;
//				if(track2.tpcNClsCrossedRows()<70) continue;
				if(track2.pt()<1.1)
				{
					if(track2.pt()<0.85)tpccut=2.2;
					else tpccut=1;
					if(fabs(track2.tpcNSigmaPr())>tpccut) continue;
				}
				else if(track2.hasTOF())
				{
					if((fabs(track2.tpcNSigmaPr())>tpccut)||(fabs(track2.tofNSigmaPr())>tofcut)) continue;
				}
				else
				continue;

//--------------------------------------*/								
								
//				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};

struct crossr2p24pik{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP1","pt for +ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM1","pt for -ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptP2","pt for +ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM2","pt for -ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP1","rho_1 for +ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM1","rho_1 for -ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiP2","rho_1 for +ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM2","rho_1 for -ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12","rho_2 for +ve_1 -ve_2 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21","rho_2 for +ve_2 -ve_1 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP1","pt for +ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM1","pt for -ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiP2","pt for +ve_2",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM2","pt for -ve_2",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12","ptpt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21","ptpt for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12","ptn for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21","ptn for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTPCPr,aod::pidTOFPr,aod::pidTPCKa,aod::pidTOFKa,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			bool flag;
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID2--------------------------------------------
			flag=true;
//--------------------KAON------------------------------------------------------
			tpccut=2;
			tofcut=2;
			if(track1.pt()<0.6)
			{
				if(track1.pt()<0.45)tpccut=2;
				else if(track1.pt()<0.55)tpccut=1;
				else tpccut=0.6;
				if(fabs(track1.tpcNSigmaKa())>tpccut)flag=false;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaKa())>tpccut)||(fabs(track1.tofNSigmaKa())>tofcut))flag=false;
			}
			else
			flag=false;
//-----------------------------------------------------------------------------------------------------------*/				
//			if(fabs(track1.tpcNSigmaEl())<1.2) flag=false;
			if(flag==true)
			{
				histos.fill(HIST("bb2"),track1.pt(),track1.beta());
				histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
				histos.fill(HIST("phi"),track1.phi());	
				histos.fill(HIST("eta"),track1.eta());
				if(track1.sign()==1)
				{	
					histos.fill(HIST("h2d_n1_etaPhiP2"),track1.eta(),track1.phi());
					histos.fill(HIST("h1d_n1_ptP2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
					histos.fill(HIST("h2d_pt_etaPhiP2"),track1.eta(),track1.phi(),track1.pt());
				}
				else
				{
					histos.fill(HIST("h2d_n1_etaPhiM2"),track1.eta(),track1.phi());
					histos.fill(HIST("h2d_pt_etaPhiM2"),track1.eta(),track1.phi(),track1.pt());
					histos.fill(HIST("h1d_n1_ptM2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				}
			}






//---------------------------------------PID1--------------------------------------------
//-----------------------------PION---------------------------------------------------
			tpccut=2.5;
			tofcut=2.5;
			if(fabs(track1.tpcNSigmaPi())>=tpccut) continue;	
			if(track1.pt()>=0.6)
			{
				if(track1.hasTOF())
				{
					if(fabs(track1.tofNSigmaPi())>=tofcut) continue;
				}
				else 
				continue;
			}
//---------------------------------------------------------------------------------*/
			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.pt(),track1.beta());
			histos.fill(HIST("bb3"),track1.pt(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP1"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP1"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM1"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM1"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID2---------------------------------------
//----------------------------Kaon------------------------------------------------------------
				tpccut=2;
				tofcut=2;
				if(track2.pt()<0.6)
				{
					if(track2.pt()<0.45)tpccut=2;
					else if(track2.pt()<0.55)tpccut=1;
					else tpccut=0.6;
					if(fabs(track2.tpcNSigmaKa())>tpccut) continue;
				}
				else if(track2.hasTOF())
				{
					if((fabs(track2.tpcNSigmaKa())>tpccut)||(fabs(track2.tofNSigmaKa())>tofcut)) continue;
				}
				else
				continue;
//-----------------------------------------------------------------------------------------------*/	
								
//				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						if(track1.sign()==1)
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
						else
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};

struct crossr2p24pip{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP1","pt for +ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM1","pt for -ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptP2","pt for +ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM2","pt for -ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP1","rho_1 for +ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM1","rho_1 for -ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiP2","rho_1 for +ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM2","rho_1 for -ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12","rho_2 for +ve_1 -ve_2 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21","rho_2 for +ve_2 -ve_1 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP1","pt for +ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM1","pt for -ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiP2","pt for +ve_2",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM2","pt for -ve_2",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12","ptpt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21","ptpt for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12","ptn for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21","ptn for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTPCPr,aod::pidTOFPr,aod::pidTPCKa,aod::pidTOFKa,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			bool flag;
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID2--------------------------------------------
			flag=true;
//-------------------Proton------------------------------------
			tpccut=2.2;
			tofcut=2;
			if(track1.pt()<1.1)
			{
				if(track1.pt()<0.85)tpccut=2.2;
				else tpccut=1;
				if(fabs(track1.tpcNSigmaPr())>tpccut) flag=false;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaPr())>tpccut)||(fabs(track1.tofNSigmaPr())>tofcut)) flag=false;
			}
			else
			flag=false;

//--------------------------------------*/
			if(fabs(track1.tpcNSigmaEl())<1.2) flag=false;
			if(flag==true)
			{
				histos.fill(HIST("bb2"),track1.pt(),track1.beta());
				histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
				histos.fill(HIST("phi"),track1.phi());	
				histos.fill(HIST("eta"),track1.eta());
				if(track1.sign()==1)
				{	
					histos.fill(HIST("h2d_n1_etaPhiP2"),track1.eta(),track1.phi());
					histos.fill(HIST("h1d_n1_ptP2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
					histos.fill(HIST("h2d_pt_etaPhiP2"),track1.eta(),track1.phi(),track1.pt());
				}
				else
				{
					histos.fill(HIST("h2d_n1_etaPhiM2"),track1.eta(),track1.phi());
					histos.fill(HIST("h2d_pt_etaPhiM2"),track1.eta(),track1.phi(),track1.pt());
					histos.fill(HIST("h1d_n1_ptM2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				}
			}






//---------------------------------------PID1--------------------------------------------
//-----------------------------PION---------------------------------------------------
			tpccut=2.5;
			tofcut=2.5;
			if(fabs(track1.tpcNSigmaPi())>=tpccut) continue;	
			if(track1.pt()>=0.6)
			{
				if(track1.hasTOF())
				{
					if(fabs(track1.tofNSigmaPi())>=tofcut) continue;
				}
				else 
				continue;
			}
//---------------------------------------------------------------------------------*/
			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.pt(),track1.beta());
			histos.fill(HIST("bb3"),track1.pt(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP1"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP1"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM1"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM1"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID2---------------------------------------
//-------------------Proton------------------------------------
				tpccut=2.2;
				tofcut=2;
				if(track2.pt()<1.1)
				{
					if(track2.pt()<0.85)tpccut=2.2;
					else tpccut=1;
					if(fabs(track2.tpcNSigmaPr())>tpccut) continue;
				}
				else if(track2.hasTOF())
				{
					if((fabs(track2.tpcNSigmaPr())>tpccut)||(fabs(track2.tofNSigmaPr())>tofcut)) continue;
				}
				else
				continue;

//--------------------------------------*/								
								
				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						if(track1.sign()==1)
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
						else
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};

struct crossr2p24pk{
	Filter col = aod::evsel::sel8 == true;
	Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
	Filter ptfilter=aod::track::pt>0.2f&&aod::track::pt<2.0f;		//till 0.5 no delta peak....peak seems to begin from 0.25.
	Filter globalfilter=requireGlobalTrackInFilter();
	//should probably also add an eta filter...

	HistogramRegistry histos{"R2P2",{},OutputObjHandlingPolicy::AnalysisObject};
	void init(InitContext const &)
	{
		const AxisSpec phi{36,0,2.0*constants::math::PI,"phi"},eta{24,-0.8,0.8,"eta"},etaphi1{864,0,864,"etaphi1"},etaphi2{864,0,864,"etaphi2"};
		histos.add("phi","phi",kTH1D,{phi});
		histos.add("eta","eta",kTH1D,{eta});
		histos.add("test","tst",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("h1d_n1_ptP1","pt for +ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM1","pt for -ve_1",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptP2","pt for +ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1d_n1_ptM2","pt for -ve_2",kTH1D,{{30,0,6,"pt"}});
		histos.add("h1i_n1_multPM","multiplicity",kTH1I,{{200,0,200,"mult"}});		
		histos.add("h2d_n1_etaPhiP1","rho_1 for +ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM1","rho_1 for -ve particle1",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiP2","rho_1 for +ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n1_etaPhiM2","rho_1 for -ve particle2",kTH2D,{eta,phi});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PP","rho_2 for +ve+ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12","rho_2 for +ve_1 -ve_2 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21","rho_2 for +ve_2 -ve_1 particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_n2_eta1Phi1Eta2Phi2MM","rho_2 for -ve-ve particle",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_pt_etaPhiP1","pt for +ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM1","pt for -ve_1",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiP2","pt for +ve_2",kTH2D,{eta,phi});
		histos.add("h2d_pt_etaPhiM2","pt for -ve_2",kTH2D,{eta,phi});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12","ptpt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21","ptpt for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP","ptpt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM","ptpt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12","ptn for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21","ptn for +ve_2 -ve_1",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP","ptn for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM","ptn for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21","npt for +ve_1 -ve_2",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2PP","npt for +ve+ve",kTH2D,{etaphi1,etaphi2});
		histos.add("h2d_npt_eta1Phi1Eta2Phi2MM","npt for -ve-ve",kTH2D,{etaphi1,etaphi2});
		histos.add("pt1","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("pt2","pt",kTH1D,{{100,0,5,"pt"}});
		histos.add("bb1","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb3","tpc",kTH2D,{{100,0,5,"pt"},{200,0,200,"de/dx"}});
		histos.add("bb2","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb4","tof",kTH2D,{{100,0,5,"pt"},{50,0.2,1.2,"beta"}});
		histos.add("bb1_b","tpcnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb1_a","tpcnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_b","tofnsigma_before",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("bb2_a","tofnsigma_after",kTH2D,{{100,0,5,"pt"},{25,-5,5,"nsigma"}});
		histos.add("proj1","tpc_proj",kTH1D,{{25,-5,5,"nsigma"}});
		histos.add("proj2","tof_proj",kTH1D,{{25,-5,5,"nsigma"}});
	}
	void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection,aod::pidTPCPi,aod::pidTOFPi,aod::pidTPCPr,aod::pidTOFPr,aod::pidTPCKa,aod::pidTOFKa,aod::pidTPCEl,aod::pidTOFbeta,aod::TracksExtra>> const& tracks)
	{
			bool flag;
			int mult=0;
			int etabin1,etabin2,phibin1,phibin2;
			float tpccut,tofcut;
		for(auto track1:tracks)
		{
//---------------------------------------PID2--------------------------------------------
			flag=true;
//--------------------KAON------------------------------------------------------
			tpccut=2;
			tofcut=2;
			if(track1.pt()<0.6)
			{
				if(track1.pt()<0.45)tpccut=2;
				else if(track1.pt()<0.55)tpccut=1;
				else tpccut=0.6;
				if(fabs(track1.tpcNSigmaKa())>tpccut)flag=false;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaKa())>tpccut)||(fabs(track1.tofNSigmaKa())>tofcut))flag=false;
			}
			else
			flag=false;
//-----------------------------------------------------------------------------------------------------------*/				
//--------------------------------------*/
			if(fabs(track1.tpcNSigmaEl())<1.2) flag=false;
			if(flag==true)
			{
				histos.fill(HIST("bb2"),track1.pt(),track1.beta());
				histos.fill(HIST("bb1"),track1.pt(),track1.tpcSignal());
				histos.fill(HIST("phi"),track1.phi());	
				histos.fill(HIST("eta"),track1.eta());
				if(track1.sign()==1)
				{	
					histos.fill(HIST("h2d_n1_etaPhiP2"),track1.eta(),track1.phi());
					histos.fill(HIST("h1d_n1_ptP2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
					histos.fill(HIST("h2d_pt_etaPhiP2"),track1.eta(),track1.phi(),track1.pt());
				}
				else
				{
					histos.fill(HIST("h2d_n1_etaPhiM2"),track1.eta(),track1.phi());
					histos.fill(HIST("h2d_pt_etaPhiM2"),track1.eta(),track1.phi(),track1.pt());
					histos.fill(HIST("h1d_n1_ptM2"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				}
			}






//---------------------------------------PID1--------------------------------------------
//-----------------------------------------------------------------------------------------------------------*/				
//-------------------Proton------------------------------------
			tpccut=2.2;
			tofcut=2;
//			if(track1.tpcNClsCrossedRows()<70) continue;
			if(track1.pt()<1.1)
			{
				if(track1.pt()<0.85)tpccut=2.2;
				else tpccut=1;
				if(fabs(track1.tpcNSigmaPr())>tpccut) continue;
			}
			else if(track1.hasTOF())
			{
				if((fabs(track1.tpcNSigmaPr())>tpccut)||(fabs(track1.tofNSigmaPr())>tofcut)) continue;
			}
			else
			continue;

//--------------------------------------*/
			if(fabs(track1.tpcNSigmaEl())<1.2) continue;

			histos.fill(HIST("bb4"),track1.pt(),track1.beta());
			histos.fill(HIST("bb3"),track1.pt(),track1.tpcSignal());

//--------------------------------------------------------END----------------------------------------------------------*/
			histos.fill(HIST("phi"),track1.phi());
			histos.fill(HIST("eta"),track1.eta());
			if(track1.sign()==1)
			{
				histos.fill(HIST("h2d_n1_etaPhiP1"),track1.eta(),track1.phi());
				histos.fill(HIST("h1d_n1_ptP1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
				histos.fill(HIST("h2d_pt_etaPhiP1"),track1.eta(),track1.phi(),track1.pt());
			}
			else
			{
				histos.fill(HIST("h2d_n1_etaPhiM1"),track1.eta(),track1.phi());
				histos.fill(HIST("h2d_pt_etaPhiM1"),track1.eta(),track1.phi(),track1.pt());
				histos.fill(HIST("h1d_n1_ptM1"),track1.pt(),1.0/(2.0*constants::math::PI*track1.pt()));
			}
				if(track1.eta()>=(0.8)||track1.eta()<=-0.8)continue;
				etabin1=(track1.eta()+0.8)*15;
				phibin1=36*track1.phi()/(2*constants::math::PI);
				mult++;
			for(auto track2:tracks)
			{
				if(track1.index()==track2.index())continue;
//---------------------------------------------------PID2---------------------------------------
//-------------------------------------------------------------------------------------------*/
//----------------------------Kaon------------------------------------------------------------
				tpccut=2;
				tofcut=2;
				if(track2.pt()<0.6)
				{
					if(track2.pt()<0.45)tpccut=2;
					else if(track2.pt()<0.55)tpccut=1;
					else tpccut=0.6;
					if(fabs(track2.tpcNSigmaKa())>tpccut) continue;
				}
				else if(track2.hasTOF())
				{
					if((fabs(track2.tpcNSigmaKa())>tpccut)||(fabs(track2.tofNSigmaKa())>tofcut)) continue;
				}
				else
				continue;
//-----------------------------------------------------------------------------------------------*/	
				if(fabs(track2.tpcNSigmaEl())<1.2) continue;

//------------------------------------------PID----END----------------------------------------
				if(track2.eta()>=(0.8)||track2.eta()<=-0.8)continue;
				etabin2=(track2.eta()+0.8)*15;
				phibin2=36*track2.phi()/(2*constants::math::PI);
				if((etabin1>=0)&&(etabin2>=0)&&(phibin1>=0)&&(phibin2>=0)&&(etabin1<24)&&(etabin2<24)&&(phibin1<36)&&(phibin2<36))
				{
					if(track1.sign()*track2.sign()==-1)
					{
						if(track1.sign()==1)
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
						else
						{
							histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
							histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
							histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
							histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
						}
					}
					else if(track1.sign()==1)
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
					else
					{
						histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5);
						histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track2.pt());
						histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt());
						histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"),36*etabin1+phibin1+0.5,36*etabin2+phibin2+0.5,track1.pt()*track2.pt());
					}
				}
			}
		}
		histos.fill(HIST("h1i_n1_multPM"),mult);
		}

};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<r2p24ch>(cfgc),	
    adaptAnalysisTask<r2p24pi>(cfgc),
    adaptAnalysisTask<r2p24ka>(cfgc),
    adaptAnalysisTask<r2p24pr>(cfgc),
/*    adaptAnalysisTask<crossr2p24pik>(cfgc),
    adaptAnalysisTask<crossr2p24pip>(cfgc),
    adaptAnalysisTask<crossr2p24pk>(cfgc),*/
    };
}

/*
Electron tpccut to be changed to 2
ADD stuff like tpcNClsCrossedRows...

*/
