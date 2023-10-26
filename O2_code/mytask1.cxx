#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include <iostream>
#include<algorithm>

using namespace o2::framework;
using namespace o2;
using namespace o2::framework::expressions;
using namespace o2::aod;

using myFilteredTracks = soa::Filtered<aod::TracksIU>;
/*
struct myfirst_task
{

	SliceCache cache;
	HistogramRegistry myhisto{"Hist_name",{},OutputObjHandlingPolicy::AnalysisObject};
	Filter etaFilter = nabs(track::eta) < 0.8f;
	Partition<myFilteredTracks> ptPart1 = track::pt > 1.f;
	Partition<myFilteredTracks> ptPart2 = track::pt < 1.f;
	
	void init(InitContext const&)
	{
		//define histograms & axes
		const AxisSpec ptaxis{100,-0.5,4,"pt"},dedxaxis{100,0,100,"de/dx"};
		myhisto.add("bloch","Bethe-Bloch graph", kTH1F,{ptaxis/*,dedxaxis*]/});
	}
	void process(aod::Collision const& collision, myFilteredTracks const& tracks)//subscribed data as arguments
	{
		auto ptlowtracks = ptPart1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
		auto pthightracks = ptPart2->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
		for(auto& ptPart : ptlowtracks)
		{
//			if(track.pt()>1)
			myhisto.fill(HIST("bloch"), ptPart.pt()/*,track.energy()*]/);  
		}
	}
};
*/
struct myfirst_task
{	
	HistogramRegistry myhistos{"duh",{},OutputObjHandlingPolicy::AnalysisObject};
	Configurable<float> up{"up", 0.5,"brief"},low{"low",0.2,"another brief"},diff{"diff",1,"diff in pair pt"};
	void init(InitContext const&)
	{
		const AxisSpec phiaxis{100,-3,5,"delta phi"};
		const AxisSpec etaaxis{500,-4,5,"delta eta"};
		myhistos.add("eh","yeild vs delta phi",kTH1F,{phiaxis/*,etaaxis*/});
		myhistos.add("1h","yeild vs phi",kTH1F,{phiaxis});
		myhistos.add("2h","yeild vs delta eta",kTH1F,{etaaxis});

	}
	void process(aod::Collision const& collision,aod::TracksIU const& tracks)
	{
		long int i=0/*,no=0*/;
		float dphi,deta;
//		for(auto coll: collision){
		for (auto track1 : tracks)
		{
			if((track1.pt()>low)&&(track1.pt()<up)){ //Technically better to just use filter.
			for(auto track2 : tracks)
			{
				if((track1.index() >=track2.index())/*||(fabs(track1.pt()-track2.pt())>diff)*/)
					continue;
				if((track2.pt()>low)&&(track2.pt()<up))
				{
					dphi=track1.phi()-track2.phi();
					if (dphi > 1.5f * PI) {
					dphi -= TwoPI;
					}
					if (dphi < -PIHalf) {
					dphi += TwoPI;
					}
					deta=track1.eta()-track2.eta();
					if (deta > 1.5f * PI) {
					deta -= TwoPI;
					}
					if (deta < -PIHalf) {
					deta += TwoPI;
					}
					myhistos.fill(HIST("eh"),dphi/*,deta*/);
					myhistos.fill(HIST("2h"),deta);
				}
				/*
				Still the plot is not being developed properly. It's been verified<?> that the phi distribution is normal. In such case the current code will give the current plot, which means
				to get the plot from the paper we need to add some other condition on the tracks. My first thought, it might work if only consider tracks that can be correlated.
				*/
			}
			}
					myhistos.fill(HIST("1h"),track1.phi());
//			if(i>1000)
//			
			break;
			i++;
//			std::cout<<" "<<track1.collisionId();
//			no++;
		} //}
//		std::cout<<" "<<no;}
//			myhistos.fill(HIST("1h"),-2);	
	}
};
struct tasc
{
	HistogramRegistry myhistos{"duh",{},OutputObjHandlingPolicy::AnalysisObject};
	Configurable<int> Id{"Id", 0,"collid"};
	void init(InitContext const&)
	{
		const AxisSpec phiaxis{100,-7,7,"delta phi"};
		myhistos.add("eh","yeild vs delta phi",kTH1F,{phiaxis});
	}
	void process(aod::Collisions const& collision,aod::TracksIU const& tracks)
	{
		for(auto coll:collision)
		{
			int i=0,id;
			for(auto track:tracks)
			{
				if(i==Id)id=track.collisionId();
				if(id!=track.collisionId())break;
				myhistos.fill(HIST("eh"),track.phi());
			}
			break;
		}
	}
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
	return WorkflowSpec
	{
		adaptAnalysisTask<myfirst_task>(cfgc),
		adaptAnalysisTask<tasc>(cfgc),
		//more tasks if needed
	};
}

