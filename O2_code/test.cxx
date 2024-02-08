/*
float tpcnsigmaarray[3]={fabs(track.tpcNSigmaPi()),fabs(track.tpcNSigmaKa()),fabs(track.tpcNSigmaPr())},
      tofnsigmaarray[3]={fabs(track.tofNSigmaPi()),fabs(track.tofNSigmaKa()),fabs(track.tofNSigmaPr())};

int p1,p2;
float nsigp1,nsigp2;
nsigp1=6;nsigp2=6;
for(int i=0;i<3;i++)
{
  if(nsigp1>tpcnsigmaarray[i])
  {
    nsigp2=nsigp1;
    p2=p1;
    p1=i;
    nsigp1=tpcnsigmaarray[i];
  }
  else if(nsigp2>tpcnsigmaarray[i])
  {
    p2=i;
    nsigp2=tpcnsigmaarray[i];
  }
}

if(nsigp1>TPCnsigMax)
  pid=0;
else if(nsigp2>TPCnsigSEP)
{
  pid=p1+1;
}
else if(track.hasTOF()) //Use TOF
{
  pid=p1+1;
  nsigp1=6;nsigp2=6;
  for(int i=0;i<3;i++)
  {
    if(nsigp1>tofnsigmaarray[i])
    {
      nsigp2=nsigp1;
      p2=p1;
      p1=i;
      nsigp1=tofnsigmaarray[i];
    }
    else if(nsigp2>tofnsigmaarray[i])
    {
      p2=i;
      nsigp2=tofnsigmaarray[i];
    }
  }
  if(nsigp1<TOFnsigMax)&&(nsigp2>TOFnsigSEP)
  {
    pid=p1+1;
  }
  else if(ifTOFreject) pid=0; //If we decide to rejct the track when TOF cant identify
  else pid=p1+1;
}
else if(ifTOFreject) pid=0;
else pid=p1+1;





//-------------------------------------------------------------------------------------------------------------------------*/


// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//-----PID Functions-----------------------------------------------------------
float tpccut,tofcut;
template <typename T>
bool NO_PID(T track)
{
  return true;
}
template <typename T>
bool PID_PION(T track)
{
   tpccut = 2.5, tofcut = 2.5;
  if (fabs(track.tpcNSigmaPi()) >= tpccut)
    return false;
  if (track.pt() >= 0.6) {
    if (track.hasTOF()) {
      if (fabs(track.tofNSigmaPi()) >= tofcut)
        return false;
    } else {
      return false;
    }
  }
  return true;
}
template <typename T>
bool PID_KAON(T track)
{
   tpccut = 2, tofcut = 2;
  if (track.pt() < 0.6) {
    if (track.pt() < 0.45)
      tpccut = 2;
    else if (track.pt() < 0.55)
      tpccut = 1;
    else
      tpccut = 0.6;
    if (fabs(track.tpcNSigmaKa()) > tpccut)
      return false;
  } else if (track.hasTOF()) {
    if ((fabs(track.tpcNSigmaKa()) > tpccut) || (fabs(track.tofNSigmaKa()) > tofcut))
      return false;
  } else {
    return false;
  }
  return true;
}
template <typename T>
bool PID_PROTON(T track)
{
   tpccut = 2.2, tofcut = 2;
  if (track.pt() < 1.1) {
    if (track.pt() < 0.85)
      tpccut = 2.2;
    else
      tpccut = 1;
    if (fabs(track.tpcNSigmaPr()) > tpccut)
      return false;
  } else if (track.hasTOF()) {
    if ((fabs(track.tpcNSigmaPr()) > tpccut) || (fabs(track.tofNSigmaPr()) > tofcut))
      return false;
  } else {
    return false;
  }
  return true;
}
//-----------------------------------------------------------------------------
template <typename T>
bool (*pidarray[4])(T) = {NO_PID, PID_PION, PID_KAON, PID_PROTON}; // Array of PID functions

namespace o2::aod
{
namespace idr2p2columns
{
DECLARE_SOA_COLUMN(BinNPIDFlag, binNpid, int8_t); // Flag tracks without proper binning as -1, and indicate type of particle 0->un-Id, 1->pion, 2->kaon, 3->proton
} // namespace idr2p2columns
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", idr2p2columns::BinNPIDFlag);
} // namespace o2::aod
struct FillFlagsTable {
  Produces<aod::Flags> ftable;
  Configurable<std::vector<std::vector<float>>>TPCnsigmacuts{"TPCnsigmacuts",{{2.5},{2,1,0.6,2},{2.2,1,2.2}},"TPC Nsigma cuts for different particle species"};
  Configurable<std::vector<std::vector<float>>>TPCpTranges{"TPCpTranges",{{2},{0.45,0,55,0.6,2},{0.85,1.1,2}},"TPC pT ranges for different particle species"};
  Configurable<std::vector<std::vector<float>>>TOFnsigmacuts{"TOFnsigmacuts",{{2.5},{2},{2}},"TOF Nsigma cuts for different particle species"};
  Configurable<std::vector<std::vector<float>>>TOFpTranges{"TOFpTranges",{{0.6},{0.6},{1.1}},"TOF pT ranges for different particle species"};

  Configurable<float> TPCnsigMax{"TPCnsigMax",5,"Maximum nsigma for identification (TPC)"};
  Configurable<float> TPCnsigSEP{"TPCnsigSEP",5,"Minimum separation from other particles (TPC)"};
  Configurable<float> TOFnsigMax{"TOFnsigMax",5,"Maximum nsigma for identification (TOF)"};
  Configurable<float> TOFnsigSEP{"TOFnsigSEP",5,"Minimum separation from other particles (TOF)"};
  Configurable<bool> ifTOFreject{"ifTOFreject",false,"Whether track should be rejected when TOF can't identify"};

  template <typename T>
  int PID_trial(T track)
  {
    float tpcnsigmaarray[3]={fabs(track.tpcNSigmaPi()),fabs(track.tpcNSigmaKa()),fabs(track.tpcNSigmaPr())},
      tofnsigmaarray[3]={fabs(track.tofNSigmaPi()),fabs(track.tofNSigmaKa()),fabs(track.tofNSigmaPr())};

int p1,p2,pid;
float nsigp1,nsigp2;
nsigp1=6;nsigp2=6;
for(int i=0;i<3;i++)
{
  if(nsigp1>tpcnsigmaarray[i])
  {
    nsigp2=nsigp1;
    p2=p1;
    p1=i;
    nsigp1=tpcnsigmaarray[i];
  }
  else if(nsigp2>tpcnsigmaarray[i])
  {
    p2=i;
    nsigp2=tpcnsigmaarray[i];
  }
}

if(nsigp1>TPCnsigMax)
  pid=0;
else if(nsigp2>TPCnsigSEP)
{
  pid=p1+1;
}
else if(track.hasTOF()) //Use TOF
{
  pid=p1+1;
  nsigp1=6;nsigp2=6;
  for(int i=0;i<3;i++)
  {
    if(nsigp1>tofnsigmaarray[i])
    {
      nsigp2=nsigp1;
      p2=p1;
      p1=i;
      nsigp1=tofnsigmaarray[i];
    }
    else if(nsigp2>tofnsigmaarray[i])
    {
      p2=i;
      nsigp2=tofnsigmaarray[i];
    }
  }
  if((nsigp1<TOFnsigMax)&&(nsigp2>TOFnsigSEP))
  {
    pid=p1+1;
  }
  else if(ifTOFreject) pid=0; //If we decide to rejct the track when TOF cant identify
  else pid=p1+1;
}
else if(ifTOFreject) pid=0;
else pid=p1+1;

return pid;
  }
  bool PID(float trackpt,float tracknsigmatpc,float tracknsigmatof,int species)
  {
    int tpcindex=-1,tofindex=-1;
    auto tpcpt=(std::vector<std::vector<float>>)TPCpTranges;
    auto tpcnsigma=(std::vector<std::vector<float>>)TPCnsigmacuts;
    auto tofpt=(std::vector<std::vector<float>>)TOFpTranges;
    auto tofnsigma=(std::vector<std::vector<float>>)TOFnsigmacuts;
    for(int i=0;i<tpcpt[species].size();i++)
      if(trackpt<tpcpt[species][i])
        tpcindex=i;
    for(int i=0;i<tofpt[species].size();i++)
      if(trackpt<tofpt[species][i])
        tofindex=i;
    if((tracknsigmatpc<tpcnsigma[species][tpcindex])&&(tofindex==-1))
      return true;
    else if((tracknsigmatpc<tpcnsigma[species][tpcindex])&&(tracknsigmatof<tofnsigma[species][tofindex]))
        return true;
    else
      return false;
  }
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    histos.add("tst","test",kTH1F,{{5,-0.5,4.5,"count"}});

    histos.add("nsigmatpc","N#sigma_{TPC}",kTH1F,{{50,0,2.1,"p_T"},{25,-6,6,"N#sigma"}});
    histos.add("nsigmatof","N#sigma_{TOF}",kTH1F,{{50,0,2.1,"p_T"},{25,-6,6,"N#sigma"}});
  }

  void processData(soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::TracksExtra> const& tracks)
  {
    int8_t etabin, phibin, binNpid;
    for (auto track : tracks) {
      etabin = (track.eta() + 0.8) * 15; // 15= 24/1.6
      phibin = 36 * track.phi() / (2 * constants::math::PI);
      if ((etabin < 0) || (etabin >= 24) || (phibin < 0) || (phibin >= 36)) {
        binNpid = -1;
      } else {
        float nsigma_array[3]={track.tpcNSigmaPi(),track.tpcNSigmaKa(),track.tpcNSigmaPr()};
        binNpid = 0;
        for (int8_t i = 0; i < 4; i++) {
          if (pidarray<decltype(track)>[i](track))
            binNpid = binNpid * 10 + i;
          if (binNpid > 10) // If a track is identified as two different tracks.
          {
            if (fabs(nsigma_array[(binNpid / 10) - 1]) < fabs(nsigma_array[(binNpid % 10) - 1])) // The track is identified as the particle whose |nsigma| is the least.
              binNpid /= 10;
            else
              binNpid %= 10;
          }
        }
        float tofnsigma_array[3]={track.tofNSigmaPi(),track.tofNSigmaKa(),track.tofNSigmaPr()};
        for (int8_t i = 1; i < 4; i++)
          if((pidarray<decltype(track)>[i](track))!=PID(track.pt(),fabs(nsigma_array[i-1]),fabs(tofnsigma_array[i-1]),i-1))
            histos.fill(HIST("tst"),1);
        histos.fill(HIST("nsigmatpc"),nsigma_array[binNpid]);
        histos.fill(HIST("nsigmatof"),tofnsigma_array[binNpid]);
      }
      ftable(binNpid);
    }
  }
  PROCESS_SWITCH(FillFlagsTable,processData,"Process Data",true);

  void processMC(soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::TracksExtra,aod::McTrackLabels> const& recotracks, aod::McParticles const& gentracks)
  {
  }
  PROCESS_SWITCH(FillFlagsTable,processMC, "Process MC Data",false);
};

struct r2p24id {

  Configurable<float> minpT{"minpT", 0.2, "Minimum pT"};
  Configurable<float> maxpT{"maxpT", 2.0, "Maximum pT"};
  Configurable<float> trackpartition{"trackpartition", 1.0, "where(in pT) to partition"};

  Configurable<int8_t> pid_particle1{"pid_particle1", 1, "Define particle1 type"}; //1->Pion, 2->Kaon, 3->Proton
  Configurable<int8_t> pid_particle2{"pid_particle2", 1, "Define particle2 type"};

  Configurable<bool> iftrackpartition{"iftrackpartition", false, "If track partition is needed"};
  Configurable<bool> ifpid{"ifpid", false, "If PID is needed"};

  struct histarray {
    std::shared_ptr<TH2> h2d_1p[2][2][2], h2d_2p[2][2][4];
    std::shared_ptr<TH1> h1d_1p[2][2];
  } hist;

  int mult1, mult2;
  uint8_t etabin1, phibin1, etabin2, phibin2;
  int8_t sign1, sign2;
  bool iftrack2;

  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  SliceCache cache;

  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > minpT&& aod::track::pt < maxpT;
  Filter globalfilter = requireGlobalTrackInFilter();
  Filter properbinfilter = aod::idr2p2columns::binNpid != -1;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>>> Tracks_set1 = (ifnode(iftrackpartition.node() == true, (aod::track::pt < trackpartition), true) && ifnode(ifpid.node() == true, (aod::idr2p2columns::binNpid == pid_particle1), true));
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>>> Tracks_set2 = (ifnode(iftrackpartition.node() == true, (aod::track::pt > trackpartition), true) && ifnode(ifpid.node() == true, (aod::idr2p2columns::binNpid == pid_particle2), true));
  //---------------------------------------------------------------------------

  void init(InitContext const&)
  {
    iftrack2 = (((int8_t)pid_particle1 != (int8_t)pid_particle2)&&((bool)ifpid)) || ((bool)iftrackpartition); //denotes whether the partiton1 is different from partition2
    //-----Defining Histograms---------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution Particle", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution Particle", kTH1D, {eta});
    histos.add("h1d_n1_ptP1", "p_T for +ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM1", "p_T for -ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptP2", "p_T for +ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM2", "p_T for -ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP1", "#rho_1 for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM1", "#rho_1 for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiP2", "#rho_1 for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM2", "#rho_1 for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12", "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21", "#rho_2 for +ve_2 -ve_1 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP1", "p_T for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM1", "p_T for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiP2", "p_T for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM2", "p_T for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12", "p_Tp_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21", "p_Tp_T for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12", "p_Tn for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21", "p_Tn for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    //---------------------------------------------------------------------------
    //-----Histogram Arrays------------------------------------------------------
    hist.h2d_1p[0][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM1"));
    hist.h2d_1p[0][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM1"));
    hist.h2d_1p[0][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP1"));
    hist.h2d_1p[0][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP1"));
    hist.h2d_1p[1][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM2"));
    hist.h2d_1p[1][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM2"));
    hist.h2d_1p[1][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP2"));
    hist.h2d_1p[1][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP2"));

    hist.h2d_2p[0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"));
    hist.h2d_2p[0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"));
    if (iftrack2) {
      hist.h2d_2p[0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"));
      hist.h2d_2p[0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"));
    } else {
      hist.h2d_2p[0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"));
      hist.h2d_2p[0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"));
    }
    hist.h2d_2p[1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"));
    hist.h2d_2p[1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"));
    hist.h2d_2p[1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"));

    hist.h1d_1p[0][0] = histos.template get<TH1>(HIST("h1d_n1_ptP1"));
    hist.h1d_1p[0][1] = histos.template get<TH1>(HIST("h1d_n1_ptM1"));
    hist.h1d_1p[1][0] = histos.template get<TH1>(HIST("h1d_n1_ptP2"));
    hist.h1d_1p[1][1] = histos.template get<TH1>(HIST("h1d_n1_ptM2"));
    //---------------------------------------------------------------------------
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    auto tracks1 = Tracks_set1->sliceByCached(aod::track::collisionId, filteredCollision.globalIndex(), cache);
    auto tracks2 = Tracks_set2->sliceByCached(aod::track::collisionId, filteredCollision.globalIndex(), cache);
    mult1 = tracks1.size();
    mult2 = tracks2.size();
    if ((iftrack2 && ((mult1 < 1) || (mult2 < 1))) || ((!iftrack2) && (mult1 < 2))) // Reject Collisions without sufficient particles
      return;

    for (auto track1 : tracks) {
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
    }
    histos.fill(HIST("h1i_n1_multPM"), mult1 + mult2);

    for (auto track1 : tracks1) {
      //---Single Particle Distribution (particle1)----------------------------------------
      sign1 = (track1.sign() + 1) / 2;
      hist.h1d_1p[0][sign1]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt())); // h1d_n1_pt*1
      hist.h2d_1p[0][sign1][0]->Fill(track1.eta(), track1.phi());                                // h2d_n1_etaphi*1
      hist.h2d_1p[0][sign1][1]->Fill(track1.eta(), track1.phi(), track1.pt());                   // h2d_pt_etaPhi*1
      //-----------------------------------------------------------------------
    }
    if (iftrack2) {
      for (auto track2 : tracks2) {
        //---Single Particle Distribution (particle2)----------------------------------------
        sign2 = (track2.sign() + 1) / 2;
        hist.h1d_1p[1][sign2]->Fill(track2.pt(), 1.0 / (2.0 * constants::math::PI * track2.pt())); // h1d_n1_pt*2
        hist.h2d_1p[1][sign2][0]->Fill(track2.eta(), track2.phi());                                // h2d_n1_etaphi*2
        hist.h2d_1p[1][sign2][1]->Fill(track2.eta(), track2.phi(), track2.pt());                   // h2d_pt_etaPhi*2
        //-----------------------------------------------------------------------
      }
    }
    for (auto track1 : tracks1) {
      sign1 = (track1.sign() + 1) / 2;
      etabin1 = (track1.eta() + 0.8) * 15; // 15= 24/1.6
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);

      for (auto track2 : tracks2) {

        if (track1.index() == track2.index())
          continue;

        etabin2 = (track2.eta() + 0.8) * 15; // 15= 24/1.6
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        sign2 = (track2.sign() + 1) / 2;

        //-----Two Particle Distribution---------------------------------------
        hist.h2d_2p[sign1][sign2][0]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);                            // h2d_n2_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][1]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());               // h2d_npt_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][2]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());               // h2d_ptn_eta1Phi1Eta2Phi2*
        hist.h2d_2p[sign1][sign2][3]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt()); // h2d_ptpt_eta1Phi1Eta2Phi2*
        //---------------------------------------------------------------------
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec{
    adaptAnalysisTask<FillFlagsTable>(cfgc)/*,
    adaptAnalysisTask<r2p24id>(cfgc),*/
    };
}
