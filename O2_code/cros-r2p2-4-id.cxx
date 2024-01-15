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
#include <string.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//-----PID Functions-----------------------------------------------------------
template <typename T>
bool NO_PID(T track)
{
  return true;
}
template <typename T>
bool PID_PION(T track)
{
  float tpccut = 2.5, tofcut = 2.5;
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
  float tpccut = 2, tofcut = 2;
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
  float tpccut = 2.2, tofcut = 2;
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

/*template <typename H,typename T1, typename T2, typename T3>
void CommonInit1(H& histos, T1& h2d_1p, T2& h2d_2p, T3& h1d_1p) // Init template for like particle correlation
{
  //-----Defining Histograms---------------------------------------------------
  const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
  histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
  histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
  histos.add("h1d_n1_ptP", "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
  histos.add("h1d_n1_ptM", "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
  histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
  histos.add("h2d_n1_etaPhiP", "#rho_1 for +ve particles", kTH2D, {eta, phi});
  histos.add("h2d_n1_etaPhiM", "#rho_1 for -ve particles", kTH2D, {eta, phi});
  histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_n2_eta1Phi1Eta2Phi2PM", "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_pt_etaPhiP", "p_T for +ve particles", kTH2D, {eta, phi});
  histos.add("h2d_pt_etaPhiM", "p_T for -ve particles", kTH2D, {eta, phi});
  histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM", "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM", "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_npt_eta1Phi1Eta2Phi2PM", "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
  histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
  histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
  histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
  histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
  histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
  //---------------------------------------------------------------------------
  //-----Histogram Arrays------------------------------------------------------
  h2d_1p[0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM"));
  h2d_1p[0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM"));
  h2d_1p[1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP"));
  h2d_1p[1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP"));
  h2d_2p[0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"));
  h2d_2p[1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"));
  h2d_2p[1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"));
  h2d_2p[1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"));
  h2d_2p[1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"));
  h2d_2p[2][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"));
  h2d_2p[2][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"));
  h2d_2p[2][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"));
  h2d_2p[2][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"));
  h1d_1p[0] = histos.template get<TH1>(HIST("h1d_n1_ptP"));
  h1d_1p[1] = histos.template get<TH1>(HIST("h1d_n1_ptM"));
  //---------------------------------------------------------------------------
}

template <typename H,typename T1, typename T2, typename T3>
void CommonInit2(H& histos,  T1& h2d_1p, T2& h2d_2p, T3& h1d_1p) // Init template for cross-correlation
{
  //-----Defining Histograms---------------------------------------------------
  const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
  histos.add("h1d_n1_phi1", "#phi distribution Particle 1", kTH1D, {phi});
  histos.add("h1d_n1_eta1", "#eta distribution Particle 1", kTH1D, {eta});
  histos.add("h1d_n1_phi2", "#phi distribution Particle 2", kTH1D, {phi});
  histos.add("h1d_n1_eta2", "#eta distribution Particle 2", kTH1D, {eta});
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
  histos.add("h2d_pt_etaPhiP1", "p_T for +ve_1", kTH2D, {eta, phi});
  histos.add("h2d_pt_etaPhiM1", "p_T for -ve_1", kTH2D, {eta, phi});
  histos.add("h2d_pt_etaPhiP2", "p_T for +ve_2", kTH2D, {eta, phi});
  histos.add("h2d_pt_etaPhiM2", "p_T for -ve_2", kTH2D, {eta, phi});
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
  histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
  histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
  histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
  histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
  //---------------------------------------------------------------------------
  //-----Histogram Arrays------------------------------------------------------
  h2d_1p[0][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM1"));
  h2d_1p[0][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM1"));
  h2d_1p[0][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP1"));
  h2d_1p[0][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP1"));
  h2d_1p[1][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM2"));
  h2d_1p[1][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM2"));
  h2d_1p[1][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP2"));
  h2d_1p[1][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP2"));
  h2d_2p[0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"));
  h2d_2p[0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"));
  h2d_2p[0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"));
  h2d_2p[0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"));
  h2d_2p[0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"));
  h2d_2p[1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"));
  h2d_2p[1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"));
  h2d_2p[1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"));
  h2d_2p[1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"));
  h2d_2p[1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"));
  h2d_2p[1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"));
  h2d_2p[1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"));
  h2d_2p[1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"));
  h1d_1p[0][0] = histos.template get<TH1>(HIST("h1d_n1_ptP1"));
  h1d_1p[0][1] = histos.template get<TH1>(HIST("h1d_n1_ptM1"));
  h1d_1p[1][0] = histos.template get<TH1>(HIST("h1d_n1_ptP2"));
  h1d_1p[1][1] = histos.template get<TH1>(HIST("h1d_n1_ptM2"));
  //---------------------------------------------------------------------------
}

template <typename TrackTable, typename Histograms,typename T1, typename T2, typename T3>
void CommonProcess1(TrackTable tracks, Histograms& histos,  T1& h2d_1p, T2& h2d_2p, T3& h1d_1p, int p) // Process template for like particle correlation
{
  int mult = 0, sign;
  int etabin1, etabin2, phibin1, phibin2;
  for (auto track1 : tracks) {
    //-----Single Particle Distribution----------------------------------------
    //-----PION PID------------------------------------------------------------
    histos.fill(HIST("h1d_n1_pt"), track1.pt());
    histos.fill(HIST("h2d_n1_tof1"), track1.p(), track1.beta());
    histos.fill(HIST("h2d_n1_tpc1"), track1.p(), track1.tpcSignal());
    if (!pidarray<decltype(track1)>[p](track1))
      continue;
    histos.fill(HIST("h2d_n1_tof2"), track1.p(), track1.beta());
    histos.fill(HIST("h2d_n1_tpc2"), track1.p(), track1.tpcSignal());
    //-------------------------------------------------------------------------
    histos.fill(HIST("h1d_n1_phi"), track1.phi());
    histos.fill(HIST("h1d_n1_eta"), track1.eta());
    sign = (track1.sign() + 1) / 2;
    h1d_1p[sign]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
    h2d_1p[sign][0]->Fill(track1.eta(), track1.phi());
    h2d_1p[sign][1]->Fill(track1.eta(), track1.phi(), track1.pt());
    etabin1 = (track1.eta() + 0.8) * 15; //15= 24/1.6
    phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
    if ((etabin1 < 0) || (etabin1 >= 24) || (phibin1 < 0) || (phibin1 >= 36))
      continue;
    mult++;
    //------------------------------------------------------------------------
    for (auto track2 : tracks) {
      //-----Two Particle Distribution------------------------------------------
      if (track1.index() == track2.index())
        continue; // No auto-correlations
      //-----PION PID----------------------------------------------------------
      if (!pidarray<decltype(track2)>[p](track2))
        continue;
      //-----------------------------------------------------------------------
      etabin2 = (track2.eta() + 0.8) * 15; //15= 24/1.6
      phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
      if ((etabin2 < 0) || (etabin2 >= 24) || (phibin2 < 0) || (phibin2 >= 36))
        continue;
      sign = (track1.sign() + track2.sign() + 2) / 2;
      h2d_2p[sign][0]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
      h2d_2p[sign][1]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
      h2d_2p[sign][2]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
      h2d_2p[sign][3]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
      //-----------------------------------------------------------------------
    }
  }
  histos.fill(HIST("h1i_n1_multPM"), mult);
}

template <typename TrackTable, typename Histograms,typename T1, typename T2, typename T3>
void CommonProcess2(TrackTable tracks, Histograms& histos,  T1& h2d_1p, T2& h2d_2p, T3& h1d_1p, int name, int p2) // Process template for cross-correlation
{
  bool flag;
  int mult = 0, sign1, sign2;
  int etabin1, etabin2, phibin1, phibin2;
  for (auto track1 : tracks) {
    //-----Single Particle Distribution----------------------------------------
    //-----KAON PID Particle2--------------------------------------------------
    flag = pidarray<decltype(track1)>[p2](track1);
    if (flag == true) {
      histos.fill(HIST("h2d_n1_tof2"), track1.pt(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc2"), track1.pt(), track1.tpcSignal());
      histos.fill(HIST("h1d_n1_phi2"), track1.phi());
      histos.fill(HIST("h1d_n1_eta2"), track1.eta());
      sign1 = (track1.sign() + 1) / 2;
      h1d_1p[1][sign1]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      h2d_1p[1][sign1][0]->Fill(track1.eta(), track1.phi());
      h2d_1p[1][sign1][1]->Fill(track1.eta(), track1.phi(), track1.pt());
    }
    //-------------------------------------------------------------------------
    //-----PROTON PID Particle1------------------------------------------------
    if (!pidarray<decltype(track1)>[name](track1))
      continue;
    histos.fill(HIST("h2d_n1_tof1"), track1.pt(), track1.beta());
    histos.fill(HIST("h2d_n1_tpc1"), track1.pt(), track1.tpcSignal());
    histos.fill(HIST("h1d_n1_phi1"), track1.phi());
    histos.fill(HIST("h1d_n1_eta1"), track1.eta());
    sign1 = (track1.sign() + 1) / 2;
    h1d_1p[0][sign1]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
    h2d_1p[0][sign1][0]->Fill(track1.eta(), track1.phi());
    h2d_1p[0][sign1][1]->Fill(track1.eta(), track1.phi(), track1.pt());
    //-------------------------------------------------------------------------
    etabin1 = (track1.eta() + 0.8) * 15; //15= 24/1.6
    phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
    if ((etabin1 < 0) || (etabin1 >= 24) || (phibin1 < 0) || (phibin1 >= 36))
      continue;
    mult++;
    //-------------------------------------------------------------------------
    for (auto track2 : tracks) {
      //-----Two Particle Distribution-----------------------------------------
      if (track1.index() == track2.index())
        continue; // No auto-correlations
      //-----KAON PID Particle 2-----------------------------------------------
      if (!pidarray<decltype(track2)>[p2](track2))
        continue;
      //-----------------------------------------------------------------------
      etabin2 = (track2.eta() + 0.8) * 15; //15= 24/1.6
      phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
      if ((etabin2 < 0) || (etabin2 >= 24) || (phibin2 < 0) || (phibin2 >= 36))
        continue;
      sign2 = (track2.sign() + 1) / 2;
      h2d_2p[sign1][sign2][0]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
      h2d_2p[sign1][sign2][1]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
      h2d_2p[sign1][sign2][2]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
      h2d_2p[sign1][sign2][3]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
      //-----------------------------------------------------------------------
    }
  }
  histos.fill(HIST("h1i_n1_multPM"), mult);
}

// All particle correlation
struct r2p24ch {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2], h2d_2p[3][4];
  std::shared_ptr<TH1> h1d_1p[2];

  void init(InitContext const&)
  {
    CommonInit1(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess1(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 0);
  }
};

// Pion pair correlation
struct r2p24pi {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2], h2d_2p[3][4];
  std::shared_ptr<TH1> h1d_1p[2];

  void init(InitContext const&)
  {
    CommonInit1(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess1(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 1);
  }
};

// Kaon pair correlation
struct r2p24ka {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2], h2d_2p[3][4];
  std::shared_ptr<TH1> h1d_1p[2];

  void init(InitContext const&)
  {
    CommonInit1(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess1(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 2);
  }
};

// Proton pair correlation
struct r2p24pr {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2], h2d_2p[3][4];
  std::shared_ptr<TH1> h1d_1p[2];

  void init(InitContext const&)
  {
    CommonInit1(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess1(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 3);
  }
};

// Pion-Kaon pair correlation
struct r2p24pik {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2][2], h2d_2p[2][2][4];
  std::shared_ptr<TH1> h1d_1p[2][2];

  void init(InitContext const&)
  {
    CommonInit2(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess2(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 1, 2);
  }
};

// Pion-Proton pair correlation
struct r2p24pip {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2][2], h2d_2p[2][2][4];
  std::shared_ptr<TH1> h1d_1p[2][2];

  void init(InitContext const&)
  {
    CommonInit2(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess2(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 1, 3);
  }
};

// Proton-Kaon pair correlation
struct r2p24pk {
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TH2> h2d_1p[2][2][2], h2d_2p[2][2][4];
  std::shared_ptr<TH1> h1d_1p[2][2];

  void init(InitContext const&)
  {
    CommonInit2(histos, h2d_1p, h2d_2p, h1d_1p);
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    CommonProcess2(tracks, histos, h2d_1p, h2d_2p, h1d_1p, 3, 2);
  }
};
*/
namespace o2::aod
{
namespace columns
{
DECLARE_SOA_COLUMN(BinFlag, binflag, bool); // To flag tracks without proper binning
}
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", columns::BinFlag);
} // namespace o2::aod
struct FillFlagsTable {
  Produces<aod::Flags> ftable;
  void process(aod::Tracks const& tracks)
  {
    int etabin, phibin;
    for (auto track : tracks) {
      etabin = (track.eta() + 0.8) * 15; // 15= 24/1.6
      phibin = 36 * track.phi() / (2 * constants::math::PI);
      if ((etabin < 0) || (etabin >= 24) || (phibin < 0) || (phibin >= 36))
        ftable(false);
      else
        ftable(true);
    }
  }
};
struct r2p24id {
  Configurable<float> minpT{"minpT", 0.2, "Minimum pT"};
  Configurable<float> maxpT{"maxpT", 2.0, "Maximum pT"};
  //-----Track&Event Selection-------------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > minpT&& aod::track::pt < maxpT;
  Filter globalfilter = requireGlobalTrackInFilter();
  Filter properbinfilter = aod::columns::binflag == true; //'aod::Flags::binflag == true;' did not work
  //---------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct histarray {
    std::shared_ptr<TH2> h2d_1p[4][2][2], h2d_2p[7][2][2][4];
    std::shared_ptr<TH1> h1d_1p[4][2];
  } hist;

  void init(InitContext const&)
  {
    //-----Defining Histograms-------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    char name[4][3] = {"ch", "pi", "ka", "pr"}, histname[50];
    histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    for (int i = 0; i < 7; i++) {
      if (i < 4) { // Single particle distribution histograms
        snprintf(histname, sizeof(histname), "h1d_n1_ptP_%s", name[i]);
        histos.add(histname, "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
        snprintf(histname, sizeof(histname), "h1d_n1_ptM_%s", name[i]);
        histos.add(histname, "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
        snprintf(histname, sizeof(histname), "h2d_n1_etaPhiP_%s", name[i]);
        histos.add(histname, "#rho_1 for +ve particles", kTH2D, {eta, phi});
        snprintf(histname, sizeof(histname), "h2d_n1_etaPhiM_%s", name[i]);
        histos.add(histname, "#rho_1 for -ve particles", kTH2D, {eta, phi});
        snprintf(histname, sizeof(histname), "h2d_pt_etaPhiP_%s", name[i]);
        histos.add(histname, "p_T for +ve particles", kTH2D, {eta, phi});
        snprintf(histname, sizeof(histname), "h2d_pt_etaPhiM_%s", name[i]);
        histos.add(histname, "p_T for -ve particles", kTH2D, {eta, phi});
      }
      //---Two patricle distribution histograms--------------------------------
      int e = i, f = i;
      if (i > 3) {
        e = 1;
        f = 2;
        if (i == 5)
          f = 3;
        else if (i == 6)
          e = 3;
        // There are two kinds of  "********PM_*" histograms here, one is for when track1 is +ve the other is for when track2 is +ve.
        snprintf(histname, sizeof(histname), "h2d_n2_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptpt_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "p_Tp_T for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptn_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "p_Tn for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_npt_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "np_T for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_n2_eta1Phi1Eta2Phi2PM_%s%s", name[f], name[e]);
        histos.add(histname, "#rho_2 for -ve_1 +ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptpt_eta1Phi1Eta2Phi2PM_%s%s", name[f], name[e]);
        histos.add(histname, "p_Tp_T for -ve_1 +ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptn_eta1Phi1Eta2Phi2PM_%s%s", name[f], name[e]);
        histos.add(histname, "p_Tn for -ve_1 +ve_2 particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_npt_eta1Phi1Eta2Phi2PM_%s%s", name[f], name[e]);
        histos.add(histname, "np_T for -ve_1 +ve_2 particles", kTH2D, {etaphi1, etaphi2});
      } else {
        snprintf(histname, sizeof(histname), "h2d_n2_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptpt_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_ptn_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
        snprintf(histname, sizeof(histname), "h2d_npt_eta1Phi1Eta2Phi2PM_%s%s", name[e], name[f]);
        histos.add(histname, "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
      }
      snprintf(histname, sizeof(histname), "h2d_n2_eta1Phi1Eta2Phi2PP_%s%s", name[e], name[f]);
      histos.add(histname, "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_n2_eta1Phi1Eta2Phi2MM_%s%s", name[e], name[f]);
      histos.add(histname, "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_ptpt_eta1Phi1Eta2Phi2PP_%s%s", name[e], name[f]);
      histos.add(histname, "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_ptpt_eta1Phi1Eta2Phi2MM_%s%s", name[e], name[f]);
      histos.add(histname, "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_ptn_eta1Phi1Eta2Phi2PP_%s%s", name[e], name[f]);
      histos.add(histname, "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_ptn_eta1Phi1Eta2Phi2MM_%s%s", name[e], name[f]);
      histos.add(histname, "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_npt_eta1Phi1Eta2Phi2PP_%s%s", name[e], name[f]);
      histos.add(histname, "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
      snprintf(histname, sizeof(histname), "h2d_npt_eta1Phi1Eta2Phi2MM_%s%s", name[e], name[f]);
      histos.add(histname, "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
      //-----------------------------------------------------------------------
    }
    //-------------------------------------------------------------------------
    //-----Assigning histogram pointers to an Array----------------------------
    //------Single Particle..........------------------------------------------
    hist.h1d_1p[0][0] = histos.template get<TH1>(HIST("h1d_n1_ptM_ch"));
    hist.h1d_1p[0][1] = histos.template get<TH1>(HIST("h1d_n1_ptP_ch"));
    hist.h2d_1p[0][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM_ch"));
    hist.h2d_1p[0][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM_ch"));
    hist.h2d_1p[0][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP_ch"));
    hist.h2d_1p[0][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP_ch"));
    hist.h1d_1p[1][0] = histos.template get<TH1>(HIST("h1d_n1_ptM_pi"));
    hist.h1d_1p[1][1] = histos.template get<TH1>(HIST("h1d_n1_ptP_pi"));
    hist.h2d_1p[1][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM_pi"));
    hist.h2d_1p[1][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM_pi"));
    hist.h2d_1p[1][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP_pi"));
    hist.h2d_1p[1][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP_pi"));
    hist.h1d_1p[2][0] = histos.template get<TH1>(HIST("h1d_n1_ptM_ka"));
    hist.h1d_1p[2][1] = histos.template get<TH1>(HIST("h1d_n1_ptP_ka"));
    hist.h2d_1p[2][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM_ka"));
    hist.h2d_1p[2][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM_ka"));
    hist.h2d_1p[2][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP_ka"));
    hist.h2d_1p[2][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP_ka"));
    hist.h1d_1p[3][0] = histos.template get<TH1>(HIST("h1d_n1_ptM_pr"));
    hist.h1d_1p[3][1] = histos.template get<TH1>(HIST("h1d_n1_ptP_pr"));
    hist.h2d_1p[3][0][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiM_pr"));
    hist.h2d_1p[3][0][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiM_pr"));
    hist.h2d_1p[3][1][0] = histos.template get<TH2>(HIST("h2d_n1_etaPhiP_pr"));
    hist.h2d_1p[3][1][1] = histos.template get<TH2>(HIST("h2d_pt_etaPhiP_pr"));
    //-------------------------------------------------------------------------
    //------Two Particle............-------------------------------------------
    hist.h2d_2p[0][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_chch"));
    hist.h2d_2p[0][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_chch"));
    hist.h2d_2p[0][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_chch"));
    hist.h2d_2p[0][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_chch"));
    hist.h2d_2p[0][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_chch"));
    hist.h2d_2p[0][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_chch"));
    hist.h2d_2p[0][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_chch"));
    hist.h2d_2p[0][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_chch"));
    hist.h2d_2p[0][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_chch"));
    hist.h2d_2p[1][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_pipi"));
    hist.h2d_2p[1][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_pipi"));
    hist.h2d_2p[1][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_pipi"));
    hist.h2d_2p[1][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_pipi"));
    hist.h2d_2p[1][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_pipi"));
    hist.h2d_2p[1][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_pipi"));
    hist.h2d_2p[1][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_pipi"));
    hist.h2d_2p[1][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_pipi"));
    hist.h2d_2p[1][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_pipi"));
    hist.h2d_2p[2][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_kaka"));
    hist.h2d_2p[2][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_kaka"));
    hist.h2d_2p[2][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_kaka"));
    hist.h2d_2p[2][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_kaka"));
    hist.h2d_2p[2][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_kaka"));
    hist.h2d_2p[2][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_kaka"));
    hist.h2d_2p[2][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_kaka"));
    hist.h2d_2p[2][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_kaka"));
    hist.h2d_2p[2][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_kaka"));
    hist.h2d_2p[3][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_prpr"));
    hist.h2d_2p[3][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_prpr"));
    hist.h2d_2p[3][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_prpr"));
    hist.h2d_2p[3][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_prpr"));
    hist.h2d_2p[3][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_prpr"));
    hist.h2d_2p[3][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_prpr"));
    hist.h2d_2p[3][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_prpr"));
    hist.h2d_2p[3][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_prpr"));
    hist.h2d_2p[3][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_prpr"));
    hist.h2d_2p[4][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_pika"));
    hist.h2d_2p[4][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_pika"));
    hist.h2d_2p[4][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_pika"));
    hist.h2d_2p[4][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_pika"));
    hist.h2d_2p[4][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_kapi"));
    hist.h2d_2p[4][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_kapi"));
    hist.h2d_2p[4][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_kapi"));
    hist.h2d_2p[4][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_kapi"));
    hist.h2d_2p[4][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_pika"));
    hist.h2d_2p[4][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_pika"));
    hist.h2d_2p[4][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_pika"));
    hist.h2d_2p[4][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_pika"));
    hist.h2d_2p[4][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_pika"));
    hist.h2d_2p[4][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_pika"));
    hist.h2d_2p[4][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_pika"));
    hist.h2d_2p[4][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_pika"));
    hist.h2d_2p[5][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_pipr"));
    hist.h2d_2p[5][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_pipr"));
    hist.h2d_2p[5][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_pipr"));
    hist.h2d_2p[5][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_pipr"));
    hist.h2d_2p[5][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_prpi"));
    hist.h2d_2p[5][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_prpi"));
    hist.h2d_2p[5][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_prpi"));
    hist.h2d_2p[5][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_prpi"));
    hist.h2d_2p[5][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_pipr"));
    hist.h2d_2p[5][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_pipr"));
    hist.h2d_2p[5][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_pipr"));
    hist.h2d_2p[5][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_pipr"));
    hist.h2d_2p[5][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_pipr"));
    hist.h2d_2p[5][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_pipr"));
    hist.h2d_2p[5][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_pipr"));
    hist.h2d_2p[5][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_pipr"));
    hist.h2d_2p[6][0][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2MM_prka"));
    hist.h2d_2p[6][0][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2MM_prka"));
    hist.h2d_2p[6][0][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM_prka"));
    hist.h2d_2p[6][0][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM_prka"));
    hist.h2d_2p[6][0][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_kapr"));
    hist.h2d_2p[6][0][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_kapr"));
    hist.h2d_2p[6][0][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_kapr"));
    hist.h2d_2p[6][0][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_kapr"));
    hist.h2d_2p[6][1][0][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PM_prka"));
    hist.h2d_2p[6][1][0][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PM_prka"));
    hist.h2d_2p[6][1][0][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM_prka"));
    hist.h2d_2p[6][1][0][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM_prka"));
    hist.h2d_2p[6][1][1][0] = histos.template get<TH2>(HIST("h2d_n2_eta1Phi1Eta2Phi2PP_prka"));
    hist.h2d_2p[6][1][1][1] = histos.template get<TH2>(HIST("h2d_npt_eta1Phi1Eta2Phi2PP_prka"));
    hist.h2d_2p[6][1][1][2] = histos.template get<TH2>(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP_prka"));
    hist.h2d_2p[6][1][1][3] = histos.template get<TH2>(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP_prka"));
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    int sign1, sign2, mult = 0;
    int etabin1, phibin1, etabin2, phibin2;

    for (auto track1 : tracks) {
      mult++;
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
      //---Single Particle Distribution----------------------------------------
      for (int i = 0; i < 4; i++) // Filling histograms for different particle species
      {
        if (pidarray<decltype(track1)>[i](track1)) // PID function
        {
          sign1 = (track1.sign() + 1) / 2;
          hist.h1d_1p[i][sign1]->Fill(track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt())); // h1d_n1_pt*
          hist.h2d_1p[i][sign1][0]->Fill(track1.eta(), track1.phi());                                // h2d_n1_etaphi*
          hist.h2d_1p[i][sign1][1]->Fill(track1.eta(), track1.phi(), track1.pt());                   // h2d_pt_etaPhi*
        }
      }
      //-----------------------------------------------------------------------
      etabin1 = (track1.eta() + 0.8) * 15; // 15= 24/1.6
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);

      for (auto track2 : tracks) {
        if (track1.index() == track2.index())
          continue;
        etabin2 = (track2.eta() + 0.8) * 15; // 15= 24/1.6
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        sign2 = (track2.sign() + 1) / 2;
        //-----Two Particle Distribution---------------------------------------
        for (int i = 0; i < 7; i++) // iterating through different combination of particles
        {
          int e = i, f = i;
          if (i > 3) { // combinations for cross correlation
            e = 1;
            f = 2;
            if (i == 5)
              f = 3;
            else if (i == 6)
              e = 3;
          }
          if ((pidarray<decltype(track1)>[e](track1)) && (pidarray<decltype(track2)>[f](track2))) {                                       // This If statement can be avoided by using a better PID Technique, where a PID flag can be assigned to the track and that flag can be used as an index for the histogram array. I am working on a better PID, but for now I have to stick with this.
            hist.h2d_2p[i][sign1][sign2][0]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);                            // h2d_n2_eta1Phi1Eta2Phi2*
            hist.h2d_2p[i][sign1][sign2][1]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());               // h2d_npt_eta1Phi1Eta2Phi2*
            hist.h2d_2p[i][sign1][sign2][2]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());               // h2d_ptn_eta1Phi1Eta2Phi2*
            hist.h2d_2p[i][sign1][sign2][3]->Fill(36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt()); // h2d_ptpt_eta1Phi1Eta2Phi2*
          }
        }
        //---------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    /*adaptAnalysisTask<r2p24ch>(cfgc),
    adaptAnalysisTask<r2p24pi>(cfgc),
    adaptAnalysisTask<r2p24ka>(cfgc),
    adaptAnalysisTask<r2p24pr>(cfgc),
    adaptAnalysisTask<r2p24pik>(cfgc),
    adaptAnalysisTask<r2p24pip>(cfgc),
    adaptAnalysisTask<r2p24pk>(cfgc),*/
    adaptAnalysisTask<FillFlagsTable>(cfgc),
    adaptAnalysisTask<r2p24id>(cfgc),
  };
}
