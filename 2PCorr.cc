// main05.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; jet finding; slowjet; celljet;

// This is a simple test program.
// It studies jet production at the LHC, using SlowJet and CellJet.
// Note: the two finders are intended to construct approximately the same
// jet properties, but provides output in slightly different format,
// and have here not been optimized to show maximum possible agreement.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"

#include <vector>

#define PI 3.1415926

using namespace std;

double GetDeltaEta(double eta_trg, double eta_ass)
{
        double deltaEta=eta_ass-eta_trg;
        return deltaEta;
}

double GetDeltaPhi(double phi_trg, double phi_ass)
{     
        double deltaPhi=phi_ass-phi_trg;
        if(deltaPhi>PI)
          deltaPhi=deltaPhi-2*PI;
        if(deltaPhi<-PI/2)
          deltaPhi=deltaPhi+2*PI;
        return deltaPhi;
}

int main() {

  // Number of events, generated and listed ones.
  int nEvent    = 50000;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  Event& event = pythia.event;

  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 60.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Common parameters for the two jet finders.
  double etaMax   = 2.5;
  double radius   = 0.4;
  double pTjetMin = 60.;
  // Exclude neutrinos (and other invisible) from study.
  int    nSel     = 2;
  // Range and granularity of CellJet jet finder.
  int    nEta     = 80;
  int    nPhi     = 64;

  // Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  // Set up CellJet jet finder.
  CellJet cellJet( etaMax, nEta, nPhi, nSel);

  // Histograms. Note similarity in names, even when the two jet finders
  // do not calculate identically the same property (pT vs. ET, y vs. eta).
  TH1D* nJetsS = new TH1D("number of jets, SlowJet", "", 50, -0.5, 49.5);
  TH1D* nJetsC = new TH1D("number of jets, CellJet", "", 50, -0.5, 49.5);
  TH1D* nJetsD = new TH1D("number of jets, CellJet - SlowJet", "", 45, -22.5, 22.5);
  TH1D* eTjetsS = new TH1D("pT for jets, SlowJet", "", 100, 0., 500.);
  TH1D* eTjetsC = new TH1D("eT for jets, CellJet", "", 100, 0., 500.);
  TH1D* etaJetsS = new TH1D("y for jets, SlowJet", "", 100, -5., 5.);
  TH1D* etaJetsC = new TH1D("eta for jets, CellJet", "", 100, -5., 5.);
  TH1D* phiJetsS = new TH1D("phi for jets, SlowJwt", "", 100, -M_PI, M_PI);
  TH1D* phiJetsC = new TH1D("phi for jets, CellJet", "", 100, -M_PI, M_PI);
  TH1D* distJetsS = new TH1D("R distance between jets, SlowJet", "", 100, 0., 10.);
  TH1D* distJetsC = new TH1D("R distance between jets, CellJet", "", 100, 0., 10.);
  TH1D* eTdiffS = new TH1D("pT difference, SlowJet", "", 100, -100., 400.);
  TH1D* eTdiffC = new TH1D("eT difference, CellJet", "", 100, -100., 400.);

  TH1D* Nchhist = new TH1D("Nch","N_{ch}", 500, 0., 500.);
  TH2D* ptetahist = new TH2D("ptetahist",";#eta;p_{T} [GeV]", 50, -5., 5., 200, 0., 20.);

  double etamin=-2.4;
  double etamax=2.4;
  int NEtaBins=48;
  int NPhiBins=32;
  double etabinwidth=(etamax-etamin)/NEtaBins;
  double phibinwidth=PI/16.;

  TH2D* hSignal = new TH2D("signal",";#Delta#eta;#Delta#phi",
                                     NEtaBins+1,etamin-etamax-etabinwidth/2.,etamax-etamin+etabinwidth/2.,                    
	                             NPhiBins-1,-(PI-phibinwidth)/2.0,(PI*3.0-phibinwidth)/2.0);

  TH2D* hBackground = new TH2D("background",";#Delta#eta;#Delta#phi",
                                     NEtaBins+1,etamin-etamax-etabinwidth/2.,etamax-etamin+etabinwidth/2.,
                                     NPhiBins-1,-(PI-phibinwidth)/2.0,(PI*3.0-phibinwidth)/2.0);

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open("pythia_jets.root","recreate");

  vector<vector<TVector3>> pvectvect;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze Slowet jet properties. List first few.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    nJetsS->Fill( slowJet.sizeJet() );
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      eTjetsS->Fill( slowJet.pT(i) );
      etaJetsS->Fill( slowJet.y(i) );
      phiJetsS->Fill( slowJet.phi(i) );
    }

    // Fill SlowJet distance between jets.
    for (int i = 0; i < slowJet.sizeJet() - 1; ++i)
    for (int j = i +1; j < slowJet.sizeJet(); ++j) {
      double dEta = slowJet.y(i) - slowJet.y(j);
      double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsS->Fill( dR );
    }

    // Fill SlowJet pT-difference between jets (to check ordering of list).
    for (int i = 1; i < slowJet.sizeJet(); ++i)
      eTdiffS->Fill( slowJet.pT(i-1) - slowJet.pT(i) );

    // Analyze CellJet jet properties. List first few.
    cellJet. analyze( pythia.event, pTjetMin, radius );
    if (iEvent < nListJets) cellJet.list();

    // Fill CellJet inclusive jet distributions.
    nJetsC->Fill( cellJet.size() );
    for (int i = 0; i < cellJet.size(); ++i) {
      eTjetsC->Fill( cellJet.eT(i) );
      etaJetsC->Fill( cellJet.etaWeighted(i) );
      phiJetsC->Fill( cellJet.phiWeighted(i) );
    }

    // Fill CellJet distance between jets.
    for (int i = 0; i < cellJet.size() - 1; ++i)
    for (int j = i +1; j < cellJet.size(); ++j) {
      double dEta = cellJet.etaWeighted(i)
        - cellJet.etaWeighted(j);
      double dPhi = abs( cellJet.phiWeighted(i)
        - cellJet.phiWeighted(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsC->Fill( dR );
    }

    // Fill CellJet ET-difference between jets (to check ordering of list).
    for (int i = 1; i < cellJet.size(); ++i)
      eTdiffC->Fill( cellJet.eT(i-1) - cellJet.eT(i) );

    // Compare number of jets for the two finders.
    nJetsD->Fill( cellJet.size() - slowJet.sizeJet() );

    // Loop over particles 
    int nMult = 0;
    vector<TVector3> pvect;
    for (int i = 0; i < event.size(); ++i)
    {
      if (event[i].statusAbs() <=89 && event[i].statusAbs() >=81) {
        if(event[i].charge()==0) continue;
        ptetahist->Fill( event[i].eta(), event[i].pT() );
        if(fabs(event[i].eta())>2.5) continue; 
        if(event[i].pT()<0.5 || event[i].pT()>5) continue;
        nMult++;

        TVector3 pp;
        pp.SetPtEtaPhi(event[i].pT(),event[i].eta(),event[i].phi());
        pvect.push_back(pp);
      }
    }
    Nchhist->Fill(nMult);

    if(nMult<50) continue;

    for(int i=0; i<pvect.size(); i++)
    {
      TVector3 pp1 = pvect.at(i);

      // signal pairs
      for(int j=i+1; j<pvect.size(); j++)
      {
        TVector3 pp2 = pvect.at(j);

        double deltaPhi = pp1.DeltaPhi(pp2);
        double deltaEta = pp1.Eta()-pp2.Eta();

        hSignal->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
        hSignal->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
        hSignal->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
        hSignal->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
        hSignal->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
        hSignal->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
      }

      // background pairs
      for(int m=0; m<pvectvect.size(); m++)
      {
        vector<TVector3> pvect_bkg = pvectvect.at(m);
        for(int j=0; j<pvect_bkg.size(); j++)
        { 
          TVector3 pp2 = pvect_bkg.at(j);
        
          double deltaPhi = pp1.DeltaPhi(pp2);
          double deltaEta = pp1.Eta()-pp2.Eta();
        
          hBackground->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
          hBackground->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
          hBackground->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
          hBackground->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
          hBackground->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
          hBackground->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
        }
      }
    }

    if(pvectvect.size()==10) pvectvect.erase(pvectvect.begin());
    pvectvect.push_back(pvect);
/*
      for (int j = i+1; j < event.size(); ++j)
      {
        if (event[j].statusAbs() <=89 && event[j].statusAbs() >=81) {
          if(event[j].charge()==0) continue;
          if(fabs(event[i].eta())>2.5) continue; 
        
          double deltaPhi=GetDeltaPhi(event[i].phi(),event[j].phi());
          double deltaEta=GetDeltaEta(event[i].eta(),event[j].eta());

          hSignal->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
          hSignal->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0);
          hSignal->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
          hSignal->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0);
          hSignal->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
          hSignal->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0);
        }
      }
    }
*/

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << nJetsS << nJetsC << nJetsD << eTjetsS << eTjetsC
       << etaJetsS << etaJetsC << phiJetsS << phiJetsC
       << distJetsS << distJetsC << eTdiffS << eTdiffC;

  TH2D* hCorr = (TH2D*)hSignal->Clone("hCorr");
  hCorr->Divide(hBackground);

  //  Write tree.
  nJetsS->Write();
  nJetsC->Write();
  nJetsD->Write();
  eTjetsS->Write();
  eTjetsC->Write();
  etaJetsS->Write();
  etaJetsC->Write();
  phiJetsS->Write();
  phiJetsC->Write();
  distJetsS->Write();
  distJetsC->Write();
  eTdiffS->Write();
  eTdiffC->Write();

  Nchhist->Write();
  ptetahist->Write();

  hSignal->Write();
  hBackground->Write();
  hCorr->Write();

  delete file;


  // Done.
  return 0;
}
