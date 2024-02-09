//#include "/Users/austinbaty/Desktop/FastJet/fastjet-install/include/fastjet/ClusterSequence.hh"
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <vector>

using namespace Pythia8;
//using namespace fastjet;
#include <iostream>

int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  Event& event = pythia.event;


  pythia.readString("Beams:frameType = 2");
  //electron
  pythia.readString("Beams:idA = 11");
  //muon
  //pythia.readString("Beams:idA = 13");
  //proton
  pythia.readString("Beams:idB = 2212");

  //specify the beam energies
  //electron: 5-30 GeV for EIC
  //proton:   50-250 GeV for EIC
  //muon ion collider: 960 GeV muon, 275 GeV proton
  pythia.settings.parm("Beams:eA = 30");
  pythia.settings.parm("Beams:eB = 250");
  
  //neutral current
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  //charged current
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  
  //minimum Q2 cut
  pythia.settings.parm("PhaseSpace:Q2Min", 1.0);

  //needed setting for DIS + showering
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  pythia.init();

  
  float genQScale = 0;
  float genJetPt = 0;
  float genJetEta = 0;
  float genJetPhi = 0;
  int   genJetChargedMultiplicity = 0;
  std::vector< std::vector<int> > gendau_chg;
  std::vector< std::vector<int> > gendau_pid;
  std::vector< std::vector<float> > gendau_pt;
  std::vector< std::vector<float> > gendau_eta;
  std::vector< std::vector<float> > gendau_phi;
  std::vector< std::vector<int> > gendau_mom;

  TFile * f = TFile::Open("EIC.root","recreate");
  TTree * trackTree = new TTree("trackTree","v1");
 
  trackTree->Branch("genQScale",&genQScale);
  trackTree->Branch("genJetEta",&genJetEta);
  trackTree->Branch("genJetPt",&genJetPt);
  trackTree->Branch("genJetPhi",&genJetPhi);
  trackTree->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);

  trackTree->Branch("genDau_chg",		&gendau_chg); 
  trackTree->Branch("genDau_pid",		&gendau_pid);	 
  trackTree->Branch("genDau_pt",		&gendau_pt);
  trackTree->Branch("genDau_eta",		&gendau_eta);	 
  trackTree->Branch("genDau_phi",		&gendau_phi );
  trackTree->Branch("genDau_mom",		&gendau_mom );


  // Begin event loop.
  int nEvent = 100;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn    = event[4].p();
    Vec4 peOut   = event[6].p();
    Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);


   
    
//    vector<PseudoJet> particles;
    int multiplicity = 0;
    //start loop at i=7 because we want to exclude the scattered lepton (which is i=6 index)
    for (int i = 7; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()){
        //basic kinematic cuts
        if( TMath::Abs( pythia.event[i].eta() ) > 4 || pythia.event[i].pT() < 0.2) continue;
        multiplicity++;
//        particles.push_back( PseudoJet(   pythia.event[i].px(),  pythia.event[i].py(),  pythia.event[i].pz(), pythia.event[i].e() ) ); 
      }

      //fill and clear branches
      trackTree->Fill();
      genQScale = 0;
      genJetPt = 0;
      genJetEta = 0;
      genJetPhi = 0;
      genJetChargedMultiplicity = 0;
      gendau_chg.clear();
      gendau_pid.clear();
      gendau_pt.clear();
      gendau_eta.clear();
      gendau_phi.clear();
      gendau_mom.clear();
    }
    //if(multiplicity<50) continue;
  
/*
    //jet clustering
    // choose a jet definition
    double R = 0.8;
    JetDefinition jet_def(antikt_algorithm, R);

    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
  
    


    // print out some infos
    cout << "Multiplicity: " << multiplicity << " Clustering with " << jet_def.description() << endl;

    // print the jets
    cout <<   "        pt y phi" << endl;
    for (unsigned i = 0; i < jets.size(); i++) {
      cout << "jet " << i << ": "<< jets[i].pt() << " " 
                     << jets[i].rap() << " " << jets[i].phi() << endl;
      vector<PseudoJet> constituents = jets[i].constituents();
      for (unsigned j = 0; j < constituents.size(); j++) {
        cout << "    constituent " << j << "'s pt: " << constituents[j].pt()
             << endl;
      }
    }
    */

  }


  return 0;
}
