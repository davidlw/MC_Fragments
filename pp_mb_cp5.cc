//external includes (root, fastjet, pythia)
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TInterpreter.h"

//standard cpp libraries
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

//namespaces
using namespace Pythia8;
using namespace fastjet;


int main(int argc, char *argv[]) {

//need to generate these dictionaries to output a vector of vectors...	
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

//get input job number (used for labeling outputs in batch submission if needed)
  if(argc != 2)
  {
    std::cout << "Usage: ./pp_highMultGen <job number>" << std::endl;
    return 1;
  }  
  int jobNumber = std::atoi(argv[1]);

  //cut on the minimum number of charged constituents in jets before the jet is output (0 is no cut)
  const float NchCut = 0;

//*********************************************** PYTHIA SETTINGS **************************************************
  Pythia pythia;
  Event& event = pythia.event;

  //for frame type 2 we specify both beams individually
  pythia.readString("Beams:frameType = 2");
  //proton for both beams
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  //specify the beam energies (13TeV pp here)
  pythia.readString("Beams:eA = 6500");
  pythia.readString("Beams:eB = 6500");

  //process  
  pythia.readString("HardQCD:all = on");
//  pythia.settings.parm("PhaseSpace:pTHatMin", 0.1);

  //(optional) general CMS settings
  pythia.readString("Check:epTolErr = 0.01");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");
  
  //(optional) CP5 tune settings if desired
  pythia.readString("Tune:pp 14");
  pythia.readString("Tune:ee 7");
  pythia.readString("MultipartonInteractions:bProfile=2");
  pythia.readString("MultipartonInteractions:pT0Ref=1.41");
  pythia.readString("MultipartonInteractions:coreRadius=0.7634");
  pythia.readString("MultipartonInteractions:coreFraction=0.63");
  pythia.readString("ColourReconnection:range=5.176");
  pythia.readString("SigmaTotal:zeroAXB=off");
  pythia.readString("SpaceShower:alphaSorder=2");
  pythia.readString("SpaceShower:alphaSvalue=0.118");
  pythia.readString("SigmaProcess:alphaSvalue=0.118");
  pythia.readString("SigmaProcess:alphaSorder=2");
  pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
  pythia.readString("MultipartonInteractions:alphaSorder=2");
  pythia.readString("TimeShower:alphaSorder=2");
  pythia.readString("TimeShower:alphaSvalue=0.118");
  pythia.readString("SigmaTotal:mode = 0");
  pythia.readString("SigmaTotal:sigmaEl = 21.89");
  pythia.readString("SigmaTotal:sigmaTot = 100.309");

  //Set the RNG seed to be based on the system clock (seed = 0 means use clock)
  //Should be set to another value if you want reproducable MC events
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
 
  //initialize
  pythia.init();

  
  //*********************************************** ROOT and Tree SETUP  **************************************************
  clock_t now = clock();

  //basic stuff
  float NchHadrons;  
  float weight;
  int processCode;
  int nMPI;
  std::vector< float > px;
  std::vector< float > py;
  std::vector< float > pz;
  std::vector< float > m;
  std::vector< int > pid;
  std::vector< int > chg;
  //DIS stuff
  float Q2,W2,x,y;  

  TFile * f = TFile::Open(Form("pp_mb_cp5_%d.root",jobNumber),"recreate");
  TTree * trackTree = new TTree("trackTree","v1");
 
  trackTree->Branch("NchHadrons",&NchHadrons);
  trackTree->Branch("weight",&weight);
  trackTree->Branch("processCode",&processCode);
  trackTree->Branch("nMPI",&nMPI);
  trackTree->Branch("px",&px);
  trackTree->Branch("py",&py);
  trackTree->Branch("pz",&pz);
  trackTree->Branch("m",&m);
  trackTree->Branch("pid",&pid);
  trackTree->Branch("chg",&chg);
  
  trackTree->Branch("Q2",&Q2);
  trackTree->Branch("W2",&W2);
  trackTree->Branch("x",&x);
  trackTree->Branch("y",&y);
  
  //******************************************** ANALYZER ********************************************************
  // Begin event loop.
  int nEvent = 100;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if( iEvent%1000 == 0 ) std::cout << iEvent << std::endl;
    if (!pythia.next()) continue;
    
    //dumping basic particle info
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal()){
	weight = pythia.info.weight();
	processCode = pythia.info.code();
	nMPI = pythia.info.nMPI();
        if(pythia.event[i].isCharged() && pythia.event[i].isHadron()) NchHadrons++;
        px.push_back(pythia.event[i].px());
        py.push_back(pythia.event[i].py());
        pz.push_back(pythia.event[i].pz());
        m.push_back(pythia.event[i].m());
        pid.push_back(pythia.event[i].id());
        chg.push_back(pythia.event[i].charge());
      }

    //used for filtering on multiplicity if desired 
    if(NchHadrons<NchCut) continue;

    //fill and clear branches
    trackTree->Fill();
    NchHadrons = 0;
    weight = 0;
    processCode = 0;
    nMPI = 0;
    px.clear();
    py.clear();
    pz.clear();
    m.clear();
    pid.clear();
    chg.clear();
    x=0;
    y=0;
    Q2=0;
    W2=0;
  }

  //end of job stuff
  pythia.stat();
  std::cout << ((float)(clock() - now))/CLOCKS_PER_SEC << " seconds" << std::endl;

  trackTree->Write();
  f->Close();
  return 0;
}
