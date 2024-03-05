#include <string>
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TCanvas.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include "coordinateTools.h"

void analyze(){

  float NchHadrons;
  std::vector< float > * px = 0;
  std::vector< float > * py = 0;
  std::vector< float > * pz = 0;
  std::vector< float > * m = 0;
  std::vector< int > * pid = 0;
  std::vector< int > * chg = 0;

  TFile * inputFile = TFile::Open("/Users/davidlw/pythia8310/examples/pp_mb_cp5_1.root","read");
  TTree * t = (TTree*)inputFile->Get("trackTree");

  t->SetBranchAddress("px",&px);
  t->SetBranchAddress("py",&py);
  t->SetBranchAddress("pz",&pz);
  t->SetBranchAddress("m",&m);
  t->SetBranchAddress("pid",&pid); 
  t->SetBranchAddress("chg",&chg); 

  //event loop
  for(int i = 0; i<t->GetEntries(); i++){

    t->GetEntry(i);

    std::cout<<"There are "<<px->size()<<" particles in this event"<<std::endl;
    for(unsigned int j = 0; j<px->size(); j++) std::cout<<"px of "<<j<<"th particle is: " << px->at(j) << std::endl;
  }


  inputFile->Close();
}
