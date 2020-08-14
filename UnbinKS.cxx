#include "TROOT.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCut.h"
#include "TChain.h"
#include "TApplication.h"
#include "math.h"
#include "TLegendEntry.h"
#include "TGraph.h"


#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
using namespace std;

//------PARAMETERS-----//
const int N_Nu = 5;  // number of Nu bins
const double E_max = 2.5;
const double E_min = 0.;
const double limit_xf = 0.1;
const int nshift_E = 99;
const double step_E = 1.1/1000.0;
const int nbins = 100;

int Nu_bin;      //number of Nu bin
TString Nuclei_Type;


//---------FUNCTIONS-------//

//-------MODIFIED FEYNMAN X--------//
double Calculate_Modified_Xf(Float_t Shift, Float_t Nu, Float_t P, Float_t Pt, Float_t Q2, Float_t W, Float_t Zh)
{
  //cout << "Shift " << Shift << endl;
  //double xf = ((Nu + 0.9385)*(TMath::Sqrt(P*P-Pt*Pt) - TMath::Sqrt(Q2+Nu*Nu)*Zh*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0  .1395*0.1395*W*W)/2./W));
  double xf = ((Nu + 0.9385)*(TMath::Sqrt((P+Shift)*(P+Shift)-(P+Shift)*(Pt/P)*(P+Shift)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+Shift/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W));
  //cout<< "Xf " << xf << endl;
  return xf;
}


int main(int argc, char *argv[]){
 
  Nuclei_Type = (TString) argv[1];
  Nu_bin = atoi(argv[2]);

  double Nu_min = 3.2 + Nu_bin*0.2; 
  double Nu_max = Nu_min + 0.2;;

  cout << "The Nuclei type studied is " << Nuclei_Type << endl;
 
  cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;

  ostringstream nubin_label;
  nubin_label << Nu_min  <<" < #nu < " << Nu_max << " GeV ";

  cout<< "The cut on Xf is " << limit_xf << endl;

  //------Opening data files-----//
  TFile *file = new TFile(Form("/Users/ghamiltonm/Documents/data/" + Nuclei_Type + "_data.root"));
  //-----Creating output file-----//
  TFile* fout = new TFile(Form("OUTPUT/test3_"+Nuclei_Type+"_nubin%d.root",Nu_bin),"RECREATE");
  //-----Opening TTree----//
  TTree* tree = (TTree*)file->Get("ntuple_data");
  
  //Reading Branches with appropiate variables. 
  Float_t TargType;
  Float_t Q2;
  Float_t Nu;
  Float_t Xb;
  Float_t W;
  Float_t SectorEl;
  Float_t ThetaPQ;
  Float_t PhiPQ;
  Float_t Zh;
  Float_t Pt;
  Float_t W2p;
  Float_t Xf;
  Float_t T;
  Float_t P;
  Float_t T4;
  Float_t deltaZ;
  Float_t NmbPion;
  tree->SetBranchAddress("TargType",&TargType);
  tree->SetBranchAddress("Q2",&Q2);
  tree->SetBranchAddress("Nu",&Nu);
  tree->SetBranchAddress("Xb",&Xb);
  tree->SetBranchAddress("W",&W);
  tree->SetBranchAddress("SectorEl",&SectorEl);
  tree->SetBranchAddress("ThetaPQ",&ThetaPQ);
  tree->SetBranchAddress("PhiPQ",&PhiPQ);
  tree->SetBranchAddress("Zh",&Zh);
  tree->SetBranchAddress("Pt",&Pt);
  tree->SetBranchAddress("W2p",&W2p);
  tree->SetBranchAddress("Xf",&Xf);
  tree->SetBranchAddress("T",&T);
  tree->SetBranchAddress("P",&P);
  tree->SetBranchAddress("T4",&T4);
  tree->SetBranchAddress("deltaZ",&deltaZ);
  tree->SetBranchAddress("NmbPion",&NmbPion);

  Int_t nentries = tree->GetEntries();
  double energy_shift;
  int i_shift=0;


  vector<vector<double>> x_gaussS(100);
  vector<double> x_gaussD;

  //double* x_gaussD_array = new double[tree->GetEntries()];
  //double* x_gaussS_array = new double[tree->GetEntries()];

  TH1F* Hist_ProbKS  = new TH1F("Hist_ProbKS" , "Hist_ProbKS",   50,  0.0, 1.0);

  std::cout<<"Entering the loop over entries "<< nentries <<std::endl;

  for(Int_t j=1; j<= nentries; j++){
    if(j%5000000==0) std::cout<< "Processing event number " << j/1000000.0 << "M "<< std::endl;
    tree->GetEntry(j);
    //Apply Cuts bin in Nu
    if(Nu> Nu_max || Nu < Nu_min) continue;
    //Apply cuts for Deuterium, fill with Energy = Zh*Nu
    if(TargType==1 && Xf> limit_xf && Zh*Nu < E_max){
      x_gaussD.push_back(Zh*Nu);
     }
    //Apply cuts for Nuclei, fill with Energy = Zh*Nu + energy shift. 
    //cout << "OK" <<endl;
    if(TargType!=2) continue;

    //x_gaussS.push_back(std::vector<double>());
    for (int i_shift = 0; i_shift < nshift_E; ++i_shift)
    {
      energy_shift = step_E*i_shift; //In GeV // 5 MeV steps.                                                                                                         
      double Xf_Nuclei = Calculate_Modified_Xf(energy_shift, Nu, P, Pt,  Q2,  W, Zh);
      if(Xf_Nuclei> limit_xf && Zh*Nu + energy_shift < E_max ){
        x_gaussS[i_shift].push_back((Zh*Nu)+energy_shift);
      }
    }
  
  }
  std::cout<< "END LOOP OVER ENTRIES" << std::endl;
  //double pKS[100];
  TGraph* gr = new TGraph();
  //double x[], y[];
  for (int i_shift = 0; i_shift < nshift_E; ++i_shift)
  {
    double* x_gaussD_array = new double[tree->GetEntries()];
    double* x_gaussS_array = new double[tree->GetEntries()];

    sort(x_gaussD.begin(), x_gaussD.end());
    sort(x_gaussS[i_shift].begin(), x_gaussS[i_shift].end());

    copy(x_gaussD.begin(), x_gaussD.end(), x_gaussD_array);
    copy(x_gaussS[i_shift].begin(), x_gaussS[i_shift].end(), x_gaussS_array);

    double pKS = TMath::KolmogorovTest(nentries, x_gaussD_array, nentries, x_gaussS_array, "");
    Hist_ProbKS->Fill(pKS);
    //energy_shift = nshift_E*i_shift;
    gr->SetPoint(i_shift, nshift_E*i_shift*1000, -1*TMath::Log10(pKS));
  }
  //TGraph* gr = new TGraph(nshift_E,energy_shift,pKS);

  fout->cd();
  Hist_ProbKS->Write();
  gr->SetName("G");
  gr->Write();
  //Nuclei->Write();
  //fout->Print();
  std::cout<<" ABOUT TO CLOSE " << std::endl;
  fout->Close();
  
  std::cout<< " BYE BYE " << std::endl;
  return 1;
}
