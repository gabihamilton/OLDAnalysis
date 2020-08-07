// ----------------------------------------------------------------------------------- //
/*
  ROOT macro for illustrating the power of the Kolmogorov-Smirnov test.

  The Kolmogorov-Smirnov test (KS-test) is a general test for two distributions in 1D
  (for example to histograms) are the same. This program applies both a binned and an
  unbinned KS test, and compares it to a Chi2 test.

  References:
    http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

  Author: Troels C. Petersen (CERN)
  Email:  Troels.Petersen@cern.ch
  Date:   6th of October 2011
*/
// ----------------------------------------------------------------------------------- //

#include <vector>
#include <algorithm>
#include "TMath.h"

//------------ MIGUELS -------------//
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
#include "TRandom3.h"
#include "TStyle.h"

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


const int nbins_nu = 9;
double Nu_min;
double Nu_max;
//Double_t xbins_nu[nbins_nu+1] = {2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2};
const double limit_energy = 2.5;
int Nu_bin;

double Calculate_Modified_Xf(Float_t Shift, Float_t Nu, Float_t P, Float_t Pt, Float_t Q2, Float_t W, Float_t Zh)
{
  //cout << "Shift " << Shift << endl;
  //double xf = ((Nu + 0.9385)*(TMath::Sqrt(P*P-Pt*Pt) - TMath::Sqrt(Q2+Nu*Nu)*Zh*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0  .1395*0\.1395*W*W)/2./W));
double xf = ((Nu + 0.9385)*(TMath::Sqrt((P+Shift)*(P+Shift)-(P+Shift)*(Pt/P)*(P+Shift)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+Shift/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W));
//cout<< "Xf " << xf << endl;
return xf;
}

//-----------------------------------------------//


double sqr(double a) {
  return a*a;
}


// ----------------------------------------------------------------------------------- //
//KolmogorovSmirnovTest() {
int main(int argc, char *argv[]) {

  // ----------------------------------------------------------------------------------- //
  gROOT->Reset();

  gStyle->SetOptStat(1111);
  // gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  // gStyle->SetOptFit(0);

  gStyle->SetStatBorderSize(1);
  gStyle->SetStatFontSize(0.055);
  gStyle->SetCanvasColor(4);
  gStyle->SetPalette(1);

  // Random numbers from the Mersenne-Twister:
  TRandom3 r;
  // r.SetSeed(0);   // Initializes with the current time (i.e. random!)
  r.SetSeed(1);   // Initializes with 1, "fixing" random number array.

  //gSystem->Load("libMathCore");

  bool verbose = true;

  const int Nexp = 100;
  const Int_t Nevents = 10000; // 1000000 is the limit to segmentation fault

  double meanA  =  0.1;
  double widthA =  1.0;
  double meanB  =  0.0;
  double widthB =  1.0;


  // ------------------------------------------------------------------ //
  // Make Histograms and vectors:
  // ------------------------------------------------------------------ //

  TH1F* Hist_gaussA  = new TH1F("Hist_gaussA", "Hist_gaussA",  100, -5.0, 5.0);
  TH1F* Hist_gaussB  = new TH1F("Hist_gaussB", "Hist_gaussB",  100, -5.0, 5.0);

  TH1F* Hist_ProbMu  = new TH1F("Hist_ProbMu",  "Hist_ProbMu",   50,  0.0, 1.0);
  TH1F* Hist_ProbCSb = new TH1F("Hist_ProbCSb", "Hist_ProbCSb",  50,  0.0, 1.0);
  TH1F* Hist_ProbKSb = new TH1F("Hist_ProbKSb", "Hist_ProbKSb",  50,  0.0, 1.0);
  TH1F* Hist_ProbKS  = new TH1F("Hist_ProbKS" , "Hist_ProbKS",   50,  0.0, 1.0);


  // ------------------------------------------------------------------ //
  // Loop over generating and fitting data:
  // ------------------------------------------------------------------ //

  Int_t Ndof, igood;
  Double_t chi2, res;

  TString Nuclei_Type = argv[1];      //     TARGET TYPE
  Nu_bin = atoi(argv[2]);             //     NUMBER OF NU BIN
  Nu_min =  3.2 + Nu_bin*0.1;         //     
  Nu_max = Nu_min + 0.1;

  double limit_xf = 0.1;
  int nshifts_E = 100;
  double step_E = 1.0/1000.0;

  if(Nuclei_Type =="C" || Nuclei_Type == "Fe" || Nuclei_Type == "Pb"){
    cout << "The Nuclei type studied is " << Nuclei_Type << endl;
  }
  else std::cout << "Please choose as a first argument either C, Fe or Pb " << std::endl;

  cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;

  ostringstream nubin_label;
  nubin_label << Nu_min  <<" < #nu < " << Nu_max << " GeV ";

  //Reading Tree and Ntuple
  TFile *file;
  if(Nuclei_Type =="Pb") file = new TFile("/user/h/hamilton/ThesisProj/data/Pb_data.root");
  else if(Nuclei_Type =="Fe") file = new TFile("/user/h/hamilton/ThesisProj/data/Fe_data.root");
  else if(Nuclei_Type =="C") file = new TFile("/user/h/hamilton/ThesisProj/data/C_data.root");
  else std::cout<< "Problem, please choose either Pb, Fe or C " << std::endl;
  //  if(!file) return;
  cout<< "The cut on Xf is " << limit_xf << endl;


  TFile* fout = new TFile(Form("OUTPUT/fout_"+Nuclei_Type+"_nubin%d.root",Nu_bin),"RECREATE");
  //cout<< "OK " << endl;

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

  const Int_t nentries = tree->GetEntries();
  double energy_shift = 0;
  ////////////////////////////////////////////////////////////
  ////////FILLING, THE HISTOGRAMS.
  ///////////////////////////////////////////////////////////
  //-------------------------------------------------------
  //  cout<< "OK " << endl;

  //  Nexp = nshifts_E;
  //Nevents = nentries;

  for (int iexp=0; iexp < Nexp; iexp++) {

    vector<double> x_gaussA;
    vector<double> x_gaussB;
    //cout<< "OK " << endl;

    Double_t* x_gaussA_array = new Double_t[nentries];
    Double_t* x_gaussB_array = new Double_t[nentries];
    cout<< "OK " << endl;

    double sumA[3] = {0.0, 0.0, 0.0};
    double sumB[3] = {0.0, 0.0, 0.0};

    // Generate data:
    // --------------
    for (int i=0; i < nentries; i++) {
      if(i%5000000==0) std::cout<< "Processing event number " << i/1000000.0 << "M "<< std::endl;
      tree->GetEntry(i);
      //Apply Cuts bin in Nu
      if(Nu > Nu_max || Nu < Nu_min) continue;
      //Apply cuts for Deuterium, fill with Energy = Zh*Nu
      if(TargType==1 && Xf> limit_xf && Zh*Nu < limit_energy){
	double xA = (Zh*Nu);
	sumA[1] += xA;
	sumA[2] += xA*xA;
	x_gaussA.push_back(xA);
	Hist_gaussA->Fill(xA);
      }

      //Apply cuts for Nuclei, fill with Energy = Zh*Nu + energy shift.

      if(TargType!=2) continue;
      energy_shift = step_E*iexp; 
      double Xf_Nuclei = Calculate_Modified_Xf( energy_shift, Nu, P, Pt,  Q2,  W, Zh);
      if(Xf_Nuclei> limit_xf && Zh*Nu + energy_shift < limit_energy ){
	double xB = Zh*Nu+energy_shift;
	sumB[1] += xB;
	sumB[2] += xB*xB;
	x_gaussB.push_back(xB);
	Hist_gaussB->Fill(xB);
      }
    }//END LOOP OVER ENTRIES
  std::cout<< "END LOOP OVER ENTRIES" << std::endl;

 
  //  cout<< "OK " << endl;

    // Test distributions:
    // -------------------

    // Test mean and width:
  sumA[1] = sumA[1] / double(Nevents);
  sumA[2] = sqrt(sumA[2] / double(Nevents) - sumA[1]*sumA[1]);
  sumB[1] = sumB[1] / double(Nevents);
  sumB[2] = sqrt(sumB[2] / double(Nevents) - sumB[1]*sumB[1]);
  double dMean = sumA[1] - sumB[1];
  double dMeanError = sqrt(sqr(sumA[2])/double(Nevents) + sqr(sumB[2])/double(Nevents));
  double pMean = TMath::Erfc(dMean / dMeanError / sqrt(2.0)) / 2.0;
  Hist_ProbMu->Fill(pMean);
  //  cout<< "OK " << endl;

  // Unbinned Kolmogorov-Smirnov Test (requires sorting):
  sort(x_gaussA.begin(), x_gaussA.end());
  sort(x_gaussB.begin(), x_gaussB.end());
  cout<< "OK " << endl;    // OK
  copy(x_gaussA.begin(), x_gaussA.end(), x_gaussA_array);
  copy(x_gaussB.begin(), x_gaussB.end(), x_gaussB_array);
  //x_gaussA_array = &x_gaussA[0];
  //x_gaussB_array = &x_gaussB[0];
  cout<< "OK " << endl;
  double pKS = TMath::KolmogorovTest(Nevents, x_gaussA_array, Nevents, x_gaussB_array, "");
  cout<< "OK " << endl;
  Hist_ProbKS->Fill(pKS);
  
  // Binned Chi2 and Kolmogorov-Smirnov Test:
  double pCSbinned = Hist_gaussA->Chi2Test(Hist_gaussB, "UU");
  double pKSbinned = Hist_gaussA->KolmogorovTest(Hist_gaussB, "UU");
  Hist_ProbCSb->Fill(pCSbinned);
  Hist_ProbKSb->Fill(pKSbinned);
  cout<< "OK " << endl;
  
  if (verbose)
    printf(" %4d:  pMean: %6.4f   pCSbinned: %6.4f   pKSbinned: %6.4f   pKS: %6.4f \n",
	   iexp, pMean, pCSbinned, pKSbinned, pKS);
  cout<< "OK " << endl;
  
  if (Nexp > 1) Hist_gaussA->Reset("ICE");
  if (Nexp > 1) Hist_gaussB->Reset("ICE");
  }  //END LOOP OVER SHIFTS



  // ------------------------------------------------------------------ //
  // Show the distribution of fitting results:
  // ------------------------------------------------------------------ //

  if (Nexp == 1) {

    // Make a white canvas:
    TCanvas *c0 = new TCanvas("c0","",220,20,500,300);
    c0->SetFillColor(0);
    
    Hist_gaussA->SetLineWidth(2);
    Hist_gaussA->SetLineColor(2);
    Hist_gaussA->Draw();
    
    Hist_gaussB->SetLineWidth(2);
    Hist_gaussB->SetLineColor(4);
    Hist_gaussB->Draw("same");
    
    c0->Update();
    c0->SaveAs("Dist.pdf");
  }


  // ------------------------------------------------------------------ //
  // Show the distribution of fitting results:
  // ------------------------------------------------------------------ //

  else {

    // Make a white canvas:
    TCanvas *c1 = new TCanvas("c1","",250,70,500,800);
    c1->SetFillColor(0);
    c1->Divide(1,4);
    
    c1->cd(1);
    Hist_ProbMu->SetMinimum(0.0);
    Hist_ProbMu->SetLineWidth(2);
    Hist_ProbMu->SetLineColor(4);
    Hist_ProbMu->Draw();
    
    c1->cd(2);
    Hist_ProbCSb->SetMinimum(0.0);
    Hist_ProbCSb->SetLineWidth(2);
    Hist_ProbCSb->SetLineColor(4);
    Hist_ProbCSb->Draw();
    
    c1->cd(3);
    Hist_ProbKSb->SetMinimum(0.0);
    Hist_ProbKSb->SetLineWidth(2);
    Hist_ProbKSb->SetLineColor(4);
    Hist_ProbKSb->Draw();
    
    c1->cd(4);
    Hist_ProbKS->SetMinimum(0.0);
    Hist_ProbKS->SetLineWidth(2);
    Hist_ProbKS->SetLineColor(4);
    Hist_ProbKS->Draw();
    
    c1->Update();
    c1->SaveAs("TestDistNEW.pdf");
  }

}


//---------------------------------------------------------------------------------- 
/*
Questions:
----------
 1) First run the program to display the two distributions A and B, when
    - They are the same.
    - The mean of A is increased.
    - The width of A is enlarged.
    Get a feel for how much you need to change the distribution, before you can
    by eye see that they are not the same.
    Could you for the test of the means have calculated this? Do so, and see if
    it somewhat matches you number from above!

 2) Now run the tests 100 times, where A and B are unit Gaussians and thus identical.
    How should the distributions of the test probabilities come out? And is this the
    case?

 3) Repeat the changes in question 1), and see which tests "reacts" most to these
    modifications.
    How much of a change in the mean is required for 95% of the tests (of each kind)
    to give a probability below 5%? How much is required for the width?

 4) Obviously, the test of the means is not sensitive the a change in the width.
    Make such a test yourself by calculating the widths and the uncertainty on the
    widths (see Cowan, end of chapter 5).
    Note that the uncertainty on the widths is close to that of the means! That is
    generally good to know :-)

 5) The unbinned Kolmogorov-Smirnov test has the great advantage that it can handle
    ANY distribution (even the Cauchy distribution - remind yourself of that one!).
    Try to test different distributions than the Gaussian one (e.g. exponential,
    binomial, etc.), and see how the tests performs.


Advanced questions:
-------------------
 1) Implement in ROOT the following tests:
     - Lilliefors test
     - Shapiro-Wilk test
     - Anderson-Darling test
     - Cram
*/
