#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TCut.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TTree.h"
#include "TNtuple.h"
#include "math.h"
#include "TMath.h"
#include "TLegendEntry.h"
#include "TGraph.h"

using namespace std;

Double_t Q2_min;
Double_t Q2_max;
Double_t Nu_min;
Double_t Nu_max;
Double_t Phi_min;
Double_t Phi_max;
Double_t Pt2_min;
Double_t Pt2_max;
Double_t E_min;
Double_t E_max;

Double_t delta_Q2;
Double_t delta_Nu;
Double_t delta_Phi;
Double_t delta_Pt2;

Int_t N_Q2;
Int_t N_Nu;
Int_t N_Phi;
Int_t N_Pt2;
Int_t nbins;

TString f_location;
TString fd_ext;
TString fs_ext;
TString Target;
TString Nuclei_Type;

Int_t Nu_bin;
Int_t Bundle;
Int_t bundle_size;


/*
string toString(int i)
{
	stringstream ss;
	ss.str("");
	ss<<i;
	return ss.str();
}
*/


int main(int argc, char **argv)
{
  const Double_t Q2_min = 1.;
  const Double_t Q2_max = 4.;
  const Int_t N_Q2= 6.;//6
  
  const Double_t Phi_min = -180.;
  const Double_t Phi_max = 180.;
  const Int_t N_Phi = 12.; //12
  
  const Double_t Pt2_min = 0.;
  const Double_t Pt2_max = 1.5;
  const Int_t N_Pt2 = 6.;//6
    

  const Double_t E_min = 0.;	
  const Double_t E_max = 2.5;  

  const Int_t nbins = 100;      // number of energy bins
  
  const Int_t nSimuFiles = 4;

  TString f_location = "/user/h/hamilton/ThesisProj/data/"; // The location of data and simulation files with ntuples inside
  TString fd_ext = "_data.root";
  TString fs_ext = "_simul.root";

  //  Target =  (TString) argv[1];   // 1 for liquid or 2 for solid

  Nuclei_Type = (TString) argv[1]; // C for Carbon, Fe for Iron and Pb for Lead

  Nu_bin = (Int_t) std::stoi(argv[2]);

  Bundle = 0;//(Int_t) std::stoi(argv[4]); // in number of bundle

  bundle_size = 100;//(Int_t) std::stoi(argv[5]);


  Double_t Nu_min = 3.2 + Nu_bin*0.1;
  Double_t Nu_max = Nu_min + 0.1;

  int shift_min = Bundle*bundle_size;
  int shift_max = Bundle*bundle_size + bundle_size - 1;

  double limit_xf = 0.1;
  int nshift_E = 99;         // total number of shifts in Energy
  double step_E = 1.0/1000.0;   // size of step of energy shift

  
  delta_Q2 = (Q2_max-Q2_min)/N_Q2;
  delta_Phi = (Phi_max-Phi_min)/N_Phi;
  delta_Pt2 = (Pt2_max-Pt2_min)/N_Pt2;
 
  
  /*---------------------------CUTS-----------------------------------------*/

  TCut Target_cutD = "TargType==1"; //cut for Deuterium
  TCut Target_cutS = "TargType==2"; //cut for Solid Target
  //Simulation Cuts
  TCut Q2_cut_S = Form("Q2>%f && Q2<%f", Q2_min, Q2_max);                 
  TCut Nu_cut_S = Form("Nu>%f && Nu<%f", Nu_min, Nu_max);
  TCut Phi_cut_S = Form("PhiPQ>%f && PhiPQ<%f", Phi_min, Phi_max);
  TCut Pt2_cut_S = Form("Pt*Pt>%f && Pt*Pt<%f", Pt2_min, Pt2_max);

  TCut cuts_simulD = Target_cutD&&Q2_cut_S&&Nu_cut_S&&Phi_cut_S&&Pt2_cut_S;
  TCut cuts_simulS = Target_cutS&&Q2_cut_S&&Nu_cut_S&&Phi_cut_S&&Pt2_cut_S;
  
  TCut xf_cut = "Xf>0.1"; //  Typical xF cut
  TCut Q2_cut, Nu_cut, Phi_cut, Pt2_cut, cuts_loop, xf_mod;  //Loops Cuts                           

  //-----------------------OBTAINING THE NTUPLES-------------------------//

  TChain *data = new TChain("ntuple_data");
  data->Add(Form(f_location + Nuclei_Type + fd_ext));
  data->SetBranchStatus("*",0);
  data->SetBranchStatus("Q2",1);
  data->SetBranchStatus("Nu",1);
  data->SetBranchStatus("PhiPQ",1);
  data->SetBranchStatus("Pt",1);
  data->SetBranchStatus("Xf",1);
  data->SetBranchStatus("Zh",1);
  data->SetBranchStatus("TargType",1);
  data->SetBranchStatus("W",1);
  data->SetBranchStatus("P",1);

  TChain *reconstructed_D = new TChain("ntuple_accept");
  for(Int_t q = 0; q < nSimuFiles; q++){
    reconstructed_D->Add(Form(f_location + "D%d"+ fs_ext, q+1));
  }
  reconstructed_D->SetBranchStatus("*",0);
  reconstructed_D->SetBranchStatus("Q2",1);
  reconstructed_D->SetBranchStatus("Nu",1);
  reconstructed_D->SetBranchStatus("PhiPQ",1);
  reconstructed_D->SetBranchStatus("Pt",1);
  reconstructed_D->SetBranchStatus("Xf",1);
  reconstructed_D->SetBranchStatus("Zh",1);
  reconstructed_D->SetBranchStatus("TargType",1);
  reconstructed_D->SetBranchStatus("W",1);
  reconstructed_D->SetBranchStatus("P",1);

  reconstructed_D->Draw(">>list_accD",cuts_simulD,"goff");
  reconstructed_D->SetEventList((TEventList*)gDirectory->Get("list_accD"));

  TChain *reconstructed_S = new TChain("ntuple_accept");
  for(Int_t e = 0; e < nSimuFiles; e++){
    reconstructed_S->Add(Form(f_location + Nuclei_Type +"%d" + fs_ext, e+1));
    //    reconstructed_S->Add(f_location + Nuclei_Type + toString(e+1) + fs_ext);
  }
  reconstructed_S->SetBranchStatus("*",0);
  reconstructed_S->SetBranchStatus("Q2",1);
  reconstructed_S->SetBranchStatus("Nu",1);
  reconstructed_S->SetBranchStatus("PhiPQ",1);
  reconstructed_S->SetBranchStatus("Pt",1);
  reconstructed_S->SetBranchStatus("Xf",1);
  reconstructed_S->SetBranchStatus("Zh",1);
  reconstructed_S->SetBranchStatus("TargType",1);
  reconstructed_S->SetBranchStatus("W",1);
  reconstructed_S->SetBranchStatus("P",1);

  reconstructed_S->Draw(">>list_accS",cuts_simulS,"goff");
  reconstructed_S->SetEventList((TEventList*)gDirectory->Get("list_accS"));


  TChain *thrown_D = new TChain("ntuple_thrown");
  for(Int_t w = 0; w < nSimuFiles; w++){
    thrown_D->Add(Form(f_location+"D%d" + fs_ext, w+1));
//thrown_D->Add(f_location + "D" + toString(w+1) + fs_ext);
  }
  thrown_D->SetBranchStatus("*",0);
  thrown_D->SetBranchStatus("Q2",1);
  thrown_D->SetBranchStatus("Nu",1);
  thrown_D->SetBranchStatus("PhiPQ",1);
  thrown_D->SetBranchStatus("Pt",1);
  thrown_D->SetBranchStatus("Xf",1);
  thrown_D->SetBranchStatus("Zh",1);
  thrown_D->SetBranchStatus("TargType",1);
  thrown_D->SetBranchStatus("W",1);
  thrown_D->SetBranchStatus("P",1);

  thrown_D->Draw(">>list_thrD",cuts_simulD,"goff");
  thrown_D->SetEventList((TEventList*)gDirectory->Get("list_thrD"));
  
  TChain *thrown_S = new TChain("ntuple_thrown");
  for(Int_t r = 0; r < nSimuFiles; r++){
    thrown_S->Add(Form(f_location + Nuclei_Type + "%d" + fs_ext, r+1));
		       //    thrown_S->Add(f_location + Nuclei_Type  + toString(r+1) + fs_ext);
  }
  thrown_S->SetBranchStatus("*",0);
  thrown_S->SetBranchStatus("Q2",1);
  thrown_S->SetBranchStatus("Nu",1);
  thrown_S->SetBranchStatus("PhiPQ",1);
  thrown_S->SetBranchStatus("Pt",1);
  thrown_S->SetBranchStatus("Xf",1);
  thrown_S->SetBranchStatus("Zh",1);
  thrown_S->SetBranchStatus("TargType",1);
  thrown_S->SetBranchStatus("W",1);
  thrown_S->SetBranchStatus("P",1);

  thrown_S->Draw(">>list_thrS",cuts_simulS,"goff");
  thrown_S->SetEventList((TEventList*)gDirectory->Get("list_thrS"));


  //  CREATING THE FILE
  TFile *plots = new TFile(Form("output/fout_"+Nuclei_Type+"_nubin%d.root",Nu_bin),"RECREATE");
  //cout << "OK" << endl;

  //--------CREATING HISTOGRAMS--------//

  TH1F *D = new TH1F("D", "D", nbins, E_min, E_max);     // Deuterium Histrograms

  double energy_shift;

  std::map<int,TH1F*> histograms;

  for(int i = shift_min; i<=shift_max; i++){
    histograms[i] = new TH1F(Form("Nuclei_shift%d",i),Form("Nuclei_shift%d",i),nbins,E_min,E_max);
    histograms[i]->Sumw2();
  }

  //cout << "OK" << endl;


//------------ACCEPTANCE CORRECTION------------//

Nu_cut = Form("Nu>%f && Nu<%f", Nu_min, Nu_max);  //  Cut for Nu bin

 for (int j = 0; j < N_Q2; j++){

   Q2_cut = Form("Q2>%f && Q2<%f", Q2_min + j*delta_Q2 , Q2_min + (j+1)*delta_Q2);  //  Cut for Q2 bin			

   //----ACCEPTANCE FOR DEUTERIUM-----/
   cuts_loop=Q2_cut&&Nu_cut&&Target_cutD&&xf_cut;

   TH1F *data_histo = new TH1F("data_histo","",nbins,E_min,E_max);
   TH1F *thrown_histo = new TH1F("thrown_histo","",nbins,E_min,E_max);
   TH1F *reconstructed_histo = new TH1F("reconstructed_histo","",nbins,E_min,E_max);

   data->Draw("Zh*Nu>>data_histo",cuts_loop,"goff");
   data_histo->Sumw2();

   thrown_D->Draw("Zh*Nu>>thrown_histo",cuts_loop,"goff");
   thrown_histo->Sumw2();
   TH1F *hDT = (TH1F*) thrown_histo->Clone();
   hDT->SetName(Form("thrown_histoD%d%d", Nu_bin, j));
   hDT->Write();


   //cout << "OK" << endl;

   reconstructed_D->Draw("Zh*Nu>>reconstructed_histo",cuts_loop,"goff");
   reconstructed_histo->Sumw2();
   TH1F *hDR = (TH1F*) reconstructed_histo->Clone();
   hDR->SetName(Form("reconstructed_histoD%d%d", Nu_bin, j));
   hDR->Write();
      
   TH1F *acceptance_histo = new TH1F("acceptance_histo","",nbins,E_min,E_max);
   acceptance_histo->Divide(reconstructed_histo,thrown_histo,1,1,"B");
   TH1F *hDA = (TH1F*) acceptance_histo->Clone();
   hDA->SetName(Form("acceptance_histoD%d%d", Nu_bin, j));
   hDA->Write();

   TH1F *acceptance_correction_histo = new TH1F("acceptance_correction_histo","",nbins,E_min,E_max);
   acceptance_correction_histo->Divide(data_histo,acceptance_histo,1,1);
   TH1F *hDC = (TH1F*) acceptance_correction_histo->Clone();
   hDC->SetName(Form("acc_correctedD%d%d", Nu_bin, j));
   hDC->Write();

   // Filling the Deuterium histogram
   D->Add(acceptance_correction_histo, 1);

   delete acceptance_correction_histo;
   delete acceptance_histo;
   delete data_histo;
   delete thrown_histo;
   delete reconstructed_histo;
   delete hDT;
   delete hDR;
   delete hDA;
   delete hDC;

   //---ACCEPTANCE FOR SOLID TARGET::::ENERGY SHIFT LOOP----//
   for(int shift = shift_min; shift <= shift_max; shift++){
     //cout << "OK" << endl;     
     energy_shift = step_E*shift;

     TCut xf_mod = Form("((Nu + 0.9385)*(TMath::Sqrt((P+%f)*(P+%f)-(P+%f)*(Pt/P)*(P+%f)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+%f/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W))>0.1", energy_shift, energy_shift, energy_shift, energy_shift, energy_shift); \

     cuts_loop=Q2_cut&&Nu_cut&&Phi_cut&&Pt2_cut&&Target_cutS&&xf_mod;   //NO OLVIDARSE EL XF MOD CUT
     TCut cut = Q2_cut&&Nu_cut&&Target_cutS;

     TH1F *data_histo = new TH1F("data_histo","",nbins,E_min,E_max);
     TH1F *thrown_histo = new TH1F("thrown_histo","",nbins,E_min,E_max);
     TH1F *reconstructed_histo = new TH1F("reconstructed_histo","",nbins,E_min,E_max);
     //cout<<"OK"<<endl;
     TH1F *Data_xF = new TH1F("Data_xF", "", nbins, -1.5,1.5);
     data->Draw(Form("((Nu + 0.9385)*(TMath::Sqrt((P+%f)*(P+%f)-(P+%f)*(Pt/P)*(P+%f)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+%f/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W))>>Data_xF", energy_shift, energy_shift, energy_shift, energy_shift, energy_shift), cut, "goff");
     Data_xF->Sumw2();
     TH1F *Xf = (TH1F*) Data_xF->Clone();
     Xf->SetName(Form("Xf_shift%d_%d%d", shift, Nu_bin, j));
     Xf->Write();
     
     data->Draw(Form("(Nu*Zh)+%f>>data_histo", energy_shift), cuts_loop, "goff");
     data_histo->Sumw2();

     thrown_S->Draw(Form("(Nu*Zh)+%f>>thrown_histo", energy_shift), cuts_loop, "goff");
     thrown_histo->Sumw2();
     TH1F *hT = (TH1F*) thrown_histo->Clone();
     hT->SetName(Form("thrown_histo"+Nuclei_Type+"_shift%d_%d%d", shift, Nu_bin, j));
     hT->Write();
     
     reconstructed_S->Draw(Form("(Nu*Zh)+%f>>reconstructed_histo", energy_shift), cuts_loop, "goff");
     reconstructed_histo->Sumw2();
     TH1F *hR = (TH1F*) reconstructed_histo->Clone();
     hR->SetName(Form("reconstructed_histo"+Nuclei_Type+"_shift%d_%d%d", shift, Nu_bin, j));
     hR->Write();
     

     //cout <<"OK" << endl;
     // HISTOGRAMAS WITH ACCEPTANCE

     TH1F *acc_histo = new TH1F("acc_histo", "", nbins, E_min, E_max);
     acc_histo->Divide(reconstructed_histo, thrown_histo, 1, 1, "B");
     TH1F *hA = (TH1F*) acceptance_histo->Clone();
     hA->SetName(Form("acceptance_histo"+Nuclei_Type+"_shift%d_%d%d", Nu_bin, j));
     hA->Write();
     
     TH1F *acc_corr_histo = new TH1F("acc_corr_histo", "", nbins, E_min, E_max);
     acc_corr_histo->Divide(data_histo, acc_histo, 1, 1);
     TH1F *hC = (TH1F*) acc_corr_histo->Clone();
     hC->SetName(Form("acc_corrected"+Nuclei_Type+"_shift%d_%d%d", shift,  Nu_bin, j));
     hC->Write();
     

     histograms[shift]->Add(acc_corr_histo, 1);


     delete data_histo;
     delete thrown_histo;
     delete reconstructed_histo;
     delete acc_histo;
     delete acc_corr_histo;
     delete hT;
     delete hR;
     delete hA;
     delete Xf;
     delete Data_xF;
     delete hC;
     cout << "OK" <<endl;
   }
 }
 


 //WRITING HISTOS
 
 D->Write();

 for(int i = shift_min; i<= shift_max; i++){
   histograms[i]->Write();
 }


 //Saving the Nu and dE bins used in analysis
 TGraph* g_Ebins = new TGraph();
 TGraph* g_Nubins = new TGraph();
 TGraph* g_Xfcut = new TGraph();

 g_Nubins->SetPoint(0,0,Nu_min);
 g_Nubins->SetPoint(1,1,Nu_max);
 g_Xfcut->SetPoint(0,0,limit_xf);

 for(Int_t i_shift = 0; i_shift <= nshift_E; i_shift++){
   energy_shift = step_E*i_shift; //In GeV // 1 MeV steps.
   g_Ebins->SetPoint(i_shift,i_shift,energy_shift);
 }


 //----------------WRITE FINAL HISTOS-----------------------//
/* TH1F *data_histo_final = new TH1F("data_histo_final","Data",nbins,E_min,E_max);
TH1F *reconstructed_histo_final = new TH1F("reconstructed_histo_final","Reconstruction",nbins,E_min,E_max);
 TH1F *thrown_histo_final = new TH1F("thrown_histo_final","Thrown",nbins,E_min,E_max);

 data->Draw("Zh*Nu>>data_histo_final",cuts_simulS,"goff");
 data_histo_final->Write();

 reconstructed->Draw("Zh*Nu>>reconstructed_histo_final",cuts_simulS,"goff");
 reconstructed_histo_final->Write();

 thrown->Draw("Zh*Nu>>thrown_histo_final",cuts_simulS,"goff");
 thrown_histo_final->Write();
*/

 g_Ebins->Write("g_Ebins");
 g_Nubins->Write("g_Nubins");
 g_Xfcut->Write("g_Xfcut");

 //----------------DELETING FINAL HISTOS------------------//
 //delete data_histo_final;
 //delete reconstructed_histo_final;
 //delete thrown_histo_final;
 delete g_Ebins;
 delete g_Nubins;
 delete g_Xfcut;



//-----------------CLOSING THE FILE-----------------------//
plots->Close();

return 0;
}
