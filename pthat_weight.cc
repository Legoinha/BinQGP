#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include "TArrayD.h"

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));
Double_t ScaleX(Double_t x);

void pthat_weight(){

  TString input_file_mc_pp   = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewBDTCut/BPNew/BPMC.root";

  TString input_file_mc_PbPb   = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/2018PbPb/BPPbPbMCFlat.root";

  TFile *fin_mc_pp = new TFile(input_file_mc_pp);
  TFile *fin_mc_PbPb = new TFile(input_file_mc_PbPb);

  TTree* t_mc_pp = (TTree*)fin_mc_pp->Get("ntKp");
  TTree* t_mc_PbPb = (TTree*)fin_mc_PbPb->Get("ntKp");

  float pthat_pp;
  float weight; 
  float pthatweight;
  float pthat_PbPb;

  float pt_pp;
  float pt_PbPb;

  t_mc_pp->SetBranchAddress("pthat", &pthat_pp);
  t_mc_pp->SetBranchAddress("weight", &weight);
  t_mc_pp->SetBranchAddress("Bpt", &pt_pp);

  t_mc_PbPb->SetBranchAddress("pthatweight", &pthatweight);
  t_mc_PbPb->SetBranchAddress("pthat", &pthat_PbPb);
  t_mc_PbPb->SetBranchAddress("Bpt", &pt_PbPb);

  TH1D* pthat_histo_pp = new TH1D("pthat_histo_pp", "pthat_histo_pp", 60, 0, 1200);
  TH1D* weight_histo_pp = new TH1D("weight_histo_pp", "weight_histo_pp", 60, 0, 0.0002);
  TH1D* pt_histo_pp = new TH1D("pt_histo_pp", "pt_histo_pp", 60, 0, 100);
  TH1D* pt_histo_pp_weight = new TH1D("pt_histo_pp_weight", "pt_histo_pp_weight", 60, 0, 100);

  TH1D* pthatweight_PbPb = new TH1D("pthatweight_PbPb", "pthatweight_PbPb", 60, 0, 150);
  TH1D* pthat_histo_PbPb = new TH1D("pthat_histo_PbPb", "pthat_histo_PbPb", 60, 0, 1200);
  TH1D* pt_histo_PbPb = new TH1D("pt_histo_PbPb", "pt_histo_PbPb", 60, 0, 100);
  TH1D* pt_histo_PbPb_weight = new TH1D("pt_histo_PbPb_weight", "pt_histo_PbPb_weight", 60, 0, 100);

  TH1D* histo_pp = new TH1D("histo_pp", "histo_pp", 60, 0, 400);
  TH1D* histo_PbPb = new TH1D("histo_PbPb", "histo_PbPb", 60, 0, 400);

  int entries_pp = t_mc_pp->GetEntries();

  for(int i = 0; i < entries_pp; i++){
     t_mc_pp->GetEntry(i);
//     if(pt_pp > 5){
       pthat_histo_pp->Fill(pthat_pp);
       weight_histo_pp->Fill(weight);
       histo_pp->Fill(pthat_pp,weight);
//     }
     //pt_histo_pp->Fill(pt_pp);
     //pt_histo_pp_weight->Fill(pt_pp,weight);
  }

  int entries_PbPb = t_mc_PbPb->GetEntries();

  for(int i = 0; i < entries_PbPb; i++){
     t_mc_PbPb->GetEntry(i);
//     if(pt_PbPb > 5){
       pthat_histo_PbPb->Fill(pthat_PbPb);
       pthatweight_PbPb->Fill(pthatweight);
       histo_PbPb->Fill(pthat_PbPb,pthatweight);
//     }
  }

  histo_pp->Scale(1.0/histo_pp->Integral());
  histo_PbPb->Scale(1.0/histo_PbPb->Integral());

  pt_histo_pp->Scale(1.0/pt_histo_pp->Integral());
  pt_histo_pp_weight->Scale(1.0/pt_histo_pp_weight->Integral());

  pt_histo_PbPb->Scale(1.0/pt_histo_PbPb->Integral());
  //pt_histo_PbPb_weight->Scale(1.0/pt_histo_PbPb_weight->Integral());

  TCanvas c;
  histo_pp->SetLineColorAlpha(kBlue, 0.35);
  histo_PbPb->SetLineColorAlpha(kRed, 0.35);
  histo_PbPb->Draw("H");
  histo_pp->Draw("H same");

  TLegend *leg2 = new TLegend (0.7, 0.6, 0.9, 0.7);
  leg2->AddEntry(histo_pp, "pp", "l");
  leg2->AddEntry(histo_PbPb, "PbPb", "l");
  leg2->Draw("same");
  c.SetLogy();
  c.SaveAs("./results/Bu/pthat_pthatweight.gif");

}

