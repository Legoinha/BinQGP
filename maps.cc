#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

void maps(){

  TFile* f = new TFile("/eos/cms/store/group/phys_heavyions/zshi/SamplesForMariaSemiFinal/BPMC.root");
  TTree* t = (TTree*)f->Get("ntKp");

  float bpt, by, bmu1eta, bmu2eta, bmumueta;

  t->SetBranchAddress("Bpt", &bpt);
  t->SetBranchAddress("By", &by);
  t->SetBranchAddress("Bmu1eta", &bmu1eta);
  t->SetBranchAddress("Bmu2eta", &bmu2eta);
  t->SetBranchAddress("Bmumueta", &bmumueta);

  TH2F* histo_y = new TH2F("histo_y","histo_y",200,0,100,200,-2.4,2.4);
  TH2F* histo_eta1 = new TH2F("histo_eta1","histo_eta1",200,0,100,200,-2.4,2.4);
  TH2F* histo_eta2 = new TH2F("histo_eta2","histo_eta2",200,0,100,200,-2.4,2.4);
  TH2F* histo_eta3 = new TH2F("histo_eta3","histo_eta3",200,0,100,200,-2.4,2.4);

  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);

    histo_y->Fill(bpt,by);
    histo_eta1->Fill(bpt,bmu1eta);
    histo_eta2->Fill(bpt,bmu2eta);
    histo_eta3->Fill(bpt,bmumueta);

  }

  TCanvas c;
  c.cd();
  c.SetLogx();

  gStyle->SetOptStat(0);
  histo_y->SetTitle("");
  histo_y->GetXaxis()->SetTitle("p_{T} (B^{+}) [GeV/c]");
  histo_y->GetYaxis()->SetTitle("y (B^{+})");
  histo_y->Draw("COLZ");

  c.SaveAs("~/public/BinQGP/results/Bu/2d_map/by_map.gif");

  TCanvas c1;
  c1.cd();
  c1.SetLogx();

  gStyle->SetOptStat(0);
  histo_eta1->SetTitle("");
  histo_eta1->GetXaxis()->SetTitle("p_{T} (B^{+}) [GeV/c]");
  histo_eta1->GetYaxis()->SetTitle("#eta (#mu_{1})");
  histo_eta1->Draw("COLZ");

  c1.SaveAs("~/public/BinQGP/results/Bu/2d_map/beta1_map.gif");

  TCanvas c2;
  c2.cd();
  c2.SetLogx();

  gStyle->SetOptStat(0);
  histo_eta2->SetTitle("");
  histo_eta2->GetXaxis()->SetTitle("p_{T} (B^{+}) [GeV/c]");
  histo_eta2->GetYaxis()->SetTitle("#eta (#mu_{2})");
  histo_eta2->Draw("COLZ");

  c2.SaveAs("~/public/BinQGP/results/Bu/2d_map/beta2_map.gif");

  TCanvas c3;
  c3.cd();
  c3.SetLogx();

  gStyle->SetOptStat(0);
  histo_eta3->SetTitle("");
  histo_eta3->GetXaxis()->SetTitle("p_{T} (B^{+}) [GeV/c]");
  histo_eta3->GetYaxis()->SetTitle("#eta (#mu^{+} #mu^{-})");
  histo_eta3->Draw("COLZ");

  c3.SaveAs("~/public/BinQGP/results/Bu/2d_map/bmumueta_map.gif");

}
