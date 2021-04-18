/////////////////////////////////////////////////////////////////////////
//
// Bs and Bu mesons
//
// -Sideband subtraction and SPlot methods
// -MC comparisons
// -Data fit, fit validation and fit systematics
//
//
// August 2019
//
/////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <sstream>
#include <vector>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBifurGauss.h>
#include "TMath.h"
#include <RooGenericPdf.h>
#include "TRatioPlot.h"
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include <RooHist.h>
#include "RooStats/SPlot.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include <iostream>
#include <TF1.h>
#include <RooPolynomial.h>
#include <fstream>
#include <TGraph.h>
#include "TMultiGraph.h"
#include <TEfficiency.h>


using namespace RooStats;
using namespace RooFit;
using namespace std;

std::vector<TH1D*> sideband_subtraction(RooWorkspace w, int* n, int n_var);
std::vector<TH1D*> splot_method(RooWorkspace& w, int* n, TString* label, int n_var);

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram_mc(RooRealVar var, TTree* t, int n, TString weight); 
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n); 
void AddWeights(TTree* t);
void read_data(RooWorkspace& w, TString f_input);
void build_pdf (RooWorkspace& w, std::string choice = "nominal");
void plot_complete_fit(RooWorkspace& w);
void do_splot(RooWorkspace& w);
TH1D* make_splot(RooWorkspace& w, int n, TString label);
void validate_fit(RooWorkspace* w);
void get_ratio( std::vector<TH1D*>,  std::vector<TH1D*>,  std::vector<TString>, TString);
void pT_analysis(RooWorkspace& w,int n, TString, TString);
double get_yield_syst(RooDataSet *dt, TString syst_src);
void fit_syst_error(TString);
void fit_syst_error_bin(TString, double a, double b);


//particle
// 0 = Bu
// 1 = Bs
// 2 = B0

#define particle 2

//weights
// 1 = calculates ratio between MC and sPlot 
// 0 = does not calculate weights 

#define weights 1

//background
// 0 = does MC validation with data
// 1 = does MC validation with background

#define background 0

//add_weights
// 1 = adds weights to tree
// 0 = does not add weights to tree

# define add_weights 0


void bmesons_new(){
  
  int n_var;
  TString input_file_data;
  if(particle == 0){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewFlatSamples/BPData.root";}
  else if(particle == 1){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewFlatSamples/BsData.root";}
  else if(particle == 2){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/TrkQualCut/BZData.root";}

  TString input_file_mc;
  if(particle == 0){input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewFlatSamples/BPMC.root";}
  else if(particle == 1){input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewFlatSamples/BsMC.root";}
  else if(particle == 2){input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewFlatSamples/BsMC.root";} // change (still no MC for B0)

  TString input_file_reweighted_mc;
  if(particle == 0){input_file_reweighted_mc = "./results/Bu/mc_validation_plots/weights/tree_with_weight.root";}
  else if(particle == 1){input_file_reweighted_mc = "./results/Bs/mc_validation_plots/weights/tree_with_weight.root";}
  else if(particle == 2){input_file_reweighted_mc = "./results/B0/mc_validation_plots/weights/tree_with_weight.root";}

  std::vector<TH1D*> histos_sideband_sub;
  std::vector<TH1D*> histos_mc;
  std::vector<TH1D*> histos_splot;

#if particle == 0
  int n_bins[] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40} ;

  TString variables[] = {"PVx", "PVy", "PVz", "PVxE", "PVyE", "PVzE", "BvtxX", "BvtxY", "BvtxZtoPVZ", "BSx", "BSy", "BSz", "BSxErr", "BSyErr", "BSzErr", "BSdxdz", "BSdydz", "BSdxdzErr", "BSdydzErr", "BSWidthX", "BSWidthXErr", "BSWidthY", "BSWidthYErr", "By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr","BsvpvDistance_2D","BsvpvDisErr_2D", "Bmumumass", "Bmu1eta","Bmu2eta", "Bmu1pt", "Bmu2pt","Bmu1dxyPV","Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV","Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50", "BDT_pt_50_100"} ;    

#elif particle == 1
  int n_bins[] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40} ;

  TString variables[] = {"PVx", "PVy", "PVz", "PVxE", "PVyE", "PVzE", "BvtxX", "BvtxY", "BvtxZtoPVZ", "BSx", "BSy", "BSz", "BSxErr", "BSyErr", "BSzErr", "BSdxdz", "BSdydz", "BSdxdzErr", "BSdydzErr", "BSWidthX", "BSWidthXErr", "BSWidthY", "BSWidthYErr", "By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "BsvpvDistance_2D", "BsvpvDisErr_2D", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV","Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "Btrk2Pt", "Btrk2Eta", "Btrk2PtErr", "Btrk2Dz1", "Btrk2DzError1", "Btrk2Dxy1", "Btrk2DxyError1", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50", "BDT_pt_50_100"} ;

#elif particle == 2
  int n_bins[] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};

  TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk2Pt", "Btrk1Eta", "Btrk2Eta", "Btrk1PtErr", "Btrk2PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "BsvpvDistance_2D", "BsvpvDisErr_2D", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk2Dz1", "Btrk1DzError1", "Btrk2DzError1", "Btrk1Dxy1", "Btrk2Dxy1", "Btrk1DxyError1", "Btrk2DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt"};

#endif

  int n_n_bins = sizeof(n_bins)/sizeof(n_bins[0]);
  int n_variables = sizeof(variables)/sizeof(variables[0]);

  if(n_n_bins != n_variables){
    std::cout << "Error: number of bins does not correspond to number of variables." << std::endl;
    return;}

  n_var = n_variables;
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,input_file_data);

  build_pdf(*ws, "nominal");
  plot_complete_fit(*ws);  
  
  //validate_fit(ws);
  pT_analysis(*ws,n_bins[0], "pT.root", input_file_data);
  return;

  // SIDEBAND SUBTRACTION (needs to be run after plot_complete_fit)
  histos_sideband_sub = sideband_subtraction(*ws, n_bins, n_var);
 
  // SPLOT (fixes parameters of the fit -> they need to be unfixed for pT analysis) 
  do_splot(*ws); 
  histos_splot = splot_method(*ws,n_bins,variables, n_var);
  
  // MONTE CARLO HISTOGRAMS
  TFile *fin_mc = new TFile(input_file_mc); //use this file to add the weights (to clone original tree) and make data-MC comparisons without weights
  //TFile *fin_mc = new TFile(input_file_reweighted_mc); //use this file to make data-MC comparisons with weights

  TTree* t1_mc = particle ? (TTree*)fin_mc->Get("ntphi") : (TTree*)fin_mc->Get("ntKp");
  std::vector<TString> names;

  for(int i=0; i<n_var; i++){
    //std::cout<< "Var names: "<< histos_sideband_sub[i]->GetName()<<std::endl;
    TString weight = "weight";
    histos_mc.push_back(create_histogram_mc((*ws->var(variables[i])), t1_mc, n_bins[i], weight));
    names.push_back(TString(variables[i]));
  }
  
  // RATIO BETWENE DATA (SPLOT) AND MC
  if (weights == 1){get_ratio(histos_splot, histos_mc,names,"weights.root");}
  
  // ADDS WEIGHTS TO MC TREE (use to reweight MC)
  if (add_weights == 1){AddWeights(t1_mc);}  
  

  //COMPARISONS//
  
  //Sideband Subtraction vs. Monte Carlo
  vector<TH1D*> mc_comp_ss(histos_mc);

  vector<TH1D*> ss_comp_mc(histos_sideband_sub);

  for(int i=0; i<n_var; i++){
    TCanvas c;
    mc_comp_ss[i]->SetXTitle(variables[i]);
    mc_comp_ss[i]->SetStats(0);
    ss_comp_mc[i]->SetStats(0);

    if(particle == 0){
        mc_comp_ss[i]->SetTitle("B^{+}");}
    else if (particle == 1){
        mc_comp_ss[i]->SetTitle("B^{0}_{s}");}

    //normalization
    mc_comp_ss[i]->Scale(1/mc_comp_ss[i]->Integral());
    ss_comp_mc[i]->Scale(1/ss_comp_mc[i]->Integral());
    mc_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_ss[i]->GetMinimum(),1.1*mc_comp_ss[i]->GetMaximum());
    mc_comp_ss[i]->Draw();
    ss_comp_mc[i]->Draw("same");

    //y axis : maximum and minimum
    if((mc_comp_ss[i]->GetMaximum() > ss_comp_mc[i]->GetMaximum())){
      mc_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_ss[i]->GetMinimum(), 1.1*mc_comp_ss[i]->GetMaximum());}
    else if((ss_comp_mc[i]->GetMaximum() > mc_comp_ss[i]->GetMaximum())){
      mc_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*ss_comp_mc[i]->GetMinimum(), 1.1*ss_comp_mc[i]->GetMaximum());}

    //--TRATIO--//
    auto rp = new TRatioPlot(ss_comp_mc[i] ,mc_comp_ss[i], "divsym");
    c.SetTicks(0, 1);
    rp->SetH1DrawOpt("E");
    rp->Draw("nogrid");
    rp->GetLowerRefYaxis()->SetTitle("Data(ss)/MC");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    c.Update();

    TLegend* leg;
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(ss_comp_mc[i]->GetName(), "S. Subtraction", "l");
    leg->AddEntry(mc_comp_ss[i]->GetName(), "Monte Carlo", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
   
    if(particle == 0){
      c.SaveAs("./results/Bu/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bu.pdf");
      c.SaveAs("./results/Bu/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bu.gif");} 
    else if(particle == 1){
      c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.pdf");
      c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.gif");}
    //leg->Delete();
    }
  

  //SPlot vs. Sideband subtraction
  vector<TH1D*> sp_comp_ss(histos_splot);
  vector<TH1D*> ss_comp_sp(histos_sideband_sub);

  for(int i=0; i<n_var; i++){
    TCanvas a;
    ss_comp_sp[i]->SetYTitle("normalized entries");
    sp_comp_ss[i]->SetXTitle(variables[i]);
    ss_comp_sp[i]->SetStats(0);
    sp_comp_ss[i]->SetStats(0);

    if(particle == 0){
      sp_comp_ss[i]->SetTitle("B^{+}");}
    else if (particle == 1){
      sp_comp_ss[i]->SetTitle("B^{0}_{s}");}

    //normalization
    ss_comp_sp[i]->Scale(1/ss_comp_sp[i]->Integral());
    sp_comp_ss[i]->Scale(1/sp_comp_ss[i]->Integral());
    ss_comp_sp[i]->GetYaxis()->SetRangeUser(0.1*ss_comp_sp[i]->GetMinimum(),1.1*ss_comp_sp[i]->GetMaximum());
    ss_comp_sp[i]->Draw();
    sp_comp_ss[i]->Draw("same");

    //y axis: maximum and minimum
    if((ss_comp_sp[i]->GetMaximum() > sp_comp_ss[i]->GetMaximum())){
      sp_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*ss_comp_sp[i]->GetMinimum(), 1.1*ss_comp_sp[i]->GetMaximum());}
    else if((sp_comp_ss[i]->GetMaximum() > ss_comp_sp[i]->GetMaximum())){
      sp_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*sp_comp_ss[i]->GetMinimum(), 1.1*sp_comp_ss[i]->GetMaximum());}
 
    //--TRATIO--//
    auto rp = new TRatioPlot(ss_comp_sp[i], sp_comp_ss[i], "divsym");
    a.SetTicks(0, 1);
    rp->SetH1DrawOpt("E");
    rp->Draw("nogrid");
    rp->GetLowerRefYaxis()->SetTitle("Data(ss)/Data(sp)");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    a.Update();
 
    TLegend* leg;
    leg = new TLegend(0.6, 0.9, 0.9, 1.0);
    leg->AddEntry(ss_comp_sp[i]->GetName(), "Sideband Subtraction", "l");
    leg->AddEntry(sp_comp_ss[i]->GetName(), "SPlot", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
 
    if(particle == 0){
      a.SaveAs("./results/Bu/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bu.pdf");
      a.SaveAs("./results/Bu/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bu.gif");}
    else if(particle == 1){
      a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.pdf");
      a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.gif");}
    //leg->Delete();
  }
       
  
  //SPlot vs. Monte Carlo
  vector<TH1D*> sp_comp_mc(histos_splot);
  vector<TH1D*> mc_comp_sp(histos_mc);

  for(int i=0; i<n_var; i++){
    TCanvas a;
    mc_comp_sp[i]->SetXTitle(variables[i]);
    mc_comp_sp[i]->SetYTitle("normalized entries");
    mc_comp_sp[i]->SetStats(0);
    sp_comp_mc[i]->SetStats(0);

    if(particle == 0){
      sp_comp_mc[i]->SetTitle("B^{+}");}
    else if (particle == 1){
      sp_comp_mc[i]->SetTitle("B^{0}_{s}");}
	
    //normalization
    mc_comp_sp[i]->Scale(1/mc_comp_sp[i]->Integral());
    sp_comp_mc[i]->Scale(1/sp_comp_mc[i]->Integral());
    mc_comp_sp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_sp[i]->GetMinimum(),1.1*mc_comp_sp[i]->GetMaximum());
    mc_comp_sp[i]->Draw();
    sp_comp_mc[i]->Draw("same");

    //y axis: maximum and minimum 
    if((mc_comp_sp[i]->GetMaximum() > sp_comp_mc[i]->GetMaximum())){
      sp_comp_mc[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_sp[i]->GetMinimum(), 1.1*mc_comp_sp[i]->GetMaximum());}
    else if((sp_comp_mc[i]->GetMaximum() > mc_comp_sp[i]->GetMaximum())){
      sp_comp_mc[i]->GetYaxis()->SetRangeUser(0.1*sp_comp_mc[i]->GetMinimum(), 1.1*sp_comp_mc[i]->GetMaximum());}
	
    //--TRATIO--//	
    auto rp = new TRatioPlot(sp_comp_mc[i], mc_comp_sp[i], "divsym");
    a.SetTicks(0, 1);
    rp->SetH1DrawOpt("E");
    rp->Draw("nogrid");
    rp->GetLowerRefYaxis()->SetTitle("Data(sp)/MC");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    a.Update();
     
    TLegend* leg;	
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(mc_comp_sp[i]->GetName(), "Monte Carlo", "l");
    leg->AddEntry(sp_comp_mc[i]->GetName(), "SPlot", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
	
    if(particle == 0){
      a.SaveAs("./results/Bu/mc_validation_plots/mc_sp/" + names[i]+"_mc_validation_Bu.pdf");
      a.SaveAs("./results/Bu/mc_validation_plots/mc_sp/" + names[i]+"_mc_validation_Bu.gif");}
    else if(particle == 1){
      a.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
      a.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.gif");}
	
    //leg->Delete();
  }
 
  //Sideband subtraction vs. Monte Carlo vs SPlot
  vector<TH1D*> sp_comp(histos_splot);
  vector<TH1D*> mc_comp(histos_mc);
  vector<TH1D*> ss_comp(histos_sideband_sub);

  for(int i=0; i<n_var; i++){
    TCanvas a;	
    mc_comp[i]->SetXTitle(variables[i]);
    mc_comp[i]->SetYTitle("normalized entries");
    mc_comp[i]->SetStats(0);
    sp_comp[i]->SetStats(0);
    ss_comp[i]->SetStats(0);

    if(particle == 0){
      mc_comp[i]->SetTitle("B^{+}");}
    else if (particle == 1){
      mc_comp[i]->SetTitle("B^{0}_{s}");}
	
    //normalization
    mc_comp[i]->Scale(1/mc_comp[i]->Integral());
    sp_comp[i]->Scale(1/sp_comp[i]->Integral());
    ss_comp[i]->Scale(1/ss_comp[i]->Integral());	
	
    //y axis: maximum and minimum 
    if((mc_comp[i]->GetMaximum() > sp_comp[i]->GetMaximum()) && (mc_comp[i]->GetMaximum() > ss_comp[i]->GetMaximum())){
      mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*mc_comp[i]->GetMaximum());}
    else if((sp_comp[i]->GetMaximum() > ss_comp[i]->GetMaximum()) && (sp_comp[i]->GetMaximum() > mc_comp[i]->GetMaximum())){
      mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*sp_comp[i]->GetMaximum());}
    else{mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*ss_comp[i]->GetMaximum());}
	
    mc_comp[i]->Draw();
    sp_comp[i]->Draw("same");
    ss_comp[i]->Draw("same");
		
    TLegend* leg;
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);	
    leg->AddEntry(ss_comp[i]->GetName(), "S. Subtraction", "l");
    leg->AddEntry(mc_comp[i]->GetName(), "Monte Carlo", "l");
    leg->AddEntry(sp_comp[i]->GetName(), "SPlot", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
	
    if(particle == 0){
      a.SaveAs("./results/Bu/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bu.pdf");
      a.SaveAs("./results/Bu/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bu.gif");}
    else if(particle == 1){
      a.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
      a.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.gif");}

    //leg->Delete();
  }
//comparisons end
}
//main function ends


void pT_analysis(RooWorkspace& w, int n, TString ptfile, TString datafile){

  TString dir_name;
  if(particle == 0){dir_name = "./results/Bu/Bpt/";}
  else if(particle == 1){dir_name = "./results/Bs/Bpt/";}
  else if(particle == 2){dir_name = "./results/B0/Bpt/";}

  TFile* f_wei = new TFile(dir_name + "/"+ ptfile, "recreate"); 

  RooAbsPdf*  model = w.pdf("model");
  RooRealVar* Bpt  = w.var("Bpt");
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooRealVar Bmass = *(w.var("Bmass"));

#if particle == 0
  const int n_pt_bins = 7;
  double pt_bins [n_pt_bins + 1] = {3,5,7,10,15,20,50,100};  
#elif particle == 1
  const int n_pt_bins = 7;
  double pt_bins[n_pt_bins + 1] = {3,5,7,10,15,20,50,100};
#elif particle == 2
  const int n_pt_bins = 7;
  double pt_bins[n_pt_bins + 1] = {3,5,7,10,15,20,50,100};
#endif  

  double pt_mean[n_pt_bins];
  double pt_low[n_pt_bins];
  double pt_high[n_pt_bins];

  double yield[n_pt_bins];
  double yield_err_low[n_pt_bins];
  double yield_err_high[n_pt_bins];
  double yield_err_syst[n_pt_bins];
  double m_yield_err_syst[n_pt_bins];

  for(int k = 1; k<n_pt_bins; k++){
    yield_err_syst[k]=0;
    m_yield_err_syst[k]=0;
  } 

  RooDataSet* data_pt, data_w, data_wp;
  RooFitResult* fit_pt;
  RooRealVar* n_sig_pt;
  RooRealVar* n_comb_pt;
 
  const int n_pdf_syst=4;
  //number of pdfs
  TString syst_src[n_pdf_syst]={"nominal","signal1gauss","triple_gauss","crystal_ball"};//bkg_range, bkg_poly (change table label)
  //array of pdfs to evaluate the systematic uncertainty of the N signal
  double yield_syst[n_pt_bins][n_pdf_syst]; 
  double yield_syst_rel[n_pt_bins][n_pdf_syst];
  //value of the systematic uncertainty

  for(int i=0;i<n_pt_bins;i++){
    //select data subset corresponding to pT bin
    data_pt = (RooDataSet*) data->reduce(Form("Bpt>%lf",pt_bins[i]));
    data_pt = (RooDataSet*) data_pt->reduce(Form("Bpt<%lf",pt_bins[i+1]));
    w.import(*data_pt, Rename(Form("data_pt_%d",i)));
   
    //perform fit and save result
    fit_pt = model->fitTo(*data_pt, Minos(true), Save());

    //plots the fit result
    TCanvas b;
    b.SetTitle("");

    TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
    p1->SetTitle("");
    p1->SetBorderMode(1);
    p1->SetFrameBorderMode(0);
    p1->SetBorderSize(2);
    p1->SetBottomMargin(0.10);
    p1->Draw();

    TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
    p2->SetTitle("");
    p2->SetTopMargin(0.);
    p2->SetBottomMargin(0.2);
    p2->SetBorderMode(1);
    p2->Draw();

    p1->cd();
    
    RooPlot* massframe = Bmass.frame(Title(""));
    data_pt->plotOn(massframe,RooFit::Name("Data"));
    model->paramOn(massframe,Layout(0.60,0.99,0.95));

    model->plotOn(massframe, RooFit::Name("Fit"),Range("all"));

    //TCanvas b;
    massframe->Draw();
    
    RooHist* pull_hist = massframe->pullHist("Data","Fit");
    RooPlot* pull_plot = Bmass.frame(Title(""));

    pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
    pull_plot->SetTitle("");

    pull_plot->GetXaxis()->SetTitle("");
    pull_plot->GetXaxis()->SetTitleFont(42);
    pull_plot->GetXaxis()->SetTitleSize(0.17);
    pull_plot->GetXaxis()->SetTitleOffset(1.09);
    pull_plot->GetXaxis()->SetLabelFont(42);
    pull_plot->GetXaxis()->SetLabelSize(0.15);
    pull_plot->GetXaxis()->SetLabelOffset(0.01);
    pull_plot->GetXaxis()->SetTickLength(0.13);

    pull_plot->GetYaxis()->SetTitle("Pull hist");
    pull_plot->GetYaxis()->SetTitleFont(42);
    pull_plot->GetYaxis()->SetTitleSize(0.10);
    pull_plot->GetYaxis()->SetTitleOffset(1.09);
    pull_plot->GetYaxis()->SetLabelFont(42);
    pull_plot->GetYaxis()->SetLabelSize(0.13);
    pull_plot->GetYaxis()->SetLabelOffset(0.005);
    pull_plot->GetYaxis()->SetNdivisions(305);

    p2->cd();
    pull_plot->Draw();
    
    if(particle == 0){
      b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf_fit_plot_Bu.gif", pt_bins[i], pt_bins[i+1]));
      //b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf_fit_plot_Bu.pdf", pt_bins[i], pt_bins[i+1]));
    }
    else if(particle == 1){
      b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf_fit_plot_Bs.gif", pt_bins[i], pt_bins[i+1]));
      //b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf_fit_plot_Bs.pdf", pt_bins[i], pt_bins[i+1]));
    }
    else if(particle == 2){
      b.SaveAs(Form("./results/B0/Bpt/%lf_%lf_fit_plot_B0.gif", pt_bins[i], pt_bins[i+1]));
      //b.SaveAs(Form("./results/B0/Bpt/%lf_%lf_fit_plot_B0.pdf", pt_bins[i], pt_bins[i+1]));
    }

    //YIELD + STATISTICAL ERROR
    //floatParsFinal returns the list of floating parameters after fit
    //cout << "Value of floating parameters" << endl;
    //fit_pt->floatParsFinal().Print("s");

    n_sig_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_signal");
    n_comb_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_combinatorial");

    yield[i]= n_sig_pt->getVal();
    yield_err_low [i] = n_sig_pt->getError(); 
    yield_err_high[i] = n_sig_pt->getError(); 

    //cout << "test asym error:" << n_sig_pt->getErrorLo() << " " <<  n_sig_pt->getAsymErrorLo() << " symmetric: " <<  n_sig_pt->getError() <<  endl;

    //SYSTEMATICS
    for(int k = 1; k<n_pdf_syst; k++){
      double val = 0.; 
      double val_nominal = 0.;
      val = get_yield_syst(data_pt, syst_src[k]); // gets yield value for bin i using pdf choice k
      val_nominal = get_yield_syst(data_pt, syst_src[0]); //gets yield value for bin i using nominal pdf choice
      cout<<"syst nominal: "<<syst_src[0]<<endl;
      yield_syst_rel[i][k] = (val - val_nominal)/val_nominal;
      yield_syst[i][k] = (val - val_nominal);
    }

    yield_err_syst[i] = pow(yield_syst[i][0],2) + pow(yield_syst[i][1],2) + pow(yield_syst[i][2],2) + pow(yield_syst[i][3],2);
    m_yield_err_syst[i] = sqrt(yield_err_syst[i]);
  }

    //separate sPlot from the rest, because sPlot fixes parameters of the fit
    
    //SPLOT	
    for(int i=0;i<n_pt_bins;i++){
    data_pt = (RooDataSet*) data->reduce(Form("Bpt>%lf",pt_bins[i]));
    data_pt = (RooDataSet*) data_pt->reduce(Form("Bpt<%lf",pt_bins[i+1]));

    //sPlot technique requires model parameters (other than the yields) to be fixed
    //dosp = true -> the splot technique is applied
    bool dosp = true;
    if (dosp){
      RooRealVar* mean  = w.var("mean");
      RooRealVar* sigma1 = w.var("sigma1");
      RooRealVar* sigma2 = w.var("sigma2");
      RooRealVar* cofs = w.var("cofs");
      RooRealVar* lambda = w.var("lambda");
      
      mean->setConstant();
      sigma1->setConstant();
      sigma2->setConstant();
      cofs->setConstant();
      lambda->setConstant();
      
      SPlot("sData","An sPlot",*data_pt, model, RooArgList(*n_sig_pt,*n_comb_pt));
      
      w.import(*data_pt, Rename(Form("data_pt_WithSWeights_%d",i)));
      
      RooDataSet* data_w = (RooDataSet*) w.data(Form("data_pt_WithSWeights_%d",i));

      RooDataSet* data_wb = new RooDataSet(data_w->GetName(),data_w->GetTitle(),data_w,*data_w->get(),0,"n_signal_sw"); 
     
      //weighted average pT
      double mean_w=data_wb->mean(*Bpt);
      double mean_s=data_pt->mean(*Bpt);
      pt_mean[i] = data_wb->mean(*Bpt);
      cout<<"mean_weight:"<<mean_w<<endl;
      cout<<"mean:"<< mean_s<<endl;
      
      pt_low[i]= pt_mean[i]-pt_bins[i];
      pt_high[i]= pt_bins[i+1]-pt_mean[i];
    }else{
      pt_low [i]= 0.5*(pt_bins[i+1]-pt_bins[i]);
      pt_high[i]= 0.5*(pt_bins[i+1]-pt_bins[i]);  
    }

    //normalize yield to bin width
    double bin_width = pt_bins[i+1]-pt_bins[i];
    yield[i] = yield[i]/bin_width;
    yield_err_low[i] = yield_err_low[i]/bin_width;
    yield_err_high[i] = yield_err_high[i]/bin_width;
    m_yield_err_syst[i] = m_yield_err_syst[i]/bin_width;
  }

  //plot yield vs average pT
  TCanvas c;
  TMultiGraph* mg = new TMultiGraph();

  TGraphAsymmErrors* gr = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
  gr->SetTitle("");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(1);
  gr->SetLineColor(1);
  gr->SetTitle("Differential Signal Yield");
  gr->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  gr->GetYaxis()->SetTitle("dN_{s}/dp_{T} [GeV^{-1}]");
  //gr->Draw("AP");
  gr->Write();
  //f_wei->Close();
  //delete f_wei;

  double pt_zero[n_pt_bins];
  for (int i=0;i<n_pt_bins;i++) pt_zero[i]= 0.;

  TGraphAsymmErrors* grs = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_zero,pt_zero, m_yield_err_syst, m_yield_err_syst);
  grs->SetTitle("");
  grs->SetMarkerColor(4);
  grs->SetMarkerStyle(1);
  grs->SetLineColor(2);
  grs->SetTitle("Differential Signal Yield");
  grs->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  grs->GetYaxis()->SetTitle("dN_{s}/dp_{T} [GeV^{-1}]");
  grs->Write();
  f_wei->Close();
  //grs->Draw("same");

  mg->Add(gr);
  mg->Add(grs);
   
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(gr, "Statistical Uncertainty", "lp");
  leg->AddEntry(grs, "Systematic Uncertainty", "lp");

  mg->Draw("AP");
  leg->Draw();
  mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  mg->GetYaxis()->SetTitle("dN_{s}/dp_{T} [GeV^{-1}]");
     
  if(particle == 0){
    //c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.pdf");
    c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.gif");
  }
  else if(particle == 1){
    //c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.pdf");
    c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.gif");
  }
  else if(particle == 2){
    //c.SaveAs("./results/B0/Bpt/raw_yield_pt_B0.pdf");
    c.SaveAs("./results/B0/Bpt/raw_yield_pt_B0.gif");
  }
  //save both vectors (systematic and statistical errors) in a text file

  TCanvas l;
  //log scale
  l.SetLogy();
  TGraphAsymmErrors* grlog = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
  grlog->SetTitle("");
  grlog->SetMarkerColor(4);
  grlog->SetMarkerStyle(21);
  grlog->SetTitle("Differential Signal Yield (Logscale)");
  grlog->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  grlog->GetYaxis()->SetTitle("dN_{s}/dp_{T} [GeV^{-1}]");
  grlog->Draw("AP");

  if(particle == 0){
    //l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.pdf");
    l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.gif");
  }
  else if(particle == 1){
    //l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.pdf");
    l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.gif");
  }
  else if(particle == 2){
    //l.SaveAs("./results/B0/Bpt/raw_yield_pt_logscale_B0.pdf");
    l.SaveAs("./results/B0/Bpt/raw_yield_pt_logscale_B0.gif");
  }
  //} 
  
  //evaluates the N systematics for each bin
  //for(int i=0;i<n_pt_bins;i++){
  // fit_syst_error_bin(datafile, pt_bins[i], pt_bins[i+1]);
  // }

   cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-50" << '|' << setw(15) << "50-100" << '|' << endl;
   cout << '|' << setw(15) << "Nominal" << '|' << setw(15) << yield_syst_rel[0][0] << '|' << setw(15) << yield_syst_rel[1][0] << '|' << setw(15) << yield_syst_rel[2][0] << '|' << setw(15) << yield_syst_rel[3][0] << '|' << setw(15) << yield_syst_rel[4][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << setw(15) << yield_syst_rel[6][0]  << '|' << endl;
   cout << '|' << setw(15) << "Signal1gauss" << '|' << setw(15) << yield_syst_rel[0][1] << '|' << setw(15) << yield_syst_rel[1][1] << '|' << setw(15) << yield_syst_rel[2][1] << '|' << setw(15) << yield_syst_rel[3][1] << '|' << setw(15) << yield_syst_rel[4][1] << '|' << setw(15) << yield_syst_rel[5][1]  << '|' << setw(15) << yield_syst_rel[6][1]<< '|' << endl;
   cout << '|' << setw(15) << "Triple_gauss" << '|' << setw(15) << yield_syst_rel[0][2] << '|' << setw(15) << yield_syst_rel[1][2] << '|' << setw(15) << yield_syst_rel[2][2] << '|' << setw(15) << yield_syst_rel[3][2] << '|' << setw(15) << yield_syst_rel[4][2] << '|' << setw(15) << yield_syst_rel[5][2] << '|' << setw(15) << yield_syst_rel[6][2] << '|' << endl;
   cout << '|' << setw(15) << "Crystal_ball" << '|' << setw(15) << yield_syst_rel[0][3] << '|' << setw(15) <<  yield_syst_rel[1][3] << '|' << setw(15) <<  yield_syst_rel[2][3] << '|' << setw(15) <<  yield_syst_rel[3][3] << '|' << setw(15) << yield_syst_rel[4][3] << '|' << setw(15) << yield_syst_rel[5][3] << '|' << setw(15) << yield_syst_rel[6][3] << '|' << endl;
}
//pT_analysis end

//get the ratio between the data (splot method) and the MC and save it in a root file
void get_ratio( std::vector<TH1D*> data, std::vector<TH1D*> mc,  std::vector<TString> v_name, TString filename) {

  TString dir_name = particle ? "./results/Bs/mc_validation_plots/weights/" : "./results/Bu/mc_validation_plots/weights/";

  TFile* f_wei = new TFile(dir_name + "/"+ filename, "recreate");


  TH1D* h_aux;
  //std::vector<TH1D*> histos;
  h_aux->SetDefaultSumw2(kTRUE);


  for(int i=0; i<(int)data.size(); i++) {

    //auto rp = new TRatioPlot(histos_splot[i], histos_mc[i], "divsym");
    //rp->Write("ratioplot_"+variables[i]);

    h_aux = (TH1D*)data.at(i)->Clone("weights_"+v_name.at(i));

    h_aux->SetMaximum(6.);
    h_aux->SetMinimum(0.);
    h_aux->SetStats(0);
    h_aux->GetYaxis()->SetTitle("Data / MC");

    //normalization
    h_aux->Scale(1/h_aux->Integral());
    mc[i]->Scale(1/mc[i]->Integral());

    h_aux->Divide(mc.at(i));
    
    f_wei->cd();
    h_aux->Write();
    
    TCanvas c;
    h_aux->Draw();
    c.SaveAs(dir_name+"/"+v_name.at(i) + "_weights.gif");
    //output: a root file and plots gifs
  
  }

  f_wei->Write();
  f_wei->ls();
  f_wei->Close();


  return;
}

//get_ratio ends

void read_data(RooWorkspace& w, TString f_input){
  TFile* fin_data = new TFile(f_input);
  TTree* t1_data;

  if(particle == 0){t1_data = (TTree*)fin_data->Get("ntKp");}
  else if(particle == 1){t1_data = (TTree*)fin_data->Get("ntphi");}
  else if(particle == 2){t1_data = (TTree*)fin_data->Get("ntKstar");}

  RooArgList arg_list ("arg_list");

  arg_list.add(*(w.var("Bmass")));
  if(particle != 2){
    arg_list.add(*(w.var("PVx")));
    arg_list.add(*(w.var("PVy")));
    arg_list.add(*(w.var("PVz")));
    arg_list.add(*(w.var("PVxE")));
    arg_list.add(*(w.var("PVyE")));
    arg_list.add(*(w.var("PVzE")));
    arg_list.add(*(w.var("BvtxX")));
    arg_list.add(*(w.var("BvtxY")));
    arg_list.add(*(w.var("BvtxZtoPVZ")));
    arg_list.add(*(w.var("BSx")));
    arg_list.add(*(w.var("BSy")));
    arg_list.add(*(w.var("BSz")));
    arg_list.add(*(w.var("BSxErr")));
    arg_list.add(*(w.var("BSyErr")));
    arg_list.add(*(w.var("BSzErr")));
    arg_list.add(*(w.var("BSdxdz")));
    arg_list.add(*(w.var("BSdydz")));
    arg_list.add(*(w.var("BSdxdzErr")));
    arg_list.add(*(w.var("BSdydzErr")));
    arg_list.add(*(w.var("BSWidthX")));
    arg_list.add(*(w.var("BSWidthXErr")));
    arg_list.add(*(w.var("BSWidthY")));
    arg_list.add(*(w.var("BSWidthYErr")));
  }
  arg_list.add(*(w.var("By")));
  arg_list.add(*(w.var("Bpt")));
  arg_list.add(*(w.var("Btrk1Pt")));
  arg_list.add(*(w.var("Btrk1Eta")));
  arg_list.add(*(w.var("Btrk1PtErr")));
  arg_list.add(*(w.var("Bchi2cl")));
  arg_list.add(*(w.var("BsvpvDistance")));
  arg_list.add(*(w.var("BsvpvDisErr")));
  arg_list.add(*(w.var("BsvpvDistance_2D")));
  arg_list.add(*(w.var("BsvpvDisErr_2D")));
  arg_list.add(*(w.var("Bmumumass")));
  arg_list.add(*(w.var("Bmu1eta")));
  arg_list.add(*(w.var("Bmu2eta")));
  arg_list.add(*(w.var("Bmu1pt")));
  arg_list.add(*(w.var("Bmu2pt")));
  arg_list.add(*(w.var("Bmu1dxyPV")));
  arg_list.add(*(w.var("Bmu2dxyPV")));
  arg_list.add(*(w.var("Bmu1dzPV")));
  arg_list.add(*(w.var("Bmu2dzPV")));
  arg_list.add(*(w.var("Bd0")));
  arg_list.add(*(w.var("Bd0Err")));
  arg_list.add(*(w.var("Bdtheta")));
  arg_list.add(*(w.var("Balpha")));
  arg_list.add(*(w.var("Btrk1Dz1")));
  arg_list.add(*(w.var("Btrk1DzError1")));
  arg_list.add(*(w.var("Btrk1Dxy1")));
  arg_list.add(*(w.var("Btrk1DxyError1")));
  arg_list.add(*(w.var("Bmumueta")));
  arg_list.add(*(w.var("Bmumuphi")));
  arg_list.add(*(w.var("Bmumupt")));
  if(particle == 0){
    arg_list.add(*(w.var("BDT_pt_3_5")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
    arg_list.add(*(w.var("BDT_pt_50_100")));
  }   
  if( (particle == 1) || (particle == 2)){
    arg_list.add(*(w.var("Btrk2Pt")));
    arg_list.add(*(w.var("Btrk2Eta")));
    arg_list.add(*(w.var("Btrk2PtErr")));
    arg_list.add(*(w.var("Btrk2Dz1")));
    arg_list.add(*(w.var("Btrk2DzError1")));
    arg_list.add(*(w.var("Btrk2Dxy1")));
    arg_list.add(*(w.var("Btrk2DxyError1")));
  }
  if(particle == 1){
    arg_list.add(*(w.var("BDT_pt_3_5")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
    arg_list.add(*(w.var("BDT_pt_50_100")));
  }
  
  RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
  w.import(*data, Rename("data"));
}

void build_pdf(RooWorkspace& w, std::string choice) {

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*) w.data("data");

  //  RooDataSet* reduceddata_central;
  double left;
  double right;
  double mass_peak;

  if(particle == 0){
    left = 5.2;
    right = 5.4;
    mass_peak = 5.265;
  }
  else if(particle == 1){
    left = 5.3;
    right = 5.45;
    mass_peak = 5.366;
  }
  else if(particle == 2){
    left = 5.2;
    right = 5.4;
    mass_peak = 5.280;
  }

  //SIGNAL// (for Bs and Bu signal + RT B0 signal)
  // sum of two gaussians
  RooRealVar mean("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.1);
  RooRealVar sigma1("sigma1","sigma1",0.02,0.005,0.05);
  RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);
  RooRealVar sigma2("sigma2","sigma2",0.01,0.005,0.05);
  RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);
  RooRealVar cofs("cofs", "cofs", 0.3, 0., 1.);
  RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);

  // triple gaussian
  RooRealVar sigma3("sigma3","sigma3",0.012,0.010,0.030);
  RooGaussian signal3("signal3","signal3",Bmass, mean, sigma3);  
  RooRealVar cofs1("cofs1", "cofs1", 0.3, 0., 1.);
  RooAddPdf signal_triple("signal_triple", "signal_triple", RooArgList(signal1,signal2,signal3), RooArgList(cofs,cofs1));

  // crystal ball funciton
  RooRealVar alpha("alpha", "alpha", mass_peak-0.015/2, mass_peak-0.08, mass_peak-0.003);
  RooRealVar n("n", "n", 5.00);
  RooCBShape crystal_ball("crystal_ball", "crystal_ball", Bmass, mean, sigma1, alpha, n);


  //BACKGROUND//
  //error function (for Bu peaking background)
  RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
  RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  
  m_nonprompt_shift.setConstant(kTRUE);
  m_nonprompt_scale.setConstant(kTRUE);
  RooGenericPdf erf("erf","erf","TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Bmass,m_nonprompt_scale,m_nonprompt_shift));
 
  //exponential (for Bs, Bu and B0 combinatorial background)
  RooRealVar lambda("lambda","lambda",-2.,-5.,1.0);
  RooExponential fit_side("fit_side", "fit_side", Bmass, lambda);

  // 1st order polynomial (for Bs, Bu and B0 combinatorial background)
  RooRealVar slope("slope","slope",0,-5,5);
  RooPolynomial poly_bkg("poly_bkg", "poly_bkg", Bmass, slope);

  //jpsi_pi component (for Bu jpsi background)
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693, 5.346, 5.347);
  RooRealVar m_jpsipi_sigma1l("m_jpsipi_sigma1l","m_jpsipi_sigma1l",0.0290762,0.010,0.150);
  RooRealVar m_jpsipi_sigma1r("m_jpsipi_sigma1r","m_jpsipi_sigma1r",0.0652519,0.010,0.150);
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876);
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073);
  RooRealVar m_jpsipi_sigma2("m_jpsipi_sigma2","m_jpsipi_sigma2",0.0994712,0.020,0.500);
  RooRealVar m_jpsipi_sigma3("m_jpsipi_sigma3","m_jpsipi_sigma3",0.330152,0.020,0.500);
  RooRealVar m_jpsipi_fraction2("m_jpsipi_fraction2","m_jpsipi_fraction2",0.234646,0.0,1.0);
  RooRealVar m_jpsipi_fraction3("m_jpsipi_fraction3","m_jpsipi_fraction3",0.114338,0.0,1.0);

  m_jpsipi_mean1.setConstant(kTRUE);
  m_jpsipi_mean2.setConstant(kTRUE);
  m_jpsipi_mean3.setConstant(kTRUE);
  m_jpsipi_sigma1l.setConstant(kTRUE);
  m_jpsipi_sigma1r.setConstant(kTRUE);
  m_jpsipi_sigma2.setConstant(kTRUE);
  m_jpsipi_sigma3.setConstant(kTRUE);
  m_jpsipi_fraction2.setConstant(kTRUE);
  m_jpsipi_fraction3.setConstant(kTRUE);  

  RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1","m_jpsipi_gaussian1",Bmass,m_jpsipi_mean1,m_jpsipi_sigma1l,m_jpsipi_sigma1r);
  RooGaussian m_jpsipi_gaussian2("m_jpsipi_gaussian2","m_jpsipi_gaussian2",Bmass,m_jpsipi_mean2,m_jpsipi_sigma2);
  RooGaussian m_jpsipi_gaussian3("m_jpsipi_gaussian3","m_jpsipi_gaussian3",Bmass,m_jpsipi_mean3,m_jpsipi_sigma3);
  RooAddPdf jpsipi("jpsipi","jpsipi",RooArgList(m_jpsipi_gaussian3,m_jpsipi_gaussian2,m_jpsipi_gaussian1),RooArgList(m_jpsipi_fraction3,m_jpsipi_fraction2));
  Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);
  Bmass.setRange("peakright",left,Bmass.getMax());

  // K Pi SWAP COMPONENT ( for WT B0 signal)
  RooRealVar sigma_swapped1("sigma_swapped1","sigma_swapped1", 0.1133, 0.010, 0.150);
  RooRealVar sigma_swapped2("sigma_swapped2","sigma_swapped2", 0.01529, 0.010, 0.150);
  RooRealVar sigma_swapped3("sigma_swapped3","sigma_swapped3", 0.0424, 0.010, 0.150);
  RooRealVar alpha1("alpha1","alpha1", 1.78, -20., 20.);
  RooRealVar alpha2("alpha2","alpha2", 0.150, -20., 20.);
  RooRealVar alpha3("alpha3","alpha3", -6.802, -20., 20.);
  RooRealVar n1_parameter("n1_parameter", "n1_parameter", 32., 0., 300.);
  RooRealVar n2_parameter("n2_parameter", "n2_parameter", 98., 0., 300.);
  RooRealVar n3_parameter("n3_parameter", "n3_parameter", 179., 0., 300.);
  RooRealVar r1("r1","r1", 0.249, 0.0, 1.0);
  RooRealVar r2("r2","r2", 0.3922, 0.0, 1.0);

  sigma_swapped1.setConstant(kTRUE);
  sigma_swapped2.setConstant(kTRUE);
  sigma_swapped3.setConstant(kTRUE);
  alpha1.setConstant(kTRUE);
  alpha2.setConstant(kTRUE);
  alpha3.setConstant(kTRUE);
  n1_parameter.setConstant(kTRUE);
  n2_parameter.setConstant(kTRUE);
  n3_parameter.setConstant(kTRUE);
  r1.setConstant(kTRUE);
  r2.setConstant(kTRUE);

  RooCBShape swapped1("swapped1","swapped1", Bmass, mean, sigma_swapped1, alpha1, n1_parameter);
  RooCBShape swapped2("swapped2","swapped2", Bmass, mean, sigma_swapped2, alpha2, n2_parameter);
  RooCBShape swapped3("swapped3","swapped3", Bmass, mean, sigma_swapped3, alpha3, n3_parameter);

  RooAddPdf k_pi_swap("k_pi_swap","k_pi_swap", RooArgSet(swapped1,swapped2,swapped3), RooArgSet(r1,r2));

  // B0 FULL MODEL
  RooRealVar f_swap("f_swap","f_swap", 0.5);
  RooAddPdf B0_signal("B0_signal","B0_signal",RooArgList(k_pi_swap,signal),RooArgList(f_swap));
  RooAddPdf B0_signal1gauss("B0_signal1gauss","B0_signal1gauss",RooArgList(k_pi_swap,signal3),RooArgList(f_swap));
  RooAddPdf B0_signaltriple("B0_signaltriple","B0_signaltriple",RooArgList(k_pi_swap,signal_triple),RooArgList(f_swap));
  RooAddPdf B0_signalCB("B0_signalCB","B0_signalCB",RooArgList(k_pi_swap,crystal_ball),RooArgList(f_swap));

  // NORMALISATIONS
  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
  RooRealVar f_erf("f_erf","f_erf",0.2,0.0,1.);//2.50259e-01,0.0,1.
  RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",0.03996, 0.038, 0.040); 
  f_jpsipi.setConstant(kTRUE);
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));

  if(particle == 0){//B+
    if (choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }else if (choice=="bkg_poly"){
      RooAddPdf model("model", "model", RooArgList(signal,poly_bkg,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    } else if (choice=="bkg_range"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side,jpsipi),RooArgList(n_signal,n_combinatorial, n_jpsipi));
      w.import(model);
    }else if (choice == "signal1gauss"){
      RooAddPdf model("model", "model", RooArgList(signal3,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }else if (choice == "triple_gauss"){
      RooAddPdf model("model", "model", RooArgList(signal_triple,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }else if (choice == "crystal_ball"){
      RooAddPdf model("model", "model", RooArgList(crystal_ball,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }
  }

  else if(particle == 1){//Bs
    if (choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side), RooArgList(n_signal, n_combinatorial)); 
      w.import(model);
    }else if (choice=="bkg_poly"){
      RooAddPdf model("model", "model", RooArgList(signal,poly_bkg),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    } else if (choice=="bkg_range"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "signal1gauss"){
      RooAddPdf model("model", "model", RooArgList(signal3,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "triple_gauss"){
      RooAddPdf model("model", "model", RooArgList(signal_triple,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "crystal_ball"){
      RooAddPdf model("model", "model", RooArgList(crystal_ball,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }
  }

  else if(particle == 2){
    if (choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(B0_signal,fit_side), RooArgList(n_signal, n_combinatorial));
      w.import(model);
    }else if (choice=="bkg_poly"){
      RooAddPdf model("model", "model", RooArgList(B0_signal,poly_bkg),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    } else if (choice=="bkg_range"){
      RooAddPdf model("model", "model", RooArgList(B0_signal,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "signal1gauss"){
      RooAddPdf model("model", "model", RooArgList(B0_signal1gauss,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "triple_gauss"){
      RooAddPdf model("model", "model", RooArgList(B0_signaltriple,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }else if (choice == "crystal_ball"){
      RooAddPdf model("model", "model", RooArgList(B0_signalCB,fit_side),RooArgList(n_signal,n_combinatorial));
      w.import(model);
    }
  }
}
//build_pdf ends

void fit_syst_error(TString fname){
  //returns the yield's systematic uncertainty
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,fname);


  RooDataSet* data = (RooDataSet*) ws->data("data");

  build_pdf(*ws,"nominal");
  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres_nom = model->fitTo(*data,Save());

  build_pdf(*ws,"bkg_poly");
  model = ws->pdf("model");
  RooFitResult* fitres_bgpol = model->fitTo(*data,Save());

  build_pdf(*ws,"bkg_range");
  model = ws->pdf("model");
  RooFitResult* fitres_bgrange = model->fitTo(*data,Range("peakright"),Save());

  build_pdf(*ws,"sig1gauss");
  model = ws->pdf("model");
  RooFitResult* fitres_sig1 = model->fitTo(*data,Save());
  
  RooRealVar* n_nom = (RooRealVar*) fitres_nom  ->floatParsFinal().find("n_signal");
  RooRealVar* n_bgp = (RooRealVar*) fitres_bgpol->floatParsFinal().find("n_signal");
  RooRealVar* n_bra = (RooRealVar*) fitres_bgrange->floatParsFinal().find("n_signal");
  RooRealVar* n_sig1 = (RooRealVar*) fitres_sig1->floatParsFinal().find("n_signal");

  double n0  = n_nom->getVal();
  double n1  = n_bgp->getVal();
  double n2  = n_bra->getVal();
  double n3 = n_sig1->getVal();

  double syst1 = (n1-n0)/n0;
  double syst2 = (n2-n0)/n0;
  double syst3 = (n3-n0)/n0;
  
  cout << "syst bg pdf:" << syst1 * 100 << "\%\trange" << syst2 * 100 << "\%\tsig1gauss" << syst3 *100;
 

}

//fit_syst_error ends

void fit_syst_error_bin(TString fname, double bin_min, double bin_max){
  //prints the yield's systematic uncertainty per bin
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,fname);

  //  RooRealVar* Bpt  = ws->var("Bpt");
  RooDataSet* data = (RooDataSet*) ws->data("data");
  RooDataSet* data_bin;
  //select data subset corresponding to pT bin
  data_bin = (RooDataSet*) data->reduce(Form("Bpt>%lf",bin_min));
  data_bin = (RooDataSet*) data_bin->reduce(Form("Bpt<%lf",bin_max));
  //ws->import(*data_bin, Rename(Form("data_bin_%g_%g",bin_min,bin_max)));   

  build_pdf(*ws, "nominal");
  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres_nom = model->fitTo(*data_bin,Save());

  build_pdf(*ws,"bkg_poly");
  model = ws->pdf("model");
  RooFitResult* fitres_bgpol = model->fitTo(*data_bin,Save());

  build_pdf(*ws,"bkg_range");
  model = ws->pdf("model");
  RooFitResult* fitres_bgrange = model->fitTo(*data_bin,Range("peakright"),Save());

  build_pdf(*ws,"sig1gauss");
  model = ws->pdf("model");
  RooFitResult* fitres_sig1 = model->fitTo(*data_bin,Save());
    
  RooRealVar* n_nom = (RooRealVar*) fitres_nom  ->floatParsFinal().find("n_signal");
  RooRealVar* n_bgp = (RooRealVar*) fitres_bgpol->floatParsFinal().find("n_signal");
  RooRealVar* n_bra = (RooRealVar*) fitres_bgrange->floatParsFinal().find("n_signal");
  RooRealVar* n_sig1 = (RooRealVar*) fitres_sig1->floatParsFinal().find("n_signal");

  double n0  = n_nom->getVal();
  double n1  = n_bgp->getVal();
  double n2  = n_bra->getVal();
  double n3 = n_sig1->getVal();

  double syst1 = (n1-n0)/n0;
  double syst2 = (n2-n0)/n0;
  double syst3 = (n3-n0)/n0;
  
  cout << "bin_min"<<bin_min<< "_" << " bin_max" << bin_max << endl;
  cout << "syst bg pdf:" << std::setw(5) << syst1 * 100 << "\%\trange" << syst2 * 100 << "\%\tsig1gauss" << syst3 *100 << "\%" << endl;
 
}

//fit_syst_error_bin ends

double get_yield_syst(RooDataSet* data_bin, TString syst_src) {
  //returns the yield's value per bin
  
  //cout << "aaa 0\n";
  //data_bin->Print();
  //cout << "aaa 1\n";

  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  ws->import(*data_bin, Rename("data"));
  //  ws->Import(data_bin);

  //  TString mdl = (syst_src.Contains("range")) ? "nominal" : syst_src;
  TString rng = (syst_src.Contains("range")) ? "peakright" : "all";
  //cout << "bla 0\end ";
  build_pdf(*ws,syst_src.Data());

  //cout << "bla 1 \end ";
  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres = model->fitTo(*data_bin,Range(rng),Save());

  //build_pdf(*ws,"bkg_range");
  //model = ws->pdf("model");
  //RooFitResult* fitres_bgrange = model->fitTo(*data_bin,Range("peakright"),Save());

  RooRealVar* n1_var = (RooRealVar*) fitres ->floatParsFinal().find("n_signal");

  double n1  = n1_var->getVal();

  return n1; 
}


//get_yield_syst ends

void plot_complete_fit(RooWorkspace& w){

  RooAbsPdf*  model = w.pdf("model");
  RooDataSet* data = (RooDataSet*) w.data("data");
  data->Print();

  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* lambda   = w.var("lambda");

  model->fitTo(*data,Range("all"));

  RooPlot* massframe = Bmass.frame();

  if(particle == 0){
    data->plotOn(massframe, RooFit::Name("Data"));
    model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),Range("all"),LineColor(kOrange),LineStyle(kDashed)); 
    model->plotOn(massframe, RooFit::Name("B->J/psi X"),Components("erf"),Range("all"),LineColor(kGreen+3),LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("B->J/psi pi"),Components("jpsipi"),Range("all"),LineColor(kPink+10),LineStyle(kDashed));
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("Bmass (GeV)");
  }else if(particle == 1){
    data->plotOn(massframe, RooFit::Name("Data"));
    model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),Range("all"),LineColor(kOrange),LineStyle(kDashed));
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("Bmass (GeV)");
  }else if(particle == 2){
    data->plotOn(massframe, RooFit::Name("Data"));
    model->plotOn(massframe, RooFit::Name("Fit"), RooFit::Range("all"), RooFit::LineColor(kRed), RooFit::LineStyle(1), RooFit::LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Total Signal"), RooFit::Components("B0_signal"), RooFit::Range("all"), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("Corr Tag"), RooFit::Components("signal"), RooFit::Range("all"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Mis Tag"), RooFit::Components("k_pi_swap"), RooFit::Range("all"), RooFit::LineColor(kGray), RooFit::LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Combinatorial"), RooFit::Components("fit_side"), RooFit::Range("all"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
  }
  //model->paramOn(massframe,Layout(0.55,0.95,0.90));

  TCanvas d;
  d.SetTitle("");

  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
  p1->SetTitle("");
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw(); 
     
  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
  p2->SetTitle("");
  p2->SetTopMargin(0.); 
  p2->SetBottomMargin(0.2);   
  p2->SetBorderMode(1); 
  p2->Draw();

  p1->cd();
  massframe->Draw();

  TLatex* tex11 = new TLatex(0.6,0.8,"302.3 pb^{-1} (pp) 5.02 TeV");
  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  //tex11->Draw();
  tex11 = new TLatex(0.6,0.85,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  //tex11->Draw();
  
  double lambda_str = lambda->getVal();
  double lambda_err = lambda->getError();
  double chis = massframe->chiSquare();
 
  TLatex* tex12 = new TLatex(0.15, 0.85, Form("#lambda_{exp} = %.3lf #pm %.3lf",lambda_str,lambda_err));
  tex12->SetNDC(kTRUE);
  tex12->SetTextFont(42);
  tex12->SetTextSize(0.04);
 
  TLatex* tex13 = new TLatex(0.15, 0.8, Form("#chi/DOF = %.3lf",chis));
  tex13->SetNDC(kTRUE);
  tex13->SetTextFont(42);
  tex13->SetTextSize(0.04);
  
  TLegend *leg = new TLegend (0.7, 0.7, 0.9, 0.9);

  if(particle == 0){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/psi X", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/psi pi", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  }else if(particle == 1){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  }else if(particle == 2){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
    leg->AddEntry(massframe->findObject("Total Signal"), "Total sig", "l");
    leg->AddEntry(massframe->findObject("Corr Tag"), "Corr. sig", "l");
    leg->AddEntry(massframe->findObject("Mis Tag"), "Mis-tag sig", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Comb. bkg", "l");
  }
  leg->Draw("same");

  //pull dists

  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  RooPlot *pull_plot = Bmass.frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");

  pull_plot->GetXaxis()->SetTitle("");
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.15);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.13); 
 
  pull_plot->GetYaxis()->SetTitle("Pull hist");
  pull_plot->GetYaxis()->SetTitleFont(42);  
  pull_plot->GetYaxis()->SetTitleSize(0.10);
  pull_plot->GetYaxis()->SetTitleOffset(1.09);
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.13);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);
  pull_plot->GetYaxis()->SetNdivisions(305);

  p2->cd();
  pull_plot->Draw();
  
  if(particle == 0){
    //d.SaveAs("./results/Bu/complete_fit_Bu.pdf");
    d.SaveAs("./results/Bu/complete_fit_Bu.gif");
  }
  else if(particle == 1){
    //d.SaveAs("./results/Bs/complete_fit_Bs.pdf");
    d.SaveAs("./results/Bs/complete_fit_Bs.gif");
  }
  else if(particle == 2){
    //d.SaveAs("./results/B0/complete_fit_B0.pdf");
    d.SaveAs("./results/B0/complete_fit_B0.gif");
  }
}
//plot_complete_fit ends

//SIDEBAND SUBTRACTION//
std::vector<TH1D*> sideband_subtraction(RooWorkspace w, int* n, int n_var){
  
  RooDataSet* data = (RooDataSet*) w.data("data");

  RooAbsPdf* model = w.pdf("model");
  RooAbsPdf* BpModel = w.pdf("signal");
  RooAbsPdf* BgModel = w.pdf("fit_side");

  RooRealVar Bmass = *(w.var("Bmass"));

  vector<RooRealVar> variables;

  variables.push_back(*(w.var("Bmass"))); 
  variables.push_back(*(w.var("PVx")));
  variables.push_back(*(w.var("PVy")));
  variables.push_back(*(w.var("PVz")));
  variables.push_back(*(w.var("PVxE")));
  variables.push_back(*(w.var("PVyE")));
  variables.push_back(*(w.var("PVzE")));
  variables.push_back(*(w.var("BvtxX")));
  variables.push_back(*(w.var("BvtxY")));
  variables.push_back(*(w.var("BvtxZtoPVZ"))); 
  variables.push_back(*(w.var("BSx")));
  variables.push_back(*(w.var("BSy")));
  variables.push_back(*(w.var("BSz")));
  variables.push_back(*(w.var("BSxErr")));
  variables.push_back(*(w.var("BSyErr")));
  variables.push_back(*(w.var("BSzErr")));
  variables.push_back(*(w.var("BSdxdz")));
  variables.push_back(*(w.var("BSdydz")));
  variables.push_back(*(w.var("BSdxdzErr")));
  variables.push_back(*(w.var("BSdydzErr")));
  variables.push_back(*(w.var("BSWidthX")));
  variables.push_back(*(w.var("BSWidthXErr")));
  variables.push_back(*(w.var("BSWidthY")));
  variables.push_back(*(w.var("BSWidthYErr")));
  variables.push_back(*(w.var("By")));
  variables.push_back(*(w.var("Bpt")));
  variables.push_back(*(w.var("Btrk1Pt")));
  variables.push_back(*(w.var("Btrk1Eta")));
  variables.push_back(*(w.var("Btrk1PtErr")));
  variables.push_back(*(w.var("Bchi2cl")));
  variables.push_back(*(w.var("BsvpvDistance")));
  variables.push_back(*(w.var("BsvpvDisErr")));
  variables.push_back(*(w.var("BsvpvDistance_2D")));
  variables.push_back(*(w.var("BsvpvDisErr_2D")));
  variables.push_back(*(w.var("Bmumumass")));
  variables.push_back(*(w.var("Bmu1eta")));
  variables.push_back(*(w.var("Bmu2eta")));
  variables.push_back(*(w.var("Bmu1pt")));
  variables.push_back(*(w.var("Bmu2pt")));
  variables.push_back(*(w.var("Bmu1dxyPV")));
  variables.push_back(*(w.var("Bmu2dxyPV")));
  variables.push_back(*(w.var("Bmu1dzPV")));
  variables.push_back(*(w.var("Bmu2dzPV")));
  variables.push_back(*(w.var("Bd0")));
  variables.push_back(*(w.var("Bd0Err")));
  variables.push_back(*(w.var("Bdtheta")));
  variables.push_back(*(w.var("Balpha")));
  variables.push_back(*(w.var("Btrk1Dz1")));
  variables.push_back(*(w.var("Btrk1DzError1")));
  variables.push_back(*(w.var("Btrk1Dxy1")));
  variables.push_back(*(w.var("Btrk1DxyError1")));
  variables.push_back(*(w.var("Bmumueta")));
  variables.push_back(*(w.var("Bmumuphi")));
  variables.push_back(*(w.var("Bmumupt")));
  
  
  if(particle == 0){
    variables.push_back(*(w.var("BDT_pt_3_5")));
    variables.push_back(*(w.var("BDT_pt_5_7")));
    variables.push_back(*(w.var("BDT_pt_7_10")));
    variables.push_back(*(w.var("BDT_pt_10_15")));
    variables.push_back(*(w.var("BDT_pt_15_20")));
    variables.push_back(*(w.var("BDT_pt_20_50")));
    variables.push_back(*(w.var("BDT_pt_50_100")));
   }
  
  
  if(particle == 1){
    variables.push_back(*(w.var("Btrk2Pt")));
    variables.push_back(*(w.var("Btrk2Eta")));
    variables.push_back(*(w.var("Btrk2PtErr")));
    variables.push_back(*(w.var("Btrk2Dz1")));
    variables.push_back(*(w.var("Btrk2DzError1")));
    variables.push_back(*(w.var("Btrk2Dxy1")));
    variables.push_back(*(w.var("Btrk2DxyError1")));
    variables.push_back(*(w.var("BDT_pt_3_5")));
    variables.push_back(*(w.var("BDT_pt_5_7")));
    variables.push_back(*(w.var("BDT_pt_7_10")));
    variables.push_back(*(w.var("BDT_pt_10_15")));
    variables.push_back(*(w.var("BDT_pt_15_20")));
    variables.push_back(*(w.var("BDT_pt_20_50")));
    variables.push_back(*(w.var("BDT_pt_50_100")));
  } 
  
  
  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central; 

  double left = particle ? 5.3  : 5.2;
  double right = particle ? 5.4 : 5.4;

  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);
  Bmass.setRange("peakright",left,Bmass.getMax());
  Bmass.setRange("total", Bmass.getMin(), Bmass.getMax());
  
  reduceddata_side = particle ? (RooDataSet*)data->reduce(Form("Bmass>%lf || Bmass<%lf", right, left)) : (RooDataSet*)data->reduce(Form("Bmass>%lf", right));
  reduceddata_central = (RooDataSet*)data->reduce(Form("Bmass>%lf",left));
  reduceddata_central = (RooDataSet*)reduceddata_central->reduce(Form("Bmass<%lf",right));

  //Integrating the background distribution
   
  RooAbsReal* int_fit_side_right = BgModel->createIntegral(Bmass, Bmass, "right");
  RooAbsReal* int_fit_side_left = BgModel->createIntegral(Bmass, Bmass, "left");
  RooAbsReal* int_fit_peak = BgModel->createIntegral(Bmass, Bmass, "peak");
 
  cout << "Integral left band = " << int_fit_side_left->getVal() << endl; 
  cout << "Integral right band = " << int_fit_side_right->getVal() << endl;
  cout << "Integral peak = " << int_fit_peak->getVal() << endl; 
  cout << "normalisation = " << (BgModel->createIntegral(Bmass, Bmass, "total"))->getVal() << endl;

  /* 
  TCanvas *can = new TCanvas();
  RooPlot *plot = Bmass.frame();
  BgModel->plotOn(plot);
  plot->Draw();
  can->Draw();
  can->SaveAs("~/public/BinQGP/results/fit_side.gif"); 
  */

  double factor = particle ? (int_fit_peak->getVal())/(int_fit_side_right->getVal() + int_fit_side_left->getVal()) : (int_fit_peak->getVal())/(int_fit_side_right->getVal()); 

  std::cout << std::endl << "Factor: " << factor << std::endl;

  for(int i=0; i<n_var; i++){
    std::cout << "bins: " << n[i] << std::endl;
  } 
  std::vector<TH1D*> histos;

  if(particle == 0){
    histos.push_back(create_histogram(variables[1], "PVx",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "PVy",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "PVz",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "PVxE",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "PVyE",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "PVzE",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BvtxX",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BvtxY",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "BvtxZtoPVZ",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "BSx",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "BSy",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "BSz",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "BSxErr",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "BSyErr",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "BSzErr",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "BSdxdz",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "BSdydz",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "BSdxdzErr",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "BSdydzErr",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "BSWidthX",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "BSWidthXErr",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "BSWidthY",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "BSWidthYErr",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "By",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[29]));
    histos.push_back(create_histogram(variables[31], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[30]));
    histos.push_back(create_histogram(variables[32], "BsvpvDistance_2D",factor, reduceddata_side, reduceddata_central, data, n[31]));
    histos.push_back(create_histogram(variables[33], "BsvpvDisErr_2D",factor, reduceddata_side, reduceddata_central, data, n[32]));
    histos.push_back(create_histogram(variables[34], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[33]));
    histos.push_back(create_histogram(variables[35], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[34]));
    histos.push_back(create_histogram(variables[36], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[35]));
    histos.push_back(create_histogram(variables[37], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[36]));
    histos.push_back(create_histogram(variables[38], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[37]));
    histos.push_back(create_histogram(variables[39], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[38]));
    histos.push_back(create_histogram(variables[40], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[39]));
    histos.push_back(create_histogram(variables[41], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[40]));
    histos.push_back(create_histogram(variables[42], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[41]));
    histos.push_back(create_histogram(variables[43], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[42]));
    histos.push_back(create_histogram(variables[44], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[43]));
    histos.push_back(create_histogram(variables[45], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[44]));
    histos.push_back(create_histogram(variables[46], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[45]));
    histos.push_back(create_histogram(variables[47], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[46]));
    histos.push_back(create_histogram(variables[48], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[47]));
    histos.push_back(create_histogram(variables[49], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[48]));
    histos.push_back(create_histogram(variables[50], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[49]));
    histos.push_back(create_histogram(variables[51], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[50]));
    histos.push_back(create_histogram(variables[52], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[51]));
    histos.push_back(create_histogram(variables[53], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[52]));
    histos.push_back(create_histogram(variables[54], "BDT_pt_3_5",factor, reduceddata_side, reduceddata_central, data, n[53]));
    histos.push_back(create_histogram(variables[55], "BDT_pt_5_7",factor, reduceddata_side, reduceddata_central, data, n[54]));
    histos.push_back(create_histogram(variables[56], "BDT_pt_7_10",factor, reduceddata_side, reduceddata_central, data, n[55]));
    histos.push_back(create_histogram(variables[57], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[56]));
    histos.push_back(create_histogram(variables[58], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[57]));
    histos.push_back(create_histogram(variables[59], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[58]));
    histos.push_back(create_histogram(variables[60], "BDT_pt_50_100",factor, reduceddata_side, reduceddata_central, data, n[59]));
  }else if(particle == 1){
    histos.push_back(create_histogram(variables[1], "PVx",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "PVy",factor, reduceddata_side, reduceddata_central, data, n[1])); 
    histos.push_back(create_histogram(variables[3], "PVz",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "PVxE",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "PVyE",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "PVzE",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BvtxX",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BvtxY",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "BvtxZtoPVZ",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "BSx",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "BSy",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "BSz",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "BSxErr",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "BSyErr",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "BSzErr",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "BSdxdz",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "BSdydz",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "BSdxdzErr",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "BSdydzErr",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "BSWidthX",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "BSWidthXErr",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "BSWidthY",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "BSWidthYErr",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "By",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[29]));
    histos.push_back(create_histogram(variables[31], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[30]));
    histos.push_back(create_histogram(variables[32], "BsvpvDistance_2D",factor, reduceddata_side, reduceddata_central, data, n[31]));
    histos.push_back(create_histogram(variables[33], "BsvpvDisErr_2D",factor, reduceddata_side, reduceddata_central, data, n[32]));
    histos.push_back(create_histogram(variables[34], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[33]));
    histos.push_back(create_histogram(variables[35], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[34]));
    histos.push_back(create_histogram(variables[36], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[35]));
    histos.push_back(create_histogram(variables[37], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[36]));
    histos.push_back(create_histogram(variables[38], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[37]));
    histos.push_back(create_histogram(variables[39], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[38]));
    histos.push_back(create_histogram(variables[40], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[39]));
    histos.push_back(create_histogram(variables[41], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[40]));
    histos.push_back(create_histogram(variables[42], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[41]));
    histos.push_back(create_histogram(variables[43], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[42]));
    histos.push_back(create_histogram(variables[44], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[43]));
    histos.push_back(create_histogram(variables[45], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[44]));
    histos.push_back(create_histogram(variables[46], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[45]));
    histos.push_back(create_histogram(variables[47], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[46]));
    histos.push_back(create_histogram(variables[48], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[47]));
    histos.push_back(create_histogram(variables[49], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[48]));
    histos.push_back(create_histogram(variables[50], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[49]));
    histos.push_back(create_histogram(variables[51], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[50]));
    histos.push_back(create_histogram(variables[52], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[51]));
    histos.push_back(create_histogram(variables[53], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[52]));
    histos.push_back(create_histogram(variables[54], "Btrk2Pt",factor, reduceddata_side, reduceddata_central, data, n[53]));
    histos.push_back(create_histogram(variables[55], "Btrk2Eta",factor, reduceddata_side, reduceddata_central, data, n[54]));
    histos.push_back(create_histogram(variables[56], "Btrk2PtErr",factor, reduceddata_side, reduceddata_central, data, n[55]));
    histos.push_back(create_histogram(variables[57], "Btrk2Dz1",factor, reduceddata_side, reduceddata_central, data, n[56]));
    histos.push_back(create_histogram(variables[58], "Btrk2DzError1",factor, reduceddata_side, reduceddata_central, data, n[57]));
    histos.push_back(create_histogram(variables[59], "Btrk2Dxy1",factor, reduceddata_side, reduceddata_central, data, n[58]));
    histos.push_back(create_histogram(variables[60], "Btrk2DxyError1",factor, reduceddata_side, reduceddata_central, data, n[59]));
    histos.push_back(create_histogram(variables[61], "BDT_pt_3_5",factor, reduceddata_side, reduceddata_central, data, n[60]));
    histos.push_back(create_histogram(variables[62], "BDT_pt_5_7",factor, reduceddata_side, reduceddata_central, data, n[61]));
    histos.push_back(create_histogram(variables[63], "BDT_pt_7_10",factor, reduceddata_side, reduceddata_central, data, n[62]));
    histos.push_back(create_histogram(variables[64], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[63]));
    histos.push_back(create_histogram(variables[65], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[64]));
    histos.push_back(create_histogram(variables[66], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[65]));
    histos.push_back(create_histogram(variables[67], "BDT_pt_50_100",factor, reduceddata_side, reduceddata_central, data, n[66])); 
  }

  return histos;
  
}

//sideband_subtraction ends


TH1D* create_histogram_mc(RooRealVar var, TTree* t, int n, TString weight){

  TH1D* h = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());
  TH1D* wei = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());

  TString name_string = TString(var.GetName()) + ">>htemp(" + Form("%d",n) +"," + Form("%lf", var.getMin()) + "," + Form("%lf", var.getMax()) + ")";

  t->Draw(name_string, weight);

  h = (TH1D*)gDirectory->Get("htemp")->Clone();
  h->SetTitle("");
  h->SetMarkerStyle(29);
  h->SetMarkerColor(kGreen);
  h->SetMarkerSize(1);
  h->SetLineColor(kGreen);
  return h;
}

//create_histogram_mc ends

TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n){

  std::cout<< "n in create_histogram = "<< n <<std::endl;
  TH1D* dist_side = (TH1D*)reduced->createHistogram("dist_side",var, Binning(n, var.getMin(), var.getMax()));
  dist_side->SetMarkerColor(kRed);
  dist_side->SetLineColor(kRed);
  dist_side->SetNameTitle("dist_side", "");

  TH1D* hist_dist_peak = (TH1D*)central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
  TH1D* dist_peak = new TH1D(*hist_dist_peak);
  dist_peak->SetMarkerColor(kBlue);
  dist_peak->SetLineColor(kBlue);
  dist_peak->SetNameTitle("dist_peak", "");

  hist_dist_peak->SetMarkerColor(kBlack);
  hist_dist_peak->SetLineColor(kBlack);
  hist_dist_peak->SetNameTitle("dist_total", "");

  dist_peak->Add(dist_side, -factor);
  dist_side->Scale(factor);

  dist_peak->SetStats(0);
  dist_side->SetStats(0);
  hist_dist_peak->SetStats(0);
  TCanvas c;

  hist_dist_peak->Draw();
  dist_side->Draw("same");
  dist_peak->Draw("same");
  
  dist_peak->SetXTitle(var.GetName());
  dist_side->SetXTitle(var.GetName());
  hist_dist_peak->SetXTitle(var.GetName());

  hist_dist_peak->GetYaxis()->SetRangeUser(0, 1.3*hist_dist_peak->GetMaximum());
  TLatex* tex = new TLatex(0.6,0.8,"302.3 pb^{-1} (pp) 5.02 TeV");
  /*tex->SetNDC(kTRUE);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.04);
  tex->Draw();
  tex = new TLatex(0.68,0.85,"CMS Preliminary");
  tex->SetNDC(kTRUE);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();*/

  TLegend *leg = new TLegend (0.7, 0.9, 0.9, 1.0);
  leg->AddEntry("dist_peak", "Signal", "l");
  leg->AddEntry("dist_side", "Background", "l");
  leg->AddEntry("hist_dist_peak", "Total", "l");
  leg->Draw("same");

  std::cout<<"name: "<<var.GetName()<<std::endl;
  std::cout<<"histo name: "<<dist_peak->GetName()<<std::endl;

  if(particle == 0){
    c.SaveAs("./results/Bu/sideband_sub/"+name + "sideband_sub_Bu.pdf");
    c.SaveAs("./results/Bu/sideband_sub/"+name + "sideband_sub_Bu.gif");
  }else if(particle == 1){
    c.SaveAs("./results/Bs/sideband_sub/"+name + "sideband_sub_Bs.pdf");
    c.SaveAs("./results/Bs/sideband_sub/"+name + "sideband_sub_Bs.gif");
  }

  if(background == 0){return dist_peak;}
  else if(background == 1){return dist_side;}

}
//create_histogram ends


void do_splot(RooWorkspace& w){

  RooDataSet* data = (RooDataSet*) w.data("data");   
  RooAbsPdf* model = w.pdf("model");
  //we need the fit and the dataset previously saved in the woorkspace

  RooRealVar* BpYield = w.var("n_signal");
  RooRealVar* BgYield = w.var("n_combinatorial");
  //we need the n values previously saved in the woorkspace

  //fit the model to the data
  model->fitTo(*data,Extended());

  //sPlot technique requires model parameters (other than the yields) to be fixed

  RooRealVar* mean  = w.var("mean");
  RooRealVar* sigma1 = w.var("sigma1");
  RooRealVar* sigma2 = w.var("sigma2");
  RooRealVar* cofs = w.var("cofs");
  RooRealVar* lambda = w.var("lambda");

  mean->setConstant();
  sigma1->setConstant();
  sigma2->setConstant();
  cofs->setConstant();
  lambda->setConstant();

  RooMsgService::instance().setSilentMode(true);

  //add sWeights to dataset based on model and yield variables
  //sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
  SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*BpYield,*BgYield));


  cout << endl <<  "Yield of B+ is "
       << BpYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_signal") << endl;

  cout << "Yield of background is "
       << BgYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;

  //for(Int_t i=0; i < 10; i++) {
  //if(0)
  // cout << "y Weight   "     << sData->GetSWeight(i,"BpYield")
  // << "\tb Weight   "     << sData->GetSWeight(i,"BgYield")
  // << "\ttotal Weight   " << sData->GetSumOfEventSWeight(i)
  // << endl;

  //}

  //cout << endl

  w.import(*data, Rename("dataWithSWeights"));
  //the reweighted data is saved in the woorkspace 
}

//do_splot ends

TH1D* make_splot(RooWorkspace& w, int n, TString label){

  //saves the plots of signal distributions, background distributions and signal+background distributions
  //in the end returns the histogram of signal

  TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
  cdata->Divide(2,2);

  RooAbsPdf* model  = w.pdf("model");
  RooAbsPdf* BpModel = w.pdf("signal");
  RooAbsPdf* BgModel = w.pdf("fit_side");

  RooRealVar* Bmass  = w.var("Bmass");
  RooRealVar* variable = w.var(label);

  RooRealVar* BpYield = w.var("n_signal");
  RooRealVar* BgYield = w.var("n_combinatorial");

  double sigYield = BpYield->getVal();
  double bkgYield = BgYield->getVal();
  
  RooDataSet* data = (RooDataSet*) w.data("data");

  cdata->cd(1);
  RooPlot* mframe = Bmass->frame();
  if(particle == 0){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of B+ [GeV]"));
  }else if(particle == 1){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of Bs [GeV]"));
  }
  data->plotOn(mframe);
  model->plotOn(mframe,LineColor(kRed));
  model->plotOn(mframe,Components(*BpModel),LineStyle(kDashed),LineColor(kOrange));
  model->plotOn(mframe,Components(*BgModel),LineStyle(kDashed),LineColor(kBlue));
  mframe->SetTitle("Bmass");
  mframe->Draw();

  cdata->cd(2);
  RooPlot* ptframe = variable->frame();
  data->plotOn(ptframe);
  if(particle == 0){
    ptframe->SetTitle(label + " of B+: total sample");
  }else if(particle == 1){
    ptframe->SetTitle(label + " of Bs: total sample");
  }
  ptframe->Draw();

  //get the dataset with sWeights
  RooDataSet* dataW = (RooDataSet*) w.data("dataWithSWeights");
  RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
  w.import(*dataWBp,Rename("dataWBp"));
  RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

  RooPlot* ptframe2Bp = variable->frame();
  RooPlot* ptframe2Bg = variable->frame();

  ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/n));
  ptframe2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/n));

  if(particle == 0){
    ptframe2Bp->GetXaxis()->SetTitle(label + " of B+");
    ptframe2Bg->GetXaxis()->SetTitle(label + " of B+");
  }else if(particle == 1){
    ptframe2Bp->GetXaxis()->SetTitle(label + " of Bs");
    ptframe2Bg->GetXaxis()->SetTitle(label + " of Bs");
  }

  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));
  dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(n));

  if(particle == 0){
    ptframe2Bp->SetTitle(label+" distribution of B+ for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of B+ for background (splot)");
  }else if(particle == 1){
    ptframe2Bp->SetTitle(label+" distribution of Bs for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of Bs for background (splot)");
  }

  cdata->cd(3);  ptframe2Bp->Draw();
  cdata->cd(4);  ptframe2Bg->Draw();

  if(particle == 0){
    cdata->SaveAs("./results/Bu/splot/Bmass/"+label+"sPlot_Bu.gif");
    cdata->SaveAs("./results/Bu/splot/Bmass/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    cdata->SaveAs("./results/Bs/splot/Bmass/"+label+"sPlot_Bs.gif");
    cdata->SaveAs("./results/Bs/splot/Bmass/"+label+"sPlot_Bs.pdf");
  }

  TH1D* histo_Bp_sig = (TH1D*)dataWBp->createHistogram(label,n,0,0);
  TH1D* histo_Bp_bkg = (TH1D*)dataWBg->createHistogram(label,n,0,0);

    for (int i=1; i<=n; i++) {

      if (histo_Bp_sig->GetBinContent(i)==0) histo_Bp_sig->SetBinError(i,0.);
      if (histo_Bp_bkg->GetBinContent(i)==0) histo_Bp_bkg->SetBinError(i,0.);

       histo_Bp_sig->SetBinContent(i,histo_Bp_sig->GetBinContent(i)/sigYield);
       histo_Bp_sig->SetBinError(i,histo_Bp_sig->GetBinError(i)/sigYield);

       histo_Bp_bkg->SetBinContent(i,histo_Bp_bkg->GetBinContent(i)/bkgYield);
       histo_Bp_bkg->SetBinError(i,histo_Bp_bkg->GetBinError(i)/bkgYield);

  }

  TCanvas* prov = new TCanvas ("prov","c1",200,10,700,500);
  prov->cd();
  //histo_Bp_sig->SetMarkerStyle(20);
  histo_Bp_sig->SetMarkerSize(1);
  histo_Bp_sig->SetMarkerColor(kRed);
  histo_Bp_sig->SetLineColor(kRed);
  histo_Bp_sig->SetTitle("");
  histo_Bp_sig->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/n));
  histo_Bp_sig->GetXaxis()->SetTitle(label);

  histo_Bp_sig->SetStats(0);
  histo_Bp_sig->Draw("E");

  if(particle == 0){
    prov->SaveAs("./results/Bu/splot/sig/"+label+"sPlot_Bu.gif");
    prov->SaveAs("./results/Bu/splot/sig/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    prov->SaveAs("./results/Bs/splot/sig/"+label+"sPlot_Bs.gif");
    prov->SaveAs("./results/Bs/splot/sig/"+label+"sPlot_Bs.pdf");
  }

  TCanvas* prov_bkg = new TCanvas ("prov_bkg","c2",200,10,700,500);
  prov_bkg->cd();
  histo_Bp_bkg->SetMarkerStyle(20);
  histo_Bp_bkg->SetMarkerSize(0.);
  histo_Bp_bkg->SetMarkerColor(kBlue);
  histo_Bp_bkg->SetLineColor(kBlue);
  histo_Bp_bkg->SetTitle("");
  histo_Bp_bkg->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/n));
  histo_Bp_bkg->GetXaxis()->SetTitle(label);

  histo_Bp_bkg->SetStats(0);
  histo_Bp_bkg->Draw("E");

  if(particle == 0){
    prov_bkg->SaveAs("./results/Bu/splot/bkg/"+label+"sPlot_Bu.gif");
    prov_bkg->SaveAs("./results/Bu/splot/bkg/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    prov_bkg->SaveAs("./results/Bs/splot/bkg/"+label+"sPlot_Bs.gif");
    prov_bkg->SaveAs("./results/Bs/splot/bkg/"+label+"sPlot_Bs.pdf");
  }

  TCanvas* sig_bkg = new TCanvas ("sig_bkg","c3",200,10,700,500); 
  sig_bkg->cd();


  histo_Bp_sig->Draw();
  histo_Bp_bkg->Draw("same");

  //y axis: maximum and minimum 
  if (histo_Bp_bkg->GetMaximum() > histo_Bp_sig->GetMaximum()){
    histo_Bp_sig->GetYaxis()->SetRangeUser(0.1*histo_Bp_sig->GetMinimum(), 1.1*histo_Bp_bkg->GetMaximum());
}
  else {
    histo_Bp_sig->GetYaxis()->SetRangeUser(0.1*histo_Bp_sig->GetMinimum(), 1.1*histo_Bp_sig->GetMaximum());
  }

  TLegend* legend = new TLegend(0.7,0.9,0.9,1.0);
  legend->AddEntry(histo_Bp_sig,"Signal","lep");
  legend->AddEntry(histo_Bp_bkg,"Background","lep");
  legend->Draw();

  if(particle == 0){
    sig_bkg->SaveAs("./results/Bu/splot/sig_bkg/"+label+"sPlot_Bu.gif");
    sig_bkg->SaveAs("./results/Bu/splot/sig_bkg/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    sig_bkg->SaveAs("./results/Bs/splot/sig_bkg/"+label+"sPlot_Bs.gif");
    sig_bkg->SaveAs("./results/Bs/splot/sig_bkg/"+label+"sPlot_Bs.pdf");
  }

  //cleanup
  delete cdata;
  delete prov;
  delete prov_bkg;
  delete sig_bkg;

  if(background == 0){return histo_Bp_sig;}
  else if(background == 1){return histo_Bp_bkg;}
} 

//make_splot ends

//SPLOT_METHOD//

std::vector<TH1D*> splot_method(RooWorkspace& w, int* n, TString* label, int n_var){

  std::vector<TH1D*> histos;

  for(int i = 0;i<n_var;i++){
    histos.push_back(make_splot(w,n[i],label[i]));
  }

  return histos;
}

void validate_fit(RooWorkspace* w)
{

  RooRealVar Bmass = *(w->var("Bmass"));
  RooAbsPdf* model  = w->pdf("model");
  
  vector<RooRealVar> params;
  params.push_back(*(w->var("n_signal")));

  int params_size = params.size();  

  RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass, Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

  mcstudy->generateAndFit(5000);
  vector<RooPlot*> framesPull, framesParam;

  for(int i = 0; i < params_size; ++i)
    {
      framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(200)));//FrameRange(-5,5)
      framesPull[i]->SetTitle("");
      framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
      framesParam[i]->SetTitle("");
    }

  vector<TGraph*> h;

  for(int i = 0; i < params_size; ++i){
    h.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
  }

  gStyle->SetOptFit(0111);

  TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

  gPad->SetLeftMargin(0.15);

  for(int i = 0; i < params_size; ++i){
    c_pull->cd();
    h[i]->SetTitle("");
    h[i]->Draw();
    c_pull->Update();
    h[i]->Fit("gaus","","",-20,20);
    h[i]->GetFunction("gaus")->SetLineColor(4);
    h[i]->GetFunction("gaus")->SetLineWidth(5);
    h[i]->GetXaxis()->SetTitle("Pull");
    h[i]->GetYaxis()->SetTitle("Toy MCs");
    h[i]->GetXaxis()->SetRangeUser(-10,10);
    h[i]->Draw("same");
  }

  TCanvas* c_params = new TCanvas("params", "params", 900, 800);

  for(int i = 0; i < params_size; ++i){
    c_params->cd();
    framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
    framesParam.at(i)->Draw();
  }

  if(particle == 0){
    //c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.pdf");
    c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.gif");
    //c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.pdf");
    c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.gif");
  }else if(particle == 1){
    //c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.pdf");
    c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.gif");
    //c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.pdf");
    c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.gif");
  }else if(particle == 2){
    //c_pull->SaveAs("./results/B0/pulls/pulls_poisson_B0.pdf");
    c_pull->SaveAs("./results/B0/pulls/pulls_poisson_B0.gif");
    //c_params->SaveAs("./results/B0/pulls/pulls_params_poisson_B0.pdf");
    c_params->SaveAs("./results/B0/pulls/pulls_params_poisson_B0.gif");
  }  
}

void AddWeights(TTree* t){
  TString input_file;
  input_file = particle ? "~/public/BinQGP/results/Bs/mc_validation_plots/weights/weights.root": "~/public/BinQGP/results/Bu/mc_validation_plots/weights/weights.root";
  TFile* f_wei = new TFile(input_file, "read");
 
  TH1D* h_bdt_pt_3_5 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_3_5;1"));
  TH1D* h_bdt_pt_5_7 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_5_7;1"));
  TH1D* h_bdt_pt_7_10 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_7_10;1"));
  TH1D* h_bdt_pt_10_15 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_10_15;1"));
  TH1D* h_bdt_pt_15_20 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_15_20;1"));
  TH1D* h_bdt_pt_20_50 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_20_50;1"));
  TH1D* h_bdt_pt_50_100 = (TH1D*)f_wei->Get(TString("weights_BDT_pt_50_100;1"));

 
  double bdt_pt_3_5_min = h_bdt_pt_3_5->GetXaxis()->GetXmin();
  double bdt_pt_5_7_min = h_bdt_pt_5_7->GetXaxis()->GetXmin();
  double bdt_pt_7_10_min = h_bdt_pt_7_10->GetXaxis()->GetXmin();
  double bdt_pt_10_15_min = h_bdt_pt_10_15->GetXaxis()->GetXmin();
  double bdt_pt_15_20_min = h_bdt_pt_15_20->GetXaxis()->GetXmin();
  double bdt_pt_20_50_min = h_bdt_pt_20_50->GetXaxis()->GetXmin();
  double bdt_pt_50_100_min = h_bdt_pt_50_100->GetXaxis()->GetXmin();

  double bdt_pt_3_5_max = h_bdt_pt_3_5->GetXaxis()->GetXmax();
  double bdt_pt_5_7_max = h_bdt_pt_5_7->GetXaxis()->GetXmax();
  double bdt_pt_7_10_max = h_bdt_pt_7_10->GetXaxis()->GetXmax();
  double bdt_pt_10_15_max = h_bdt_pt_10_15->GetXaxis()->GetXmax();
  double bdt_pt_15_20_max = h_bdt_pt_15_20->GetXaxis()->GetXmax();
  double bdt_pt_20_50_max = h_bdt_pt_20_50->GetXaxis()->GetXmax();
  double bdt_pt_50_100_max = h_bdt_pt_50_100->GetXaxis()->GetXmax();


  TFile* f_tree = particle ? new TFile("./results/Bs/mc_validation_plots/weights/tree_with_weight.root", "recreate") : new TFile("./results/Bu/mc_validation_plots/weights/tree_with_weight.root", "recreate");

  f_tree->cd();
  TTree* tw = t->CloneTree();
  //TTree *tw = new TTree("tw", "tw");
  
  Float_t Bpt;
  Float_t By;
  Float_t Btrk1Pt;
  Float_t Btrk1Eta;
  Float_t Btrk1PtErr;
  Float_t Bchi2cl;
  Float_t BsvpvDistance;
  Float_t BsvpvDistance_2D;
  Float_t BsvpvDisErr_2D;
  Float_t BsvpvDisErr;
  Float_t Bmumumass;
  Float_t Bmu1eta;
  Float_t Bmu2eta;
  Float_t Bmu1pt;
  Float_t Bmu2pt;
  Float_t Bmu1dxyPV;
  Float_t Bmu2dxyPV;
  Float_t Bmu1dzPV;
  Float_t Bmu2dzPV;
  Float_t Bd0;
  Float_t Bd0Err;
  Float_t Bdtheta;
  Float_t Balpha;
  Float_t Btrk1Dz1;
  Float_t Btrk1DzError1;
  Float_t Btrk1Dxy1;
  Float_t Btrk1DxyError1;
  Float_t Bmumueta;
  Float_t Bmumuphi;
  Float_t Bmumupt;
  Float_t BDT_pt_3_5;
  Float_t BDT_pt_5_7;
  Float_t BDT_pt_7_10;
  Float_t BDT_pt_10_15;
  Float_t BDT_pt_15_20;
  Float_t BDT_pt_20_50;
  Float_t BDT_pt_50_100;

  Float_t weight_BDT_pt_3_5;
  Float_t weight_BDT_pt_5_7;
  Float_t weight_BDT_pt_7_10;
  Float_t weight_BDT_pt_10_15;
  Float_t weight_BDT_pt_15_20;
  Float_t weight_BDT_pt_20_50;
  Float_t weight_BDT_pt_50_100;

  t->SetBranchAddress("Bpt", &Bpt);
  t->SetBranchAddress("By", &By);
  t->SetBranchAddress("Btrk1Pt", &Btrk1Pt);
  t->SetBranchAddress("Btrk1Eta", &Btrk1Eta);
  t->SetBranchAddress("Btrk1PtErr", &Btrk1PtErr);
  t->SetBranchAddress("Bchi2cl", &Bchi2cl);
  t->SetBranchAddress("BsvpvDistance", &BsvpvDistance);
  t->SetBranchAddress("BsvpvDistance_2D", &BsvpvDistance_2D);
  t->SetBranchAddress("BsvpvDisErr", &BsvpvDisErr);
  t->SetBranchAddress("BsvpvDisErr_2D", &BsvpvDisErr_2D);
  t->SetBranchAddress("Bmumumass", &Bmumumass);
  t->SetBranchAddress("Bmu1eta", &Bmu1eta);
  t->SetBranchAddress("Bmu2eta", &Bmu2eta);
  t->SetBranchAddress("Bmu1pt", &Bmu1pt);
  t->SetBranchAddress("Bmu2pt", &Bmu2pt);
  t->SetBranchAddress("Bmu1dxyPV", &Bmu1dxyPV);
  t->SetBranchAddress("Bmu2dxyPV", &Bmu2dxyPV);
  t->SetBranchAddress("Bmu1dzPV", &Bmu1dzPV);
  t->SetBranchAddress("Bmu2dzPV", &Bmu2dzPV);
  t->SetBranchAddress("Bd0", &Bd0);
  t->SetBranchAddress("Bd0Err", &Bd0Err);
  t->SetBranchAddress("Bdtheta", &Bdtheta);
  t->SetBranchAddress("Balpha", &Balpha);
  t->SetBranchAddress("Btrk1Dz1", &Btrk1Dz1);
  t->SetBranchAddress("Btrk1DzError1", &Btrk1DzError1);
  t->SetBranchAddress("Btrk1Dxy1", &Btrk1Dxy1);
  t->SetBranchAddress("Btrk1DxyError1", &Btrk1DxyError1);
  t->SetBranchAddress("Bmumueta", &Bmumueta);
  t->SetBranchAddress("Bmumuphi", &Bmumuphi);
  t->SetBranchAddress("Bmumupt", &Bmumupt);
  t->SetBranchAddress("BDT_pt_3_5", &BDT_pt_3_5);
  t->SetBranchAddress("BDT_pt_5_7", &BDT_pt_5_7);
  t->SetBranchAddress("BDT_pt_7_10", &BDT_pt_7_10);
  t->SetBranchAddress("BDT_pt_10_15", &BDT_pt_10_15);
  t->SetBranchAddress("BDT_pt_15_20", &BDT_pt_15_20);
  t->SetBranchAddress("BDT_pt_20_50", &BDT_pt_20_50);
  t->SetBranchAddress("BDT_pt_50_100", &BDT_pt_50_100);

  tw->Branch("weight_BDT_pt_3_5", &weight_BDT_pt_3_5);
  tw->Branch("weight_BDT_pt_5_7", &weight_BDT_pt_5_7);
  tw->Branch("weight_BDT_pt_7_10", &weight_BDT_pt_7_10);
  tw->Branch("weight_BDT_pt_10_15", &weight_BDT_pt_10_15);
  tw->Branch("weight_BDT_pt_15_20", &weight_BDT_pt_15_20);
  tw->Branch("weight_BDT_pt_20_50", &weight_BDT_pt_20_50);
  tw->Branch("weight_BDT_pt_50_100", &weight_BDT_pt_50_100);
  
  cout << "Starting Cycle" << endl;

  for (int i = 0; i < (t->GetEntries()); i++){
    t->GetEntry(i);

    if(BDT_pt_3_5 >= bdt_pt_3_5_min && BDT_pt_3_5 <= bdt_pt_3_5_max){weight_BDT_pt_3_5 = h_bdt_pt_3_5->GetBinContent(h_bdt_pt_3_5->FindBin(BDT_pt_3_5));}
    else{weight_BDT_pt_3_5 = 1.0;}

    if(BDT_pt_5_7 >= bdt_pt_5_7_min && BDT_pt_5_7 <= bdt_pt_5_7_max){weight_BDT_pt_5_7 = h_bdt_pt_5_7->GetBinContent(h_bdt_pt_5_7->FindBin(BDT_pt_5_7));}
    else{weight_BDT_pt_5_7 = 1.0;}

    if(BDT_pt_7_10 >= bdt_pt_7_10_min && BDT_pt_7_10 <= bdt_pt_7_10_max){weight_BDT_pt_7_10 = h_bdt_pt_7_10->GetBinContent(h_bdt_pt_7_10->FindBin(BDT_pt_7_10));}
    else{weight_BDT_pt_7_10 = 1.0;}

    if(BDT_pt_10_15 >= bdt_pt_10_15_min && BDT_pt_10_15 <= bdt_pt_10_15_max){weight_BDT_pt_10_15 = h_bdt_pt_10_15->GetBinContent(h_bdt_pt_10_15->FindBin(BDT_pt_10_15));}
    else{weight_BDT_pt_10_15 = 1.0;}

    if(BDT_pt_15_20 >= bdt_pt_15_20_min && BDT_pt_15_20 <= bdt_pt_15_20_max){weight_BDT_pt_15_20 = h_bdt_pt_15_20->GetBinContent(h_bdt_pt_15_20->FindBin(BDT_pt_15_20));}
    else{weight_BDT_pt_15_20 = 1.0;}

    if(BDT_pt_20_50 >= bdt_pt_20_50_min && BDT_pt_20_50 <= bdt_pt_20_50_max){weight_BDT_pt_20_50 = h_bdt_pt_20_50->GetBinContent(h_bdt_pt_20_50->FindBin(BDT_pt_20_50));}
    else{weight_BDT_pt_20_50 = 1.0;}

    if(BDT_pt_50_100 >= bdt_pt_50_100_min && BDT_pt_50_100 <= bdt_pt_50_100_max){weight_BDT_pt_50_100 = h_bdt_pt_50_100->GetBinContent(h_bdt_pt_50_100->FindBin(BDT_pt_50_100));}
    else{weight_BDT_pt_50_100 = 1.0;}

    tw->Fill();
    tw->Show();
  }  

  cout << "Ending Cycle" << endl;

  f_tree->Write();
  //tw->Show();
  f_tree->Close();
  f_wei->Close();
}
//Weights were added to TTree

void set_up_workspace_variables(RooWorkspace& w)
{

  if(particle == 0){
    double PVx_min, PVx_max;
    double PVy_min, PVy_max;
    double PVz_min, PVz_max;
    double PVxE_min, PVxE_max;
    double PVyE_min, PVyE_max;
    double PVzE_min, PVzE_max;
    double BvtxX_min, BvtxX_max;
    double BvtxY_min, BvtxY_max;
    double BvtxZtoPVZ_min, BvtxZtoPVZ_max;
    double BSx_min, BSx_max;
    double BSy_min, BSy_max;
    double BSz_min, BSz_max;
    double BSxErr_min, BSxErr_max;
    double BSyErr_min, BSyErr_max;
    double BSzErr_min, BSzErr_max;
    double BSdxdz_min, BSdxdz_max;
    double BSdydz_min, BSdydz_max;
    double BSdxdzErr_min, BSdxdzErr_max;
    double BSdydzErr_min, BSdydzErr_max;
    double BSWidthX_min, BSWidthX_max;
    double BSWidthXErr_min, BSWidthXErr_max;
    double BSWidthY_min, BSWidthY_max;
    double BSWidthYErr_min, BSWidthYErr_max;
    double mass_min, mass_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double trk1pt_min, trk1pt_max;
    double trk1eta_min, trk1eta_max;
    double trk1pterr_min, trk1pterr_max;
    double chi2cl_min, chi2cl_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDisErr_min, svpvDisErr_max;
    double svpvDistance2D_min, svpvDistance2D_max;
    double svpvDisErr2D_min, svpvDisErr2D_max;
    double mumumass_min, mumumass_max; 
    double mu1eta_min, mu1eta_max;
    double mu2eta_min, mu2eta_max;
    double mu1pt_min, mu1pt_max;
    double mu2pt_min, mu2pt_max;
    double mu1dxyPV_min, mu1dxyPV_max;
    double mu2dxyPV_min, mu2dxyPV_max;
    double mu1dzPV_min, mu1dzPV_max;
    double mu2dzPV_min, mu2dzPV_max;
    double d0_min, d0_max;
    double d0err_min, d0err_max;
    double dtheta_min, dtheta_max;
    double alpha_min, alpha_max;
    double trk1Dz1_min, trk1Dz1_max;
    double trk1DzError1_min, trk1DzError1_max;
    double trk1Dxy1_min, trk1Dxy1_max;
    double trk1DxyError1_min, trk1DxyError1_max;
    double mumueta_min, mumueta_max;
    double mumuphi_min, mumuphi_max;
    double mumupt_min, mumupt_max;
    double BDT_3_5_min, BDT_3_5_max;
    double BDT_5_7_min, BDT_5_7_max;
    double BDT_7_10_min, BDT_7_10_max;
    double BDT_10_15_min, BDT_10_15_max;
    double BDT_15_20_min, BDT_15_20_max;
    double BDT_20_50_min, BDT_20_50_max;
    double BDT_50_100_min, BDT_50_100_max;

    PVx_min = -0.1;
    PVx_max = 0.1;

    PVy_min = -0.1;
    PVy_max = 0.1;

    PVz_min = 0.;
    PVz_max = 1.5;

    PVxE_min = -0.001;
    PVxE_max = 0.001;

    PVyE_min = -0.001;
    PVyE_max = 0.001;

    PVzE_min = 0.;
    PVzE_max = 0.3;

    BvtxX_min = -1.;
    BvtxX_max = 1.;

    BvtxY_min = -1.;
    BvtxY_max = 1.;

    BvtxZtoPVZ_min = 0.;
    BvtxZtoPVZ_max = 25.;

    BSx_min = -0.1;
    BSx_max = 0.1;

    BSy_min = -0.1;
    BSy_max = 0.1;

    BSz_min = 0.;
    BSz_max = 1.5;    

    BSxErr_min = -0.001;
    BSxErr_max = 0.001;
    
    BSyErr_min = -0.001;
    BSyErr_max = 0.001;
    
    BSzErr_min = 0.;
    BSzErr_max = 0.3;
    
    BSdxdz_min = -0.001;
    BSdxdz_max = 0.002;
    
    BSdydz_min = -0.004;
    BSdydz_max = 0.004;
    
    BSdxdzErr_min = -0.001;
    BSdxdzErr_max = 0.004;
    
    BSdydzErr_min = -0.001;
    BSdydzErr_max = 0.004;
    
    BSWidthX_min = -0.001;
    BSWidthX_max = 0.005;
    
    BSWidthXErr_min = -0.001;
    BSWidthXErr_max = 0.003;
    
    BSWidthY_min = -0.001;
    BSWidthY_max = 0.005;
    
    BSWidthYErr_min = -0.001;
    BSWidthYErr_max = 0.003;

    mass_min = 5.0;
    mass_max = 6.0;
    
    y_min = -2.4;
    y_max = 2.4;

    pt_min = 0.;
    pt_max = 100.;

    trk1pt_min = 0.;
    trk1pt_max = 30.;

    trk1eta_min = -2.4;
    trk1eta_max = 2.4;

    trk1pterr_min = 0.;  
    trk1pterr_max = 0.4;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 20.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.2;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 1.;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.03;

    mumumass_min = 2.95;
    mumumass_max = 3.25; 
 
    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 0.;
    mu1pt_max = 30.;

    mu2pt_min = 0.;
    mu2pt_max = 30.;

    mu1dxyPV_min = -0.2; 
    mu1dxyPV_max = 0.2;

    mu2dxyPV_min = -0.2;
    mu2dxyPV_max = 0.22;

    mu1dzPV_min = -20.; 
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.; 
    d0_max = 1.;

    d0err_min = 0.;
    d0err_max = 0.0002;

    dtheta_min = 0.;
    dtheta_max = 1.5;

    alpha_min = 0.;
    alpha_max = 3.1;

    trk1Dz1_min = -20.;
    trk1Dz1_max = 20.;

    trk1DzError1_min = -0.1;
    trk1DzError1_max = 0.1;

    trk1Dxy1_min = -0.4;
    trk1Dxy1_max = 0.4;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.03;

    mumueta_min = -3.0;
    mumueta_max = 3.0;

    mumuphi_min = -3.1;
    mumuphi_max = 3.1;

    mumupt_min = 0.;
    mumupt_max = 50.;

    BDT_3_5_min = -0.4;
    BDT_3_5_max = 0.1;

    BDT_5_7_min = -0.3;
    BDT_5_7_max = 0.2;

    BDT_7_10_min = -0.2;
    BDT_7_10_max = 0.3;

    BDT_10_15_min = -0.1;
    BDT_10_15_max = 0.3;

    BDT_15_20_min = -0.09;
    BDT_15_20_max = 0.4;

    BDT_20_50_min = -0.05;
    BDT_20_50_max = 0.4;

    BDT_50_100_min = 0.3;
    BDT_50_100_max = 0.8;

    RooRealVar PVx("PVx","PVx",PVx_min,PVx_max);
    RooRealVar PVy("PVy","PVy",PVy_min,PVy_max);
    RooRealVar PVz("PVz","PVz",PVz_min,PVz_max);
    RooRealVar PVxE("PVxE","PVxE",PVxE_min,PVxE_max);
    RooRealVar PVyE("PVyE","PVyE",PVyE_min,PVyE_max);
    RooRealVar PVzE("PVzE","PVzE",PVzE_min,PVzE_max);
    RooRealVar BvtxX("BvtxX","BvtxX",BvtxX_min,BvtxX_max);
    RooRealVar BvtxY("BvtxY","BvtxY",BvtxY_min,BvtxY_max);
    RooRealVar BvtxZtoPVZ("BvtxZtoPVZ","BvtxZtoPVZ",BvtxZtoPVZ_min,BvtxZtoPVZ_max);
    RooRealVar BSx("BSx","BSx",BSx_min,BSx_max);
    RooRealVar BSy("BSy","BSy",BSy_min,BSy_max);
    RooRealVar BSz("BSz","BSz",BSz_min,BSz_max);
    RooRealVar BSxErr("BSxErr","BSxErr",BSxErr_min,BSxErr_max);
    RooRealVar BSyErr("BSyErr","BSyErr",BSyErr_min,BSyErr_max);
    RooRealVar BSzErr("BSzErr","BSzErr",BSzErr_min,BSzErr_max);
    RooRealVar BSdxdz("BSdxdz","BSdxdz",BSdxdz_min,BSdxdz_max);
    RooRealVar BSdydz("BSdydz","BSdydz",BSdydz_min,BSdydz_max);
    RooRealVar BSdxdzErr("BSdxdzErr","BSdxdzErr",BSdxdzErr_min,BSdxdzErr_max);
    RooRealVar BSdydzErr("BSdydzErr","BSdydzErr",BSdydzErr_min,BSdydzErr_max);
    RooRealVar BSWidthX("BSWidthX","BSWidthX",BSWidthX_min,BSWidthX_max);
    RooRealVar BSWidthXErr("BSWidthXErr","BSWidthXErr",BSWidthXErr_min,BSWidthXErr_max);
    RooRealVar BSWidthY("BSWidthY","BSWidthY",BSWidthY_min,BSWidthY_max);
    RooRealVar BSWidthYErr("BSWidthYErr","BSWidthYErr",BSWidthYErr_min,BSWidthYErr_max);
    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);  
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
    RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDisErr_min, svpvDisErr_max);
    RooRealVar BsvpvDistance_2D("BsvpvDistance_2D", "BsvpvDistance_2D", svpvDistance2D_min, svpvDistance2D_max);
    RooRealVar BsvpvDisErr_2D("BsvpvDisErr_2D", "BsvpvDisErr_2D", svpvDisErr2D_min, svpvDisErr2D_max);
    RooRealVar Bmumumass("Bmumumass", "Bmumumass", mumumass_min, mumumass_max);
    RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
    RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
    RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
    RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
    RooRealVar Bmu1dxyPV("Bmu1dxyPV", "Bmu1dxyPV", mu1dxyPV_min, mu1dxyPV_max);
    RooRealVar Bmu2dxyPV("Bmu2dxyPV", "Bmu2dxyPV", mu1dxyPV_min, mu1dxyPV_max);
    RooRealVar Bmu1dzPV("Bmu1dzPV", "Bmu1dzPV", mu1dzPV_min, mu1dzPV_max);
    RooRealVar Bmu2dzPV("Bmu2dzPV", "Bmu2dzPV", mu1dzPV_min, mu1dzPV_max);
    RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max);
    RooRealVar Bd0Err("Bd0Err", "Bd0Err", d0err_min, d0err_max);
    RooRealVar Bdtheta("Bdtheta", "Bdtheta", dtheta_min, dtheta_max);
    RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
    RooRealVar Btrk1Dz1("Btrk1Dz1","Btrk1Dz1",trk1Dz1_min,trk1Dz1_max);
    RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
    RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
    RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_3_5", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);
    RooRealVar BDT_pt_50_100("BDT_pt_50_100", "BDT_pt_50_100", BDT_50_100_min, BDT_50_100_max);

    w.import(Bmass);
    w.import(PVx);
    w.import(PVy);
    w.import(PVz);
    w.import(PVxE);
    w.import(PVyE);
    w.import(PVzE);
    w.import(BvtxX);
    w.import(BvtxY);
    w.import(BvtxZtoPVZ);
    w.import(BSx);
    w.import(BSy);
    w.import(BSz);
    w.import(BSxErr);
    w.import(BSyErr);
    w.import(BSzErr);
    w.import(BSdxdz);
    w.import(BSdydz);
    w.import(BSdxdzErr);
    w.import(BSdydzErr);
    w.import(BSWidthX);
    w.import(BSWidthXErr);
    w.import(BSWidthY);
    w.import(BSWidthYErr);
    w.import(By);
    w.import(Bpt);
    w.import(Btrk1Pt);
    w.import(Btrk1Eta);
    w.import(Btrk1PtErr);
    w.import(Bchi2cl);
    w.import(BsvpvDistance);
    w.import(BsvpvDisErr);
    w.import(BsvpvDistance_2D);
    w.import(BsvpvDisErr_2D);
    w.import(Bmumumass);
    w.import(Bmu1eta);
    w.import(Bmu2eta);
    w.import(Bmu1pt);
    w.import(Bmu2pt);
    w.import(Bmu1dxyPV);
    w.import(Bmu2dxyPV);
    w.import(Bmu1dzPV);
    w.import(Bmu2dzPV);
    w.import(Bd0);
    w.import(Bd0Err);
    w.import(Bdtheta);
    w.import(Balpha);
    w.import(Btrk1Dz1);
    w.import(Btrk1DzError1);
    w.import(Btrk1Dxy1);
    w.import(Btrk1DxyError1);
    w.import(Bmumueta);
    w.import(Bmumuphi);
    w.import(Bmumupt);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    w.import(BDT_pt_50_100);
  }
      
 else if(particle == 1){
    double PVx_min, PVx_max;
    double PVy_min, PVy_max;
    double PVz_min, PVz_max;
    double PVxE_min, PVxE_max;
    double PVyE_min, PVyE_max;
    double PVzE_min, PVzE_max;
    double BvtxX_min, BvtxX_max;
    double BvtxY_min, BvtxY_max;
    double BvtxZtoPVZ_min, BvtxZtoPVZ_max;
    double BSx_min, BSx_max;
    double BSy_min, BSy_max;
    double BSz_min, BSz_max;
    double BSxErr_min, BSxErr_max;
    double BSyErr_min, BSyErr_max;
    double BSzErr_min, BSzErr_max;
    double BSdxdz_min, BSdxdz_max;
    double BSdydz_min, BSdydz_max;
    double BSdxdzErr_min, BSdxdzErr_max;
    double BSdydzErr_min, BSdydzErr_max;
    double BSWidthX_min, BSWidthX_max;
    double BSWidthXErr_min, BSWidthXErr_max;
    double BSWidthY_min, BSWidthY_max;
    double BSWidthYErr_min, BSWidthYErr_max;
    double mass_min, mass_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double trk1pt_min, trk1pt_max;
    double trk2pt_min, trk2pt_max;
    double trk1eta_min, trk1eta_max;
    double trk2eta_min, trk2eta_max;
    double trk1pterr_min, trk1pterr_max;
    double trk2pterr_min, trk2pterr_max;
    double chi2cl_min, chi2cl_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDisErr_min, svpvDisErr_max;
    double svpvDistance2D_min, svpvDistance2D_max;
    double svpvDisErr2D_min, svpvDisErr2D_max;
    double mumumass_min, mumumass_max;
    double mu1eta_min, mu1eta_max;
    double mu2eta_min, mu2eta_max;
    double mu1pt_min, mu1pt_max;
    double mu2pt_min, mu2pt_max;
    double mu1dxyPV_min, mu1dxyPV_max;
    double mu2dxyPV_min, mu2dxyPV_max;
    double mu1dzPV_min, mu1dzPV_max;
    double mu2dzPV_min, mu2dzPV_max;
    double d0_min, d0_max;
    double d0err_min, d0err_max;
    double dtheta_min, dtheta_max;
    double alpha_min, alpha_max;
    double trk1Dz1_min, trk1Dz1_max;
    double trk2Dz1_min, trk2Dz1_max;
    double trk1DzError1_min, trk1DzError1_max;
    double trk2DzError1_min, trk2DzError1_max;
    double trk1Dxy1_min, trk1Dxy1_max;
    double trk2Dxy1_min, trk2Dxy1_max;
    double trk1DxyError1_min, trk1DxyError1_max;
    double trk2DxyError1_min, trk2DxyError1_max;
    double mumueta_min, mumueta_max;
    double mumuphi_min, mumuphi_max;
    double mumupt_min, mumupt_max;
    double BDT_3_5_min, BDT_3_5_max;
    double BDT_5_7_min, BDT_5_7_max;
    double BDT_7_10_min, BDT_7_10_max;
    double BDT_10_15_min, BDT_10_15_max;
    double BDT_15_20_min, BDT_15_20_max;
    double BDT_20_50_min, BDT_20_50_max;
    double BDT_50_100_min, BDT_50_100_max;

    PVx_min = -0.1;
    PVx_max = 0.1;

    PVy_min = -0.1;
    PVy_max = 0.1;

    PVz_min = 0.;
    PVz_max = 1.5;

    PVxE_min = -0.001;
    PVxE_max = 0.001;

    PVyE_min = -0.001;
    PVyE_max = 0.001;

    PVzE_min = 0.;
    PVzE_max = 0.2;

    BvtxX_min = -4.;
    BvtxX_max = 4.;

    BvtxY_min = -4.;
    BvtxY_max = 4.;

    BvtxZtoPVZ_min = 0.;
    BvtxZtoPVZ_max = 25.;

    BSx_min = -0.1;
    BSx_max = 0.1;

    BSy_min = -0.1;
    BSy_max = 0.1;

    BSz_min = 0.;
    BSz_max = 1.5;

    BSxErr_min = -0.001;
    BSxErr_max = 0.001;

    BSyErr_min = -0.001;
    BSyErr_max = 0.001;

    BSzErr_min = 0.;
    BSzErr_max = 0.2;

    BSdxdz_min = -0.001;
    BSdxdz_max = 0.002;

    BSdydz_min = -0.0004;
    BSdydz_max = 0.0004;

    BSdxdzErr_min = -0.001;
    BSdxdzErr_max = 0.00001;

    BSdydzErr_min = -0.001;
    BSdydzErr_max = 0.00001;

    BSWidthX_min = -0.001;
    BSWidthX_max = 0.005;

    BSWidthXErr_min = -0.001;
    BSWidthXErr_max = 0.009;

    BSWidthY_min = -0.001;
    BSWidthY_max = 0.006;

    BSWidthYErr_min = -0.001;
    BSWidthYErr_max = 0.009;

    mass_min = 5.0;
    mass_max = 6.0;

    y_min = -2.4;
    y_max = 2.4;

    pt_min = 0.;
    pt_max = 100.;

    trk1pt_min = 0.;
    trk1pt_max = 20.;

    trk2pt_min = 0.;
    trk2pt_max = 20.;

    trk1eta_min = -2.4;
    trk1eta_max = 2.4;

    trk2eta_min = -2.4;
    trk2eta_max = 2.4;

    trk1pterr_min = 0.;
    trk1pterr_max = 0.3;

    trk2pterr_min = -0.5;
    trk2pterr_max = 0.5;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 25.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.2;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 1.0;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.035;

    mumumass_min = 2.95;
    mumumass_max = 3.25;

    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 0.;
    mu1pt_max = 70.;

    mu2pt_min = 0.;
    mu2pt_max = 70.;

    mu1dxyPV_min = -0.2;
    mu1dxyPV_max = 0.2;

    mu2dxyPV_min = -0.2;
    mu2dxyPV_max = 0.2;

    mu1dzPV_min = -25.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -25.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 1.;

    d0err_min = 0.;
    d0err_max = 0.0006;

    dtheta_min = 0.;
    dtheta_max = 3.1;

    alpha_min = 0.;
    alpha_max = 3.1;

    trk1Dz1_min = -25.;
    trk1Dz1_max = 20.;

    trk2Dz1_min = -25.;
    trk2Dz1_max = 20.;

    trk1DzError1_min = -0.5;
    trk1DzError1_max = 0.5;

    trk2DzError1_min = -0.5;
    trk2DzError1_max = 0.5;

    trk1Dxy1_min = -1.;
    trk1Dxy1_max = 1.;

    trk2Dxy1_min = -1.;
    trk2Dxy1_max = 1.;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.05;

    trk2DxyError1_min = 0.;
    trk2DxyError1_max = 0.4;

    mumueta_min = -3.1;
    mumueta_max = 3.1;

    mumuphi_min = -3.1;
    mumuphi_max = 3.1;

    mumupt_min = 0.;
    mumupt_max = 100.;

    BDT_3_5_min = -0.8;
    BDT_3_5_max = 0.2;

    BDT_5_7_min = -0.8;
    BDT_5_7_max = 0.2;

    BDT_7_10_min = -0.5;
    BDT_7_10_max = 0.5;

    BDT_10_15_min = -0.3;
    BDT_10_15_max = 0.4;

    BDT_15_20_min = -0.1;
    BDT_15_20_max = 0.5;

    BDT_20_50_min = 0.2;
    BDT_20_50_max = 0.6;

    BDT_50_100_min = 0.1;
    BDT_50_100_max = 0.6;

    RooRealVar PVx("PVx","PVx",PVx_min,PVx_max);
    RooRealVar PVy("PVy","PVy",PVy_min,PVy_max);
    RooRealVar PVz("PVz","PVz",PVz_min,PVz_max);
    RooRealVar PVxE("PVxE","PVxE",PVxE_min,PVxE_max);
    RooRealVar PVyE("PVyE","PVyE",PVyE_min,PVyE_max);
    RooRealVar PVzE("PVzE","PVzE",PVzE_min,PVzE_max);
    RooRealVar BvtxX("BvtxX","BvtxX",BvtxX_min,BvtxX_max);
    RooRealVar BvtxY("BvtxY","BvtxY",BvtxY_min,BvtxY_max);
    RooRealVar BvtxZtoPVZ("BvtxZtoPVZ","BvtxZtoPVZ",BvtxZtoPVZ_min,BvtxZtoPVZ_max);
    RooRealVar BSx("BSx","BSx",BSx_min,BSx_max);
    RooRealVar BSy("BSy","BSy",BSy_min,BSy_max);
    RooRealVar BSz("BSz","BSz",BSz_min,BSz_max);
    RooRealVar BSxErr("BSxErr","BSxErr",BSxErr_min,BSxErr_max);
    RooRealVar BSyErr("BSyErr","BSyErr",BSyErr_min,BSyErr_max);
    RooRealVar BSzErr("BSzErr","BSzErr",BSzErr_min,BSzErr_max);
    RooRealVar BSdxdz("BSdxdz","BSdxdz",BSdxdz_min,BSdxdz_max);
    RooRealVar BSdydz("BSdydz","BSdydz",BSdydz_min,BSdydz_max);
    RooRealVar BSdxdzErr("BSdxdzErr","BSdxdzErr",BSdxdzErr_min,BSdxdzErr_max);
    RooRealVar BSdydzErr("BSdydzErr","BSdydzErr",BSdydzErr_min,BSdydzErr_max);
    RooRealVar BSWidthX("BSWidthX","BSWidthX",BSWidthX_min,BSWidthX_max);
    RooRealVar BSWidthXErr("BSWidthXErr","BSWidthXErr",BSWidthXErr_min,BSWidthXErr_max);
    RooRealVar BSWidthY("BSWidthY","BSWidthY",BSWidthY_min,BSWidthY_max);
    RooRealVar BSWidthYErr("BSWidthYErr","BSWidthYErr",BSWidthYErr_min,BSWidthYErr_max);
    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
    RooRealVar Btrk2Pt("Btrk2Pt","Btrk2Pt",trk2pt_min,trk2pt_max);
    RooRealVar Btrk2Eta("Btrk2Eta","Btrk2Eta",trk2eta_min,trk2eta_max);
    RooRealVar Btrk2PtErr("Btrk2PtErr","Btrk2PtErr",trk2pterr_min,trk2pterr_max);
    RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDisErr_min, svpvDisErr_max);
    RooRealVar BsvpvDistance_2D("BsvpvDistance_2D", "BsvpvDistance_2D", svpvDistance2D_min, svpvDistance2D_max);
    RooRealVar BsvpvDisErr_2D("BsvpvDisErr_2D", "BsvpvDisErr_2D", svpvDisErr2D_min, svpvDisErr2D_max);
    RooRealVar Bmumumass("Bmumumass", "Bmumumass", mumumass_min, mumumass_max);
    RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
    RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
    RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
    RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
    RooRealVar Bmu1dxyPV("Bmu1dxyPV", "Bmu1dxyPV", mu1dxyPV_min, mu1dxyPV_max);
    RooRealVar Bmu2dxyPV("Bmu2dxyPV", "Bmu2dxyPV", mu2dxyPV_min, mu2dxyPV_max);
    RooRealVar Bmu1dzPV("Bmu1dzPV", "Bmu1dzPV", mu1dzPV_min, mu1dzPV_max);
    RooRealVar Bmu2dzPV("Bmu2dzPV", "Bmu2dzPV", mu2dzPV_min, mu2dzPV_max);
    RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max);
    RooRealVar Bd0Err("Bd0Err", "Bd0Err", d0err_min, d0err_max);
    RooRealVar Bdtheta("Bdtheta", "Bdtheta", dtheta_min, dtheta_max);
    RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
    RooRealVar Btrk1Dz1("Btrk1Dz1","Btrk1Dz1",trk1Dz1_min,trk1Dz1_max);
    RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
    RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
    RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
    RooRealVar Btrk2Dz1("Btrk2Dz1","Btrk2Dz1",trk2Dz1_min,trk2Dz1_max);
    RooRealVar Btrk2DzError1("Btrk2DzError1","Btrk2DzError1",trk2DzError1_min,trk2DzError1_max);
    RooRealVar Btrk2Dxy1("Btrk2Dxy1","Btrk2Dxy1",trk2Dxy1_min,trk2Dxy1_max);
    RooRealVar Btrk2DxyError1("Btrk2DxyError1","Btrk2DxyError1",trk2DxyError1_min,trk2DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_3_5", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);
    RooRealVar BDT_pt_50_100("BDT_pt_50_100", "BDT_pt_50_100", BDT_50_100_min, BDT_50_100_max);

    w.import(Bmass);
    w.import(PVx);
    w.import(PVy);
    w.import(PVz);
    w.import(PVxE);
    w.import(PVyE);
    w.import(PVzE);
    w.import(BvtxX);
    w.import(BvtxY);
    w.import(BvtxZtoPVZ);
    w.import(BSx);
    w.import(BSy);
    w.import(BSz);
    w.import(BSxErr);
    w.import(BSyErr);
    w.import(BSzErr);
    w.import(BSdxdz);
    w.import(BSdydz);
    w.import(BSdxdzErr);
    w.import(BSdydzErr);
    w.import(BSWidthX);
    w.import(BSWidthXErr);
    w.import(BSWidthY);
    w.import(BSWidthYErr);
    w.import(By);
    w.import(Bpt);
    w.import(Btrk1Pt);
    w.import(Btrk1Eta);
    w.import(Btrk1PtErr);
    w.import(Bchi2cl);
    w.import(BsvpvDistance);
    w.import(BsvpvDisErr);
    w.import(BsvpvDistance_2D);
    w.import(BsvpvDisErr_2D);
    w.import(Bmumumass);
    w.import(Bmu1eta);
    w.import(Bmu2eta);
    w.import(Bmu1pt);
    w.import(Bmu2pt);
    w.import(Bmu1dxyPV);
    w.import(Bmu2dxyPV);
    w.import(Bmu1dzPV);
    w.import(Bmu2dzPV);
    w.import(Bd0);
    w.import(Bd0Err);
    w.import(Bdtheta);
    w.import(Balpha);
    w.import(Btrk1Dz1);
    w.import(Btrk1DzError1);
    w.import(Btrk1Dxy1);
    w.import(Btrk1DxyError1);
    w.import(Bmumueta);
    w.import(Bmumuphi);
    w.import(Bmumupt);
    w.import(Btrk2Pt);
    w.import(Btrk2Eta);
    w.import(Btrk2PtErr);
    w.import(Btrk2Dz1);
    w.import(Btrk2DzError1);
    w.import(Btrk2Dxy1);
    w.import(Btrk2DxyError1);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    w.import(BDT_pt_50_100);       
  }

  else if(particle == 2){
    double mass_min, mass_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double trk1pt_min, trk1pt_max;
    double trk2pt_min, trk2pt_max;
    double trk1eta_min, trk1eta_max;
    double trk2eta_min, trk2eta_max;
    double trk1pterr_min, trk1pterr_max;
    double trk2pterr_min, trk2pterr_max;
    double chi2cl_min, chi2cl_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDisErr_min, svpvDisErr_max;
    double svpvDistance2D_min, svpvDistance2D_max;
    double svpvDisErr2D_min, svpvDisErr2D_max;
    double mumumass_min, mumumass_max;
    double mu1eta_min, mu1eta_max;
    double mu2eta_min, mu2eta_max;
    double mu1pt_min, mu1pt_max;
    double mu2pt_min, mu2pt_max;
    double mu1dxyPV_min, mu1dxyPV_max;
    double mu2dxyPV_min, mu2dxyPV_max;
    double mu1dzPV_min, mu1dzPV_max;
    double mu2dzPV_min, mu2dzPV_max;
    double d0_min, d0_max;
    double d0err_min, d0err_max;
    double dtheta_min, dtheta_max;
    double alpha_min, alpha_max;
    double trk1Dz1_min, trk1Dz1_max;
    double trk2Dz1_min, trk2Dz1_max;
    double trk1DzError1_min, trk1DzError1_max;
    double trk2DzError1_min, trk2DzError1_max;
    double trk1Dxy1_min, trk1Dxy1_max;
    double trk2Dxy1_min, trk2Dxy1_max;
    double trk1DxyError1_min, trk1DxyError1_max;
    double trk2DxyError1_min, trk2DxyError1_max;
    double mumueta_min, mumueta_max;
    double mumuphi_min, mumuphi_max;
    double mumupt_min, mumupt_max;

    mass_min = 5.0;
    mass_max = 5.6;

    y_min = -2.4;
    y_max = 2.4;

    pt_min = 0.;
    pt_max = 70.;

    trk1pt_min = 0.;
    trk1pt_max = 15.;

    trk2pt_min = 0.;
    trk2pt_max = 15.;

    trk1eta_min = -2.4;
    trk1eta_max = 2.4;

    trk2eta_min = -2.4;
    trk2eta_max = 2.4;

    trk1pterr_min = 0.;
    trk1pterr_max = 0.2;

    trk2pterr_min = 0.;
    trk2pterr_max = 0.2;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 25.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.2;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 0.5;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.02;

    mumumass_min = 2.95;
    mumumass_max = 3.25;

    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 0.;
    mu1pt_max = 40.;

    mu2pt_min = 0.;
    mu2pt_max = 40.;

    mu1dxyPV_min = -0.1;
    mu1dxyPV_max = 0.1;

    mu2dxyPV_min = -0.1;
    mu2dxyPV_max = 0.1;

    mu1dzPV_min = -22.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -22.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 0.5;

    d0err_min = 0.;
    d0err_max = 0.0002;

    dtheta_min = 0.;
    dtheta_max = 3.2;

    alpha_min = 0.;
    alpha_max = 3.1;

    trk1Dz1_min = -22.;
    trk1Dz1_max = 20.;

    trk2Dz1_min = -22.;
    trk2Dz1_max = 20.;

    trk1DzError1_min = 0.;
    trk1DzError1_max = 0.5;

    trk2DzError1_min = 0.;
    trk2DzError1_max = 0.5;

    trk1Dxy1_min = -0.25;
    trk1Dxy1_max = 0.25;

    trk2Dxy1_min = -0.25;
    trk2Dxy1_max = 0.25;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.1;

    trk2DxyError1_min = 0.;
    trk2DxyError1_max = 0.1;

    mumueta_min = -6.;
    mumueta_max = 6.;

    mumuphi_min = -3.1;
    mumuphi_max = 3.1;

    mumupt_min = 0.;
    mumupt_max = 60.;

    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
    RooRealVar Btrk2Pt("Btrk2Pt","Btrk2Pt",trk2pt_min,trk2pt_max);
    RooRealVar Btrk2Eta("Btrk2Eta","Btrk2Eta",trk2eta_min,trk2eta_max);
    RooRealVar Btrk2PtErr("Btrk2PtErr","Btrk2PtErr",trk2pterr_min,trk2pterr_max);
    RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDisErr_min, svpvDisErr_max);
    RooRealVar BsvpvDistance_2D("BsvpvDistance_2D", "BsvpvDistance_2D", svpvDistance2D_min, svpvDistance2D_max);
    RooRealVar BsvpvDisErr_2D("BsvpvDisErr_2D", "BsvpvDisErr_2D", svpvDisErr2D_min, svpvDisErr2D_max);
    RooRealVar Bmumumass("Bmumumass", "Bmumumass", mumumass_min, mumumass_max);
    RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
    RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
    RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
    RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
    RooRealVar Bmu1dxyPV("Bmu1dxyPV", "Bmu1dxyPV", mu1dxyPV_min, mu1dxyPV_max);
    RooRealVar Bmu2dxyPV("Bmu2dxyPV", "Bmu2dxyPV", mu2dxyPV_min, mu2dxyPV_max);
    RooRealVar Bmu1dzPV("Bmu1dzPV", "Bmu1dzPV", mu1dzPV_min, mu1dzPV_max);
    RooRealVar Bmu2dzPV("Bmu2dzPV", "Bmu2dzPV", mu2dzPV_min, mu2dzPV_max);
    RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max);
    RooRealVar Bd0Err("Bd0Err", "Bd0Err", d0err_min, d0err_max);
    RooRealVar Bdtheta("Bdtheta", "Bdtheta", dtheta_min, dtheta_max);
    RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
    RooRealVar Btrk1Dz1("Btrk1Dz1","Btrk1Dz1",trk1Dz1_min,trk1Dz1_max);
    RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
    RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
    RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
    RooRealVar Btrk2Dz1("Btrk2Dz1","Btrk2Dz1",trk2Dz1_min,trk2Dz1_max);
    RooRealVar Btrk2DzError1("Btrk2DzError1","Btrk2DzError1",trk2DzError1_min,trk2DzError1_max);
    RooRealVar Btrk2Dxy1("Btrk2Dxy1","Btrk2Dxy1",trk2Dxy1_min,trk2Dxy1_max);
    RooRealVar Btrk2DxyError1("Btrk2DxyError1","Btrk2DxyError1",trk2DxyError1_min,trk2DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);

    w.import(Bmass);
    w.import(By);
    w.import(Bpt);
    w.import(Btrk1Pt);
    w.import(Btrk1Eta);
    w.import(Btrk1PtErr);
    w.import(Bchi2cl);
    w.import(BsvpvDistance);
    w.import(BsvpvDisErr);
    w.import(BsvpvDistance_2D);
    w.import(BsvpvDisErr_2D);
    w.import(Bmumumass);
    w.import(Bmu1eta);
    w.import(Bmu2eta);
    w.import(Bmu1pt);
    w.import(Bmu2pt);
    w.import(Bmu1dxyPV);
    w.import(Bmu2dxyPV);
    w.import(Bmu1dzPV);
    w.import(Bmu2dzPV);
    w.import(Bd0);
    w.import(Bd0Err);
    w.import(Bdtheta);
    w.import(Balpha);
    w.import(Btrk1Dz1);
    w.import(Btrk1DzError1);
    w.import(Btrk1Dxy1);
    w.import(Btrk1DxyError1);
    w.import(Bmumueta);
    w.import(Bmumuphi);
    w.import(Bmumupt);
    w.import(Btrk2Pt);
    w.import(Btrk2Eta);
    w.import(Btrk2PtErr);
    w.import(Btrk2Dz1);
    w.import(Btrk2DzError1);
    w.import(Btrk2Dxy1);
    w.import(Btrk2DxyError1);
  }
}
