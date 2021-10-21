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

#include <math.h>
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooDoubleCBFast.h"
#include "RooDoubleCBFast.cc"
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
#include <TH1.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooBreitWigner.h>
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
#include <RooArgSet.h>
#include <RooFormulaVar.h>
#include <string>

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
void build_pdf (RooWorkspace& w, std::string choice, RooArgSet &c_vars);
void plot_complete_fit(RooWorkspace& w, RooArgSet &c_vars);
void do_splot(RooWorkspace& w, RooArgSet &c_vars);
TH1D* make_splot(RooWorkspace& w, int n, TString label);
void validate_fit(RooWorkspace* w, RooArgSet &c_vars);
void get_ratio( std::vector<TH1D*>,  std::vector<TH1D*>,  std::vector<TString>, TString);
void pT_analysis(RooWorkspace& w,int n, TString, TString, RooArgSet &c_vars);
double get_yield_syst(RooDataSet *dt, TString syst_src, RooArgSet &c_vars, double pt_min, double pt_max);
void fit_syst_error(TString, RooArgSet &c_vars);
void fit_syst_error_bin(TString, double a, double b, RooArgSet &c_vars);
void constrainVar(TString input_file, TString inVarName, RooArgSet &c_vars, RooArgSet &c_pdfs);
double MC_fit_result(TString input_file, TString inVarName);

//particle
// 0 = Bu
// 1 = Bs
// 2 = B0

#define particle 0

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

//MC
// 1 = does fit to MC 
// 0 = does fit to data

# define MC 0

//component
// 1 = WT
// 0 = RT

# define component 0

void bmesons_new(){
  
  int n_var;
  TString input_file_data;
  if(particle == 0){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/SemiFinalBmesonSamples/BP/BPData.root";}
  else if(particle == 1){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/MoreUpdatedSamples/Bs/BsData.root";}
  else if(particle == 2){input_file_data = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/MoreUpdatedSamples/BZ/BZData.root";}

  TString input_file_mc;
  TString input_file_mc_swap;
  if(particle == 0){input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/ForMariaNewWithTnP/BP/BPMC.root";}
  else if(particle == 1){input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/ForMariaNewWithTnP/Bs/BsMC.root";}
  else if(particle == 2){
    input_file_mc = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/ForMariaNewWithTnP/BZ/BZMC.root";
    input_file_mc_swap = "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/MoreUpdatedSamples/BZ/BZMCSwap2.root";
  }

  TString input_file_reweighted_mc;
  if(particle == 0){input_file_reweighted_mc = "./results/Bu/mc_validation_plots/weights/tree_with_weight.root";}
  else if(particle == 1){input_file_reweighted_mc = "./results/Bs/mc_validation_plots/weights/tree_with_weight.root";}
  else if(particle == 2){input_file_reweighted_mc = "./results/B0/mc_validation_plots/weights/tree_with_weight.root";}

  std::vector<TH1D*> histos_sideband_sub;
  std::vector<TH1D*> histos_mc;
  std::vector<TH1D*> histos_splot;

#if particle == 0
  int n_bins[] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};

  TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr","BsvpvDistance_2D","BsvpvDisErr_2D", "Bmumumass", "Bmu1eta","Bmu2eta", "Bmu1pt", "Bmu2pt","Bmu1dxyPV","Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV","Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "BDT_pt_0_2", "BDT_pt_0_3", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50"};  

#elif particle == 1
  int n_bins[] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};

  TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "BsvpvDistance_2D", "BsvpvDisErr_2D", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV","Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "Btrk2Pt", "Btrk2Eta", "Btrk2PtErr", "Btrk2Dz1", "Btrk2DzError1", "Btrk2Dxy1", "Btrk2DxyError1", "BDT_pt_1_2", "BDT_pt_2_3", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50"};

#elif particle == 2
  int n_bins[] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};

  TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "BsvpvDistance_2D", "BsvpvDisErr_2D", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "Btrk2Pt", "Btrk2Eta", "Btrk2PtErr", "Btrk2Dz1", "Btrk2DzError1", "Btrk2Dxy1", "Btrk2DxyError1", "BDT_pt_0_2", "BDT_pt_2_3", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50"};

#endif

  if( (particle != 2) && (MC == 1) ){return;} // only fits the MC for B0

  int n_n_bins = sizeof(n_bins)/sizeof(n_bins[0]);
  int n_variables = sizeof(variables)/sizeof(variables[0]);

  if(n_n_bins != n_variables){
    std::cout << "Error: number of bins does not correspond to number of variables." << std::endl;
    return;}

  n_var = n_variables;
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);

  if( (MC == 1) && (component == 0) ){read_data(*ws,input_file_mc);}
  else if( (MC == 1) && (component == 1)){read_data(*ws,input_file_mc_swap);}
  else if(MC == 0){read_data(*ws,input_file_data);}

  if(MC == 1){cout << "Running fit on MC" << endl;}
  else if(MC == 0){cout << "Running fit on data" << endl;}

  RooArgSet c_vars;

  build_pdf(*ws, "nominal", c_vars);
  plot_complete_fit(*ws, c_vars);  
  if(MC == 1){return;}
  
  //validate_fit(ws, c_vars);
  //pT_analysis(*ws, n_bins[0], "pT.root", input_file_data, c_vars);     

  // SIDEBAND SUBTRACTION (needs to be run after plot_complete_fit)
  histos_sideband_sub = sideband_subtraction(*ws, n_bins, n_var);
 
  // SPLOT (fixes parameters of the fit -> they need to be unfixed for pT analysis) 
  do_splot(*ws,c_vars); 
  histos_splot = splot_method(*ws, n_bins, variables, n_var);
  
  // MONTE CARLO HISTOGRAMS
  TFile *fin_mc = new TFile(input_file_mc); //use this file to add the weights (to clone original tree) and make data-MC comparisons without weights
  //TFile *fin_mc = new TFile(input_file_reweighted_mc); //use this file to make data-MC comparisons with weights

  TTree* t1_mc;

  if(particle == 0){t1_mc = (TTree*)fin_mc->Get("ntKp");}
  else if(particle == 1){t1_mc = (TTree*)fin_mc->Get("ntphi");}
  else if(particle == 2){t1_mc = (TTree*)fin_mc->Get("ntKstar");}

  std::vector<TString> names;

  for(int i=0; i<n_var; i++){
    TString weight = "weight";
    histos_mc.push_back(create_histogram_mc((*ws->var(variables[i])), t1_mc, n_bins[i], weight));
    names.push_back(TString(variables[i]));
  }
  
  // RATIO BETWEEN DATA (SPLOT) AND MC
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

    if(particle == 0){mc_comp_ss[i]->SetTitle("B^{+}");}
    else if (particle == 1){mc_comp_ss[i]->SetTitle("B^{0}_{s}");}
    else if(particle == 2){mc_comp_ss[i]->SetTitle("B^{0}");}

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
    rp->Draw();
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
      c.SaveAs("./results/Bu/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bu.gif");
    } 
    else if(particle == 1){
      c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.pdf");
      c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.gif");
    }
    else if(particle == 2){
      c.SaveAs("./results/B0/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_B0.pdf");
      c.SaveAs("./results/B0/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_B0.gif");
    } 
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

    if(particle == 0){mc_comp_ss[i]->SetTitle("B^{+}");}
    else if (particle == 1){mc_comp_ss[i]->SetTitle("B^{0}_{s}");}
    else if(particle == 2){mc_comp_ss[i]->SetTitle("B^{0}");}

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
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("Data(ss)/Data(sp)");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    a.Update();
 
    TLegend* leg;
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(ss_comp_sp[i]->GetName(), "S. Subtraction", "l");
    leg->AddEntry(sp_comp_ss[i]->GetName(), "SPlot", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
 
    if(particle == 0){
      a.SaveAs("./results/Bu/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bu.pdf");
      a.SaveAs("./results/Bu/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bu.gif");
    }
    else if(particle == 1){
      a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.pdf");
      a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.gif");
    }
    else if(particle == 2){
      a.SaveAs("./results/B0/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_B0.pdf");
      a.SaveAs("./results/B0/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_B0.gif");
    }
  }
       
  //SPlot vs. Monte Carlo
  vector<TH1D*> sp_comp_mc(histos_splot);
  vector<TH1D*> mc_comp_sp(histos_mc);

  for(int i=0; i<n_var; i++){
    TCanvas b;
    mc_comp_sp[i]->SetXTitle(variables[i]);
    mc_comp_sp[i]->SetYTitle("normalized entries");
    mc_comp_sp[i]->SetStats(0);
    sp_comp_mc[i]->SetStats(0);

    if(particle == 0){mc_comp_ss[i]->SetTitle("B^{+}");}
    else if (particle == 1){mc_comp_ss[i]->SetTitle("B^{0}_{s}");}
    else if(particle == 2){mc_comp_ss[i]->SetTitle("B^{0}");}
 
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
    b.SetTicks(0, 1);
    rp->SetH1DrawOpt("E");
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("Data(sp)/MC");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    b.Update();
     
    TLegend* leg;	
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(mc_comp_sp[i]->GetName(), "Monte Carlo", "l");
    leg->AddEntry(sp_comp_mc[i]->GetName(), "SPlot", "l");
    leg->SetTextSize(0.03);
    leg->Draw("same");
	
    if(particle == 0){
      b.SaveAs("./results/Bu/mc_validation_plots/mc_sp/" + names[i]+"_mc_validation_Bu.pdf");
      b.SaveAs("./results/Bu/mc_validation_plots/mc_sp/" + names[i]+"_mc_validation_Bu.gif");
    }
    else if(particle == 1){
      b.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
      b.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.gif");
    }
    else if(particle == 2){
      b.SaveAs("./results/B0/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_B0.pdf");
      b.SaveAs("./results/B0/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_B0.gif");
    }
  }
 
  //Sideband subtraction vs. Monte Carlo vs SPlot
  vector<TH1D*> sp_comp(histos_splot);
  vector<TH1D*> mc_comp(histos_mc);
  vector<TH1D*> ss_comp(histos_sideband_sub);

  for(int i=0; i<n_var; i++){
    TCanvas d;	
    mc_comp[i]->SetXTitle(variables[i]);
    mc_comp[i]->SetYTitle("normalized entries");
    mc_comp[i]->SetStats(0);
    sp_comp[i]->SetStats(0);
    ss_comp[i]->SetStats(0);

    if(particle == 0){ss_comp[i]->SetTitle("B^{+}");}
    else if (particle == 1){ss_comp[i]->SetTitle("B^{0}_{s}");}
    else if(particle == 2){ss_comp[i]->SetTitle("B^{0}");}

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
      d.SaveAs("./results/Bu/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bu.pdf");
      d.SaveAs("./results/Bu/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bu.gif");
    }
    else if(particle == 1){
      d.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
      d.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.gif");
    }
    else if(particle == 2){
      d.SaveAs("./results/B0/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_B0.pdf");
      d.SaveAs("./results/B0/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_B0.gif");
    }
  }
//comparisons end
}
//main function ends

double MC_fit_result(TString input_file, TString inVarName){

  TFile* f = new TFile(input_file);
  RooFitResult* fitresult;
  fitresult = (RooFitResult*)f->Get("fitresult_model_data");
  RooRealVar* var = (RooRealVar*)fitresult->floatParsFinal().find(inVarName);
 
  return var->getVal();
}

void constrainVar(TString input_file, TString inVarName, RooArgSet &c_vars, RooArgSet &c_pdfs){

  TFile* f = new TFile(input_file);

  RooFitResult* fitresult;

  fitresult = (RooFitResult*)f->Get("fitresult_model_data");

  RooRealVar* var = (RooRealVar*)fitresult->floatParsFinal().find(inVarName);

  RooGaussian* gauss_constr = new RooGaussian(Form("gauss_%s",var->GetName()),
                                              Form("gauss_%s",var->GetName()),
                                              *var,
                                              RooConst(var->getVal()),
                                              RooConst(var->getError())
                                              );
  c_vars.add(*var);
  c_pdfs.add(*gauss_constr);
}


void pT_analysis(RooWorkspace& w, int n, TString ptfile, TString datafile, RooArgSet &c_vars){

  TString dir_name;
  if(particle == 0){dir_name = "./results/Bu/Bpt/";}
  else if(particle == 1){dir_name = "./results/Bs/Bpt/";}
  else if(particle == 2){dir_name = "./results/B0/Bpt/";}

  TFile* f_wei = new TFile(dir_name + ptfile, "recreate"); 

  RooAbsPdf*  model = w.pdf("model");
  RooRealVar* Bpt  = w.var("Bpt");
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooRealVar Bmass = *(w.var("Bmass"));

#if particle == 0
  const int n_pt_bins = 9;
  double pt_bins [n_pt_bins + 1] = {0,2,3,5,7,10,15,20,50,100};  
#elif particle == 1
  const int n_pt_bins = 7;
  double pt_bins[n_pt_bins + 1] = {3,5,7,10,15,20,50,100};
#elif particle == 2
  const int n_pt_bins = 7;
  double pt_bins[n_pt_bins + 1] = {2,3,5,7,10,15,20,50};
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
 
  const int n_pdf_syst=3;
  TString syst_src[n_pdf_syst] = {"nominal","bkg_range","gauss_CB"};
  double yield_syst[n_pt_bins][n_pdf_syst]; 
  double yield_syst_rel[n_pt_bins][n_pdf_syst];

  for(int i=0;i<n_pt_bins;i++){
    //select data subset corresponding to pT bin
    data_pt = (RooDataSet*) data->reduce(Form("Bpt>%lf",pt_bins[i]));
    data_pt = (RooDataSet*) data_pt->reduce(Form("Bpt<%lf",pt_bins[i+1]));
    w.import(*data_pt, Rename(Form("data_pt_%d",i)));
   
    //perform fit and save result
    if((particle == 2) && (MC == 0)){fit_pt = model->fitTo(*data_pt, Minos(true), Save(), Constrain(c_vars));}
    else{fit_pt = model->fitTo(*data_pt, Minos(true), Save());}

    //YIELD + STATISTICAL ERROR
    //floatParsFinal returns the list of floating parameters after fit
    //cout << "Value of floating parameters" << endl;
    //fit_pt->floatParsFinal().Print("s");

    if((particle == 2) && (MC == 0)){n_sig_pt = (RooRealVar*)fit_pt->floatParsFinal().find("RT_yield");}
    else{n_sig_pt = (RooRealVar*)fit_pt->floatParsFinal().find("n_signal");}
    n_comb_pt = (RooRealVar*)fit_pt->floatParsFinal().find("n_combinatorial");

    yield[i]= n_sig_pt->getVal();
    yield_err_low [i] = n_sig_pt->getError(); 
    yield_err_high[i] = n_sig_pt->getError();  

    //cout << "test asym error:" << n_sig_pt->getErrorLo() << " " <<  n_sig_pt->getAsymErrorLo() << " symmetric: " <<  n_sig_pt->getError() <<  endl;

    //SYSTEMATICS
    for(int k = 1; k<n_pdf_syst; k++){
  cout << "PDF = " << syst_src[k] << endl;

      double val = 0.; 
      double val_nominal = 0.;
      val = get_yield_syst(data_pt, syst_src[k], c_vars, pt_bins[i], pt_bins[i+1]); // gets yield value for bin i using pdf choice k
      val_nominal = get_yield_syst(data_pt, syst_src[0], c_vars, pt_bins[i], pt_bins[i+1]); //gets yield value for bin i using nominal pdf choice
      yield_syst_rel[i][k] = (val - val_nominal)/val_nominal;
      yield_syst[i][k] = (val - val_nominal);

      if(k == 0){yield_err_syst[i] += 0;}
      else{yield_err_syst[i] += pow(yield_syst[i][k],2);}
    }
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
 
      RooRealVar* lambda = w.var("lambda");
      RooRealVar* mean  = w.var("mean");
      RooRealVar* sigma1 = w.var("sigma1");

      lambda->setConstant();
      mean->setConstant();
      sigma1->setConstant();
 
      RooRealVar* sigma2;
      RooRealVar* cofs;      

      if((particle == 0) || (particle == 2)){
        sigma2 = w.var("sigma2");
        cofs = w.var("cofs");

        sigma2->setConstant();
        cofs->setConstant();
      }
      
      RooRealVar* alpha1;
      RooRealVar* alpha2;
      RooRealVar* n1;
      RooRealVar* n2;
      RooRealVar* mean_difference;
      RooRealVar* sigma1_swp;
      RooRealVar* alpha1_swp;
      RooRealVar* alpha2_swp;
      RooRealVar* n1_swp;
      RooRealVar* n2_swp;

      if( (particle == 2) && (MC == 0) ){
        alpha1 = w.var("alpha1");
        alpha2 = w.var("alpha2");        
        n1 = w.var("n1");
        n2 = w.var("n2");
        sigma1_swp = w.var("sigma1_swp");
        alpha1_swp = w.var("alpha1_swp");
        alpha2_swp = w.var("alpha2_swp");
        n1_swp = w.var("n1_swp");
        n2_swp = w.var("n2_swp");
        mean_difference = w.var("mean_difference");

        alpha1->setConstant();
        alpha2->setConstant();
        n1->setConstant();
        n2->setConstant();
        sigma1_swp->setConstant();
        alpha1_swp->setConstant();
        alpha2_swp->setConstant();
        n1_swp->setConstant();
        n2_swp->setConstant();
        mean_difference->setConstant();
      }

      SPlot* sData = new SPlot("sData","An sPlot",*data_pt, model, RooArgList(*n_sig_pt,*n_comb_pt));
      
      w.import(*data_pt, Rename(Form("data_pt_WithSWeights_%d",i)));
      
      RooDataSet* data_w = (RooDataSet*) w.data(Form("data_pt_WithSWeights_%d",i));

      RooDataSet* data_wb;
      if( (particle == 2) && (MC == 0) ){data_wb = new RooDataSet(data_w->GetName(),data_w->GetTitle(),data_w,*data_w->get(),0,"RT_yield_sw");}
      else{data_wb = new RooDataSet(data_w->GetName(),data_w->GetTitle(),data_w,*data_w->get(),0,"n_signal_sw");} 
     
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
  f_wei->cd();
  gr->Write();

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
  f_wei->cd();
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
    c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.pdf");
    c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.gif");
  }
  else if(particle == 1){
    c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.pdf");
    c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.gif");
  }
  else if(particle == 2){
    c.SaveAs("./results/B0/Bpt/raw_yield_pt_B0.pdf");
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
    l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.pdf");
    l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.gif");
  }
  else if(particle == 1){
    l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.pdf");
    l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.gif");
  }
  else if(particle == 2){
    l.SaveAs("./results/B0/Bpt/raw_yield_pt_logscale_B0.pdf");
    l.SaveAs("./results/B0/Bpt/raw_yield_pt_logscale_B0.gif");
  }

  // B0 2,3,5,7,10,15,20,50,100
/*
   cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "2-3" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-50" << '|' << endl;
   cout << '|' << setw(15) << "Nominal" << '|' << setw(15) << yield_syst_rel[0][0] << '|' << setw(15) << yield_syst_rel[1][0] << '|' << setw(15) << yield_syst_rel[2][0] << '|' << setw(15) << yield_syst_rel[3][0] << '|' << setw(15) << yield_syst_rel[4][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << setw(15) << yield_syst_rel[6][0] << '|' <<  endl;
   cout << '|' << setw(15) << "Bkg_poly" << '|' << setw(15) << yield_syst_rel[0][1] << '|' << setw(15) << yield_syst_rel[1][1] << '|' << setw(15) << yield_syst_rel[2][1] << '|' << setw(15) << yield_syst_rel[3][1] << '|' << setw(15) << yield_syst_rel[4][1] << '|' << setw(15) << yield_syst_rel[5][1]  << '|' << setw(15) <<  yield_syst_rel[6][1] << '|' << endl;
*/
   // Bs bins 3,5,7,10,15,20,50,100
/*
   cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-40" << '|' << setw(15) << "40-100" << '|' << endl;
   cout << '|' << setw(15) << "Nominal" << '|' << setw(15) << yield_syst_rel[0][0] << '|' << setw(15) << yield_syst_rel[1][0] << '|' << setw(15) << yield_syst_rel[2][0] << '|' << setw(15) << yield_syst_rel[3][0] << '|' << setw(15) << yield_syst_rel[4][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << setw(15) << yield_syst_rel[6][0] << '|' <<  endl;
   cout << '|' << setw(15) << "Bkg_poly" << '|' << setw(15) << yield_syst_rel[0][1]*100 << '|' << setw(15) << yield_syst_rel[1][1]*100 << '|' << setw(15) << yield_syst_rel[2][1]*100 << '|' << setw(15) << yield_syst_rel[3][1]*100 << '|' << setw(15) << yield_syst_rel[4][1]*100 << '|' << setw(15) << yield_syst_rel[5][1]*100 << '|' << setw(15) << yield_syst_rel[6][1]*100 << '|' << endl;
   //cout << '|' << setw(15) << "CB" << '|' << setw(15) << yield_syst_rel[0][2] << '|' << setw(15) << yield_syst_rel[1][2] << '|' << setw(15) << yield_syst_rel[2][2] << '|' << setw(15) << yield_syst_rel[3][2]  << '|' << setw(15) << yield_syst_rel[4][2] << '|' << setw(15) << yield_syst_rel[5][2] << '|' << setw(15) << yield_syst_rel[6][2] << '|' << endl;
   //cout << '|' << setw(15) << "CB" << '|' << setw(15) << yield_syst_rel[0][3] << '|' << setw(15) <<  yield_syst_rel[1][3] << '|' << setw(15) <<  yield_syst_rel[2][3] << '|' << setw(15) <<  yield_syst_rel[3][3] << '|' << setw(15) <<  yield_syst_rel[4][3]  << '|' << endl;
*/

   // 9 bins
    cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "0-2" << '|' << setw(15) << "2-3" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-50" << '|' << setw(15) << "50-100" << '|' << endl;
    cout << '|' << setw(15) << "Nominal" << '|' << setw(15) << yield_syst_rel[0][0]*100 << '|' << setw(15) << yield_syst_rel[1][0]*100 << '|' << setw(15) << yield_syst_rel[2][0]*100 << '|' << setw(15) << yield_syst_rel[3][0]*100 << '|' << setw(15) << yield_syst_rel[4][0]*100 << '|' << setw(15) << yield_syst_rel[5][0]*100 << '|' << setw(15) << yield_syst_rel[6][0]*100 << '|' << setw(15) << yield_syst_rel[7][0]*100 << '|' << setw(15) << yield_syst_rel[8][0]*100 << '|' << endl;
    cout << '|' << setw(15) << "Bkg_range" << '|' << setw(15) << yield_syst_rel[0][1]*100 << '|' << setw(15) << yield_syst_rel[1][1]*100 << '|' << setw(15) << yield_syst_rel[2][1]*100 << '|' << setw(15) << yield_syst_rel[3][1]*100 << '|' << setw(15) << yield_syst_rel[4][1]*100 << '|' << setw(15) << yield_syst_rel[5][1]*100  << '|' << setw(15) << yield_syst_rel[6][1]*100 << '|' << setw(15) << yield_syst_rel[7][1]*100 << '|' << setw(15) << yield_syst_rel[8][1]*100 << '|' << endl;
//    cout << '|' << setw(15) << "Gauss+CB" << '|' << setw(15) << yield_syst_rel[0][2]*100 << '|' << setw(15) << yield_syst_rel[1][2]*100 << '|' << setw(15) << yield_syst_rel[2][2]*100 << '|' << setw(15) << yield_syst_rel[3][2]*100 << '|' << setw(15) << yield_syst_rel[4][2]*100 << '|' << setw(15) << yield_syst_rel[5][2]*100  << '|' << setw(15) << yield_syst_rel[6][2]*100 << '|' << setw(15) << yield_syst_rel[7][2]*100 << '|' << setw(15) << yield_syst_rel[8][2]*100 << '|' << endl;

  // 7 bins
/*
   cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "2-3" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-50" << '|' << endl;
   cout << '|' << setw(15) << "Nominal" << '|' << setw(15) << yield_syst_rel[0][0] << '|' << setw(15) << yield_syst_rel[1][0] << '|' << setw(15) << yield_syst_rel[2][0] << '|' << setw(15) << yield_syst_rel[3][0] << '|' << setw(15) << yield_syst_rel[4][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << setw(15) << yield_syst_rel[6][0]  << '|' << endl;
   cout << '|' << setw(15) << "Bkg_poly" << '|' << setw(15) << yield_syst_rel[0][1] << '|' << setw(15) << yield_syst_rel[1][1] << '|' << setw(15) << yield_syst_rel[2][1] << '|' << setw(15) << yield_syst_rel[3][1] << '|' << setw(15) << yield_syst_rel[4][1] << '|' << setw(15) << yield_syst_rel[5][1]  << '|' << setw(15) << yield_syst_rel[6][1]<< '|' << endl;
   cout << '|' << setw(15) << "Bkg_range" << '|' << setw(15) << yield_syst_rel[0][2] << '|' << setw(15) << yield_syst_rel[1][2] << '|' << setw(15) << yield_syst_rel[2][2] << '|' << setw(15) << yield_syst_rel[3][2] << '|' << setw(15) << yield_syst_rel[4][2] << '|' << setw(15) << yield_syst_rel[5][2] << '|' << setw(15) << yield_syst_rel[6][2] << '|' << endl;
   cout << '|' << setw(15) << "Triple_gauss" << '|' << setw(15) << yield_syst_rel[0][3] << '|' << setw(15) <<  yield_syst_rel[1][3] << '|' << setw(15) <<  yield_syst_rel[2][3] << '|' << setw(15) <<  yield_syst_rel[3][3] << '|' << setw(15) << yield_syst_rel[4][3] << '|' << setw(15) << yield_syst_rel[5][3] << '|' << setw(15) << yield_syst_rel[6][3] << '|' << endl;
*/

// 4 bins
/*
   cout << '|' << setw(15) << "Pdf" << '|' << setw(15) << "2-3" << '|' << setw(15) << "3-5" << '|' << setw(15) << "5-7" << '|' << setw(15) << "7-10" << '|' << setw(15) << "10-15" << '|' << setw(15) << "15-20" << '|' << setw(15) << "20-50" << '|' << endl;
   cout << '|' << setw(15) << syst_src[0] << '|' << setw(15) << yield_syst_rel[0][0] << '|' << setw(15) << yield_syst_rel[1][0] << '|' << setw(15) << yield_syst_rel[2][0] << '|' << setw(15) << yield_syst_rel[3][0]  << '|' << setw(15) << yield_syst_rel[4][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << setw(15) << yield_syst_rel[5][0] << '|' << endl;
   cout << '|' << setw(15) << syst_src[1] << '|' << setw(15) << yield_syst_rel[0][1] << '|' << setw(15) << yield_syst_rel[1][1] << '|' << setw(15) << yield_syst_rel[2][1] << '|' << setw(15) << yield_syst_rel[3][1] << '|' << setw(15) << yield_syst_rel[4][1] << '|' << setw(15) << yield_syst_rel[5][1] << '|' << setw(15) << yield_syst_rel[6][1]  << '|' << endl;
   cout << '|' << setw(15) << syst_src[2] << '|' << setw(15) << yield_syst_rel[0][2] << '|' << setw(15) <<  yield_syst_rel[1][2] << '|' << setw(15) <<  yield_syst_rel[2][2] << '|' << setw(15) <<  yield_syst_rel[3][2] << '|' << setw(15) << yield_syst_rel[4][2] << '|' << setw(15) << yield_syst_rel[5][2] << '|' << setw(15) << yield_syst_rel[6][2] << '|' << endl;
*/
}
//pT_analysis end

//get the ratio between the data (splot method) and the MC and save it in a root file
void get_ratio( std::vector<TH1D*> data, std::vector<TH1D*> mc,  std::vector<TString> v_name, TString filename){

  TString dir_name;
  if(particle == 0){dir_name = "./results/Bu/mc_validation_plots/weights/";}
  else if(particle == 1){dir_name = "./results/Bs/mc_validation_plots/weights/";}
  else if(particle == 2){dir_name = "./results/B0/mc_validation_plots/weights/";}

  TFile* f_wei = new TFile(dir_name + filename, "recreate");

  TH1D* h_aux;
  h_aux->SetDefaultSumw2(kTRUE);

  for(int i=0; i<(int)data.size(); i++) {

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
    c.SaveAs(dir_name+v_name.at(i) + "_weights.gif");
    //output: a root file and plots gifs
  
  }

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
    arg_list.add(*(w.var("BDT_pt_0_2")));
    arg_list.add(*(w.var("BDT_pt_0_3")));
    arg_list.add(*(w.var("BDT_pt_3_5")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
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
    arg_list.add(*(w.var("BDT_pt_1_2")));
    arg_list.add(*(w.var("BDT_pt_2_3")));
    arg_list.add(*(w.var("BDT_pt_3_5")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
  }
  if(particle == 2){
    arg_list.add(*(w.var("BDT_pt_0_2")));
    arg_list.add(*(w.var("BDT_pt_2_3")));
    arg_list.add(*(w.var("BDT_pt_3_5")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
  }
  
  RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
  w.import(*data, Rename("data"));
}

void build_pdf(RooWorkspace& w, std::string choice, RooArgSet &c_vars){

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*)w.data("data");

  TString input_file_RT = "./results/B0/MC/RT_fit.root";
  TString input_file_WT = "./results/B0/MC/WT_fit.root";

  RooArgSet c_pdfs_RT;
  RooArgSet c_pdfs_WT;

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
    left = 5.22;
    right = 5.32;
    mass_peak = 5.280;
  }

  if( (particle == 2) && (MC == 0) && (choice != "scale_factor") ){ 
    cout << "Applying constraints" << endl;
    constrainVar(input_file_RT, "sigma1", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "sigma2", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "alpha1", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "alpha2", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "n1", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "n2", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "cofs", c_vars, c_pdfs_RT);

    constrainVar(input_file_WT, "sigma1_swp", c_vars, c_pdfs_WT);
    constrainVar(input_file_WT, "alpha1_swp", c_vars, c_pdfs_WT);
    constrainVar(input_file_WT, "alpha2_swp", c_vars, c_pdfs_WT);
    constrainVar(input_file_WT, "n1_swp", c_vars, c_pdfs_WT);
    constrainVar(input_file_WT, "n2_swp", c_vars, c_pdfs_WT);
  }

  // Variable initialization
  RooRealVar* f_swap = 0;
  RooRealVar* mean = 0;
  RooRealVar* sigma1 = 0;
  RooRealVar* sigma2 = 0;
  RooRealVar* alpha1 = 0;
  RooRealVar* alpha1_swp = 0;
  RooRealVar* n1 = 0;
  RooRealVar* n1_swp = 0;
  RooRealVar* cofs = 0;
  RooRealVar* mean_swp = 0;
  RooRealVar* sigma1_swp = 0;
  RooRealVar* sigma2_swp = 0;
  RooRealVar* cofs_swp = 0;
  RooRealVar* cofs1_swp = 0;
  RooRealVar* sigma3 = 0;
  RooRealVar* cofs1 = 0;
  RooRealVar* m_nonprompt_scale = 0;
  RooRealVar* m_nonprompt_shift = 0;
  RooRealVar* lambda = 0;
  RooRealVar* slope = 0;
  RooRealVar* m_jpsipi_mean1 = 0;
  RooRealVar* m_jpsipi_sigma1l = 0;
  RooRealVar* m_jpsipi_sigma1r = 0;
  RooRealVar* m_jpsipi_mean2 = 0;
  RooRealVar* m_jpsipi_mean3 = 0;
  RooRealVar* m_jpsipi_sigma2 = 0;
  RooRealVar* m_jpsipi_sigma3 = 0;
  RooRealVar* m_jpsipi_fraction2 = 0;
  RooRealVar* m_jpsipi_fraction3 = 0;
  RooRealVar* alpha2 = 0;
  RooRealVar* alpha2_swp = 0;
  RooRealVar* n2 = 0;
  RooRealVar* n2_swp = 0;
  RooProduct* sigma1_fix = 0;
  RooProduct* sigma2_fix = 0;
  RooProduct* sigma1_swp_fix = 0;
  RooRealVar* scale_factor = new RooRealVar("scale_factor", "scale_factor", 1., 0., 2.);
  RooRealVar* mean_difference = 0;

  if( (particle == 2) && (MC == 0)){
    cout << "Initialisating variables with constraints" << endl;

    // mis-tag fraction
    double n_rt = MC_fit_result(input_file_RT, "n_signal");
    double n_wt = MC_fit_result(input_file_WT, "n_signal_swp");
    double fraction = n_wt / (n_wt + n_rt);
    cout << "fraction = " << fraction << endl;
    f_swap = new RooRealVar("f_swap","f_swap", fraction, 0., 1.);       
    f_swap->setConstant();
    cout << "f_swap = " << f_swap->getVal() << endl;  

    // RT component (CB + CB)
    mean = new RooRealVar("mean", "mean", MC_fit_result(input_file_RT, "mean"), MC_fit_result(input_file_RT, "mean")-0.1, MC_fit_result(input_file_RT, "mean")+0.1);
    sigma1 = new RooRealVar("sigma1", "sigma1", MC_fit_result(input_file_RT, "sigma1"), 0.005, 0.5);
    sigma2 = new RooRealVar("sigma2", "sigma2", MC_fit_result(input_file_RT, "sigma2"), 0.005 ,0.5);
    alpha1 = new RooRealVar("alpha1", "alpha1", MC_fit_result(input_file_RT, "alpha1"), 0., 20.);
    alpha2 = new RooRealVar("alpha2", "alpha2", MC_fit_result(input_file_RT, "alpha2"), 0., 20.);
    n1 = new RooRealVar("n1", "n1", MC_fit_result(input_file_RT, "n1"), 0., 300.);
    n2 = new RooRealVar("n2", "n2", MC_fit_result(input_file_RT, "n2"), 0., 300.);
    cofs = new RooRealVar("cofs", "cofs", MC_fit_result(input_file_RT, "cofs"), 0., 1.);

    // WT component (double CB)
    mean_swp = new RooRealVar("mean_swp", "mean_swp", MC_fit_result(input_file_WT, "mean_swp"), MC_fit_result(input_file_WT, "mean_swp")-0.1, MC_fit_result(input_file_WT, "mean_swp")+0.1);
    sigma1_swp = new RooRealVar("sigma1_swp", "sigma1_swp", MC_fit_result(input_file_WT, "sigma1_swp"), 0.005, 0.5);
    alpha1_swp = new RooRealVar("alpha1_swp", "alpha1_swp", MC_fit_result(input_file_WT, "alpha1_swp"), 0., 20.);
    alpha2_swp = new RooRealVar("alpha2_swp", "alpha2_swp", MC_fit_result(input_file_WT, "alpha2_swp"), 0., 20.);
    n1_swp = new RooRealVar("n1_swp", "n1_swp", MC_fit_result(input_file_WT, "n1_swp"), 0., 300.);
    n2_swp = new RooRealVar("n2_swp", "n2_swp", MC_fit_result(input_file_WT, "n2_swp"), 0., 300.);

    if(choice == "scale_factor"){
      sigma1->setConstant();
      sigma2->setConstant();
      alpha1->setConstant();
      alpha2->setConstant();
      n1->setConstant();
      n2->setConstant();
      cofs->setConstant();
      sigma1_swp->setConstant();
      alpha1_swp->setConstant();
      alpha2_swp->setConstant();
      n1_swp->setConstant();
      n2_swp->setConstant();
  
      sigma1_fix = new RooProduct("sigma1_fix", "sigma1_fix", RooArgList(*scale_factor,*sigma1));
      sigma2_fix = new RooProduct("sigma2_fix", "sigma2_fix", RooArgList(*scale_factor,*sigma2));
      sigma1_swp_fix = new RooProduct("sigma1_swp_fix", "sigma1_swp_fix", RooArgList(*scale_factor,*sigma1_swp));
    }

  }
  else{
    cout << "Initialising variables without constraints" << endl;
    mean = new RooRealVar("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.2);
    sigma1 = new RooRealVar("sigma1","sigma1",0.02,0.005,0.5);
    sigma2 = new RooRealVar("sigma2","sigma2",0.01,0.005,0.5);
    alpha1 = new RooRealVar("alpha1", "alpha1", 5., 0., 20.); 
    alpha1_swp = new RooRealVar("alpha1_swp", "alpha1_swp", 5., 0., 20.);
    alpha2_swp = new RooRealVar("alpha2_swp", "alpha2_swp", 5., 0., 20.);
    n1 = new RooRealVar("n1", "n1", 10., 0., 300.);
    n1_swp = new RooRealVar("n1_swp", "n1_swp", 10., 0., 300.);
    n2_swp = new RooRealVar("n2_swp", "n2_swp", 10., 0., 300.);
    cofs = new RooRealVar("cofs", "cofs", 0.3, 0., 1.);
    mean_swp = new RooRealVar("mean_swp","mean_swp",mass_peak,mass_peak-0.1,mass_peak+0.2);
    sigma1_swp = new RooRealVar("sigma1_swp","sigma1_swp",0.02,0.005,0.5);
    sigma2_swp = new RooRealVar("sigma2_swp","sigma2_swp",0.01,0.005,0.5);
    cofs_swp = new RooRealVar("cofs_swp", "cofs_swp", 0.3, 0., 1.);
    cofs1_swp = new RooRealVar("cofs1_swp", "cofs1_swp", 0.3, 0., 1.);
    f_swap = new RooRealVar("f_swap","f_swap",0.,0.,1.);
  }

  sigma3 = new RooRealVar("sigma3","sigma3",0.012,0.010,0.030);
  cofs1 = new RooRealVar("cofs1", "cofs1", 0.3, 0., 1.);
  m_nonprompt_scale = new RooRealVar("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
  m_nonprompt_shift = new RooRealVar("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  lambda = new RooRealVar("lambda","lambda",-2.,-5.,1.0);
  slope = new RooRealVar("slope","slope",0,-5,5);
  m_jpsipi_mean1 = new RooRealVar("m_jpsipi_mean1","m_jpsipi_mean1",5.34693, 5.346, 5.347);
  m_jpsipi_sigma1l = new RooRealVar("m_jpsipi_sigma1l","m_jpsipi_sigma1l",0.0290762,0.010,0.150);
  m_jpsipi_sigma1r = new RooRealVar("m_jpsipi_sigma1r","m_jpsipi_sigma1r",0.0652519,0.010,0.150);
  m_jpsipi_mean2 = new RooRealVar("m_jpsipi_mean2","m_jpsipi_mean2",5.46876);
  m_jpsipi_mean3 = new RooRealVar("m_jpsipi_mean3","m_jpsipi_mean3",5.48073);
  m_jpsipi_sigma2 = new RooRealVar("m_jpsipi_sigma2","m_jpsipi_sigma2",0.0994712,0.020,0.500);
  m_jpsipi_sigma3 = new RooRealVar("m_jpsipi_sigma3","m_jpsipi_sigma3",0.330152,0.020,0.500);
  m_jpsipi_fraction2 = new RooRealVar("m_jpsipi_fraction2","m_jpsipi_fraction2",0.234646,0.0,1.0);
  m_jpsipi_fraction3 = new RooRealVar("m_jpsipi_fraction3","m_jpsipi_fraction3",0.114338,0.0,1.0);
  alpha2 = new RooRealVar("alpha2", "alpha2", 5., -20., 20.);
  n2 = new RooRealVar("n2", "n2", 5., 0., 300.);

  //SIGNAL
  //sum of two gaussians
  RooGaussian* signal1 = new RooGaussian("signal1","signal_gauss1",Bmass,*mean,*sigma1);
  RooGaussian* signal2 = new RooGaussian("signal2","signal_gauss2",Bmass,*mean,*sigma2);
  RooAddPdf* signal = new RooAddPdf("signal", "signal", RooArgList(*signal1,*signal2),*cofs);

  // triple gaussian
  RooGaussian* signal3 = new RooGaussian("signal3","signal3",Bmass,*mean,*sigma3);  
  RooAddPdf* signal_triple = new RooAddPdf("signal_triple","signal_triple",RooArgList(*signal1,*signal2,*signal3),RooArgList(*cofs,*cofs1));

  // crystal ball function
  RooCBShape* CB1 = new RooCBShape("CB1","CB1",Bmass,*mean,*sigma1,*alpha1,*n1);
  RooCBShape* CB2 = new RooCBShape("CB2", "CB2",Bmass,*mean,*sigma2,*alpha2,*n2);
  RooAddPdf* sum_CB = new RooAddPdf("sum_CB","sum_CB",RooArgList(*CB1,*CB2),*cofs);

  RooAddPdf* gauss_CB = new RooAddPdf("gauss_CB","gauss_CB",RooArgList(*signal1,*CB1),*cofs);
  RooAddPdf* two_gauss_CB = new RooAddPdf("two_gauss_CB","two_gauss_CB",RooArgList(*signal,*CB1),*cofs);

  // WT component
  RooDoubleCBFast* double_CB_swp = new RooDoubleCBFast("double_CB_swp", "double_CB_swp", Bmass, *mean_swp, *sigma1_swp, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);

  //BACKGROUND//
  //error function (for JPsi X peaking background)
  m_nonprompt_shift->setConstant(kTRUE);
  m_nonprompt_scale->setConstant(kTRUE);
  RooGenericPdf* erf = new RooGenericPdf("erf","erf","TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)",RooArgList(Bmass,*m_nonprompt_scale,*m_nonprompt_shift));
 
  //exponential (for combinatorial background)
  RooExponential* fit_side = new RooExponential("fit_side","fit_side",Bmass,*lambda);

  // 1st order polynomial (combinatorial background - pdf systematics)
  RooPolynomial* poly_bkg = new RooPolynomial("poly_bkg","poly_bkg",Bmass,*slope);

  //jpsi_pi component (for jpsi background)
  m_jpsipi_mean1->setConstant(kTRUE);
  m_jpsipi_mean2->setConstant(kTRUE);
  m_jpsipi_mean3->setConstant(kTRUE);
  m_jpsipi_sigma1l->setConstant(kTRUE);
  m_jpsipi_sigma1r->setConstant(kTRUE);
  m_jpsipi_sigma2->setConstant(kTRUE);
  m_jpsipi_sigma3->setConstant(kTRUE);
  m_jpsipi_fraction2->setConstant(kTRUE);
  m_jpsipi_fraction3->setConstant(kTRUE);  

  RooBifurGauss* m_jpsipi_gaussian1 = new RooBifurGauss("m_jpsipi_gaussian1","m_jpsipi_gaussian1",Bmass,*m_jpsipi_mean1,*m_jpsipi_sigma1l,*m_jpsipi_sigma1r);
  RooGaussian* m_jpsipi_gaussian2 = new RooGaussian("m_jpsipi_gaussian2","m_jpsipi_gaussian2",Bmass,*m_jpsipi_mean2,*m_jpsipi_sigma2);
  RooGaussian* m_jpsipi_gaussian3 = new RooGaussian("m_jpsipi_gaussian3","m_jpsipi_gaussian3",Bmass,*m_jpsipi_mean3,*m_jpsipi_sigma3);
  RooAddPdf* jpsipi = new RooAddPdf("jpsipi","jpsipi",RooArgList(*m_jpsipi_gaussian3,*m_jpsipi_gaussian2,*m_jpsipi_gaussian1),RooArgList(*m_jpsipi_fraction3,*m_jpsipi_fraction2));

  Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);
  Bmass.setRange("peakright",left,Bmass.getMax());

  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  double WT_yield_initial = 0;
  double RT_yield_initial = 0;
  if((particle == 2) && (MC == 0)){
    WT_yield_initial = (f_swap->getVal())*n_signal_initial;
    RT_yield_initial = n_signal_initial-WT_yield_initial;
  }
 cout << "here 0 " << endl;
  // B0 FULL MODEL
  // RT component
  RooCBShape* cb1_rt = new RooCBShape("cb1_rt","cb1_rt",Bmass,*mean,*sigma1,*alpha1,*n1);
  RooCBShape* cb2_rt = new RooCBShape("cb2_rt","cb2_rt",Bmass,*mean,*sigma2,*alpha2,*n2);
  RooAddPdf* sum_cb_rt = new RooAddPdf("sum_cb_rt","sum_cb_rt",RooArgList(*cb1_rt,*cb2_rt),*cofs);

  RooRealVar* RT_yield = new RooRealVar("RT_yield", "RT_yield", RT_yield_initial, 0.,  data->sumEntries());

  RooArgList constr_rt_list = RooArgList(c_pdfs_RT);
  constr_rt_list.add(*sum_cb_rt);
  RooProdPdf* rt_pdf = new RooProdPdf("rt_pdf", "rt_pdf", constr_rt_list); 
cout << "here 1" << endl;
  // WT component
  mean_difference = new RooRealVar("mean_difference", "mean_difference", MC_fit_result(input_file_WT, "mean_swp") - MC_fit_result(input_file_RT, "mean"), -2, 2, "GeV");
  RooProduct* mean_difference_fix = 0;
  if((particle == 2) && (MC == 0) && (choice == "scale_factor")){
    mean_difference->setConstant();
    mean_difference_fix = new RooProduct("mean_difference_fix", "mean_difference_fix", RooArgList(*scale_factor,*mean_difference));
  }
cout << "here 2" << endl;
  RooFormulaVar* mass_swp;
  if((particle == 2) && (MC == 0) && (choice == "scale_factor")){
    mass_swp = new RooFormulaVar("mass_swp", "mass_swp", "@0+@1", RooArgList(*mean,*mean_difference_fix));
  }
  else if ((particle == 2) && (MC == 0) && (choice != "scale_factor")){mass_swp = new RooFormulaVar("mass_swp", "mass_swp", "@0+@1", RooArgList(*mean,*mean_difference));}

  RooDoubleCBFast* double_CB_wt;
  if((particle == 2)){double_CB_wt = new RooDoubleCBFast("double_CB_wt", "double_CB_wt", Bmass, *mass_swp, *sigma1_swp, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);}

  RooProduct* WT_yield = new RooProduct("WT_yield","WT_yield",RooArgList(*f_swap,*RT_yield));
cout << "here 3" << endl;
  // mean difference constraint 
  TFile* file_RT = new TFile(input_file_RT);
  TFile* file_WT = new TFile(input_file_WT);

  RooFitResult* fr = (RooFitResult*)file_RT->Get("fitresult_model_data");
  RooFitResult* fr_swp = (RooFitResult*)file_WT->Get("fitresult_model_data");
 
  RooRealVar* constrained_mean =  (RooRealVar*)fr->floatParsFinal().find("mean");
  RooRealVar* constrained_mean_swp =  (RooRealVar*)fr_swp->floatParsFinal().find("mean_swp");

  double mean_diff_val = (constrained_mean_swp->getVal()) - (constrained_mean->getVal());
  double mean_diff_error = sqrt( pow(constrained_mean->getError(),2) + pow(constrained_mean_swp->getError(),2) );

  RooGaussian* mean_constr = new RooGaussian("mean_constr", "mean_constr", *mean_difference, RooConst(mean_diff_val), RooConst(mean_diff_error));
  c_vars.add(*mean_difference);

  RooArgList constr_wt_list = RooArgList(c_pdfs_WT);
  constr_wt_list.add(*mean_constr);
  if(particle == 2){constr_wt_list.add(*double_CB_wt);}
  RooProdPdf* wt_pdf = new RooProdPdf("wt_pdf", "wt_pdf", constr_wt_list);
cout << "here 4" << endl;
  // Scale factor
  RooCBShape* cb1_rt_sf;
  RooCBShape* cb2_rt_sf; 
  RooAddPdf* sum_cb_rt_sf; 
  RooDoubleCBFast* double_CB_wt_sf;

  if((particle == 2) && (MC == 0) && (choice == "scale_factor")){
    cb1_rt_sf = new RooCBShape("cb1_rt_sf","cb1_rt_sf",Bmass,*mean,*sigma1_fix,*alpha1,*n1);
    cb2_rt_sf = new RooCBShape("cb2_rt_sf","cb2_rt_sf",Bmass,*mean,*sigma2_fix,*alpha2,*n2);
    sum_cb_rt_sf = new RooAddPdf("sum_cb_rt_sf","sum_cb_rt_sf",RooArgList(*cb1_rt_sf,*cb2_rt_sf),*cofs);
    double_CB_wt_sf = new RooDoubleCBFast("double_CB_wt_sf", "double_CB_wt_sf", Bmass, *mass_swp, *sigma1_swp_fix, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);
  }
cout << "here 5" << endl;
  // NORMALISATIONS
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar* n_signal = new RooRealVar("n_signal","n_signal",n_signal_initial,0.,(data->sumEntries())*2);
  RooRealVar* n_signal_swp = new RooRealVar("n_signal_swp","n_signal_swp",n_signal_initial,0.,(data->sumEntries())*2);
  RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
  RooRealVar* f_erf = new RooRealVar("f_erf","f_erf",0.2,0.,1.);
  RooProduct* n_erf = new RooProduct("n_erf","n_erf",RooArgList(*n_signal,*f_erf));
  RooRealVar* f_jpsipi = new RooRealVar("f_jpsipi","f_jpsipi",0.03996, 0.038, 0.040);
  f_jpsipi->setConstant(kTRUE);
  RooProduct* n_jpsipi = new RooProduct("n_jpsipi","n_jpsipi",RooArgList(*n_signal,*f_jpsipi)); 

  if(particle == 0){//B+
    if(choice == "nominal"){
      RooAddPdf model("model","model",RooArgList(*signal,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "bkg_poly"){
      RooAddPdf model("model","model",RooArgList(*signal,*poly_bkg,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    } else if(choice == "bkg_range"){
      RooAddPdf model("model","model",RooArgList(*signal,*fit_side,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_jpsipi));
      w.import(model);
    }else if(choice == "signal1gauss"){
      RooAddPdf model("model","model",RooArgList(*signal3,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "crystal_ball"){
      RooAddPdf model("model","model",RooArgList(*CB1,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "gauss_CB"){
      RooAddPdf model("model","model",RooArgList(*gauss_CB,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "triple_gauss"){
      RooAddPdf model("model","model",RooArgList(*signal_triple,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }
  }

  else if(particle == 1){//Bs 
    if(choice == "nominal"){
      RooAddPdf model("model","model",RooArgList(*signal1,*fit_side),RooArgList(*n_signal,*n_combinatorial)); 
      w.import(model);
    }else if(choice == "bkg_poly"){
      RooAddPdf model("model","model",RooArgList(*signal1,*poly_bkg),RooArgList(*n_signal,*n_combinatorial));
      w.import(model);
    }else if(choice == "signal1gauss"){
      RooAddPdf model("model","model",RooArgList(*signal3,*fit_side),RooArgList(*n_signal,*n_combinatorial));
      w.import(model);
    }else if(choice == "triple_gauss"){
      RooAddPdf model("model","model",RooArgList(*signal_triple,*fit_side),RooArgList(*n_signal,*n_combinatorial));
      w.import(model);
    }else if(choice == "crystal_ball"){
      RooAddPdf model("model","model",RooArgList(*CB1,*fit_side),RooArgList(*n_signal,*n_combinatorial));
      w.import(model);
    }
  }

  else if(particle == 2){
    if(MC == 0){
      if(choice == "nominal"){
        RooAddPdf model("model","model",RooArgList(*rt_pdf,*wt_pdf,*fit_side),RooArgList(*RT_yield,*WT_yield,*n_combinatorial));
        w.import(model);
      }else if(choice == "bkg_poly"){
        RooAddPdf model("model","model",RooArgList(*rt_pdf,*wt_pdf,*poly_bkg),RooArgList(*RT_yield,*WT_yield,*n_combinatorial));
        w.import(model);
      }else if(choice == "scale_factor"){
        RooAddPdf model("model","model",RooArgList(*sum_cb_rt_sf,*double_CB_wt_sf,*fit_side),RooArgList(*RT_yield,*WT_yield,*n_combinatorial));
        w.import(model);
      }
    }
    else if(MC == 1){
      if(component == 0){
        RooAddPdf model("model","model",RooArgList(*sum_CB),RooArgList(*n_signal));
        w.import(model);
      }
      else if(component == 1){
        RooAddPdf model("model", "model",RooArgList(*double_CB_swp),RooArgList(*n_signal_swp));
        w.import(model);
      }
    }
  }
} 
//build_pdf ends

void fit_syst_error(TString fname, RooArgSet &c_vars){
  //returns the yield's systematic uncertainty
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,fname);

  RooDataSet* data = (RooDataSet*) ws->data("data");

  build_pdf(*ws,"nominal", c_vars);
  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres_nom = model->fitTo(*data,Save());

  build_pdf(*ws,"bkg_poly", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_bgpol = model->fitTo(*data,Save());

  build_pdf(*ws,"bkg_range", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_bgrange = model->fitTo(*data,Range("peakright"),Save());

  build_pdf(*ws,"signal1gauss", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_sig1 = model->fitTo(*data,Save());
 
  build_pdf(*ws,"triple_gauss", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_sig2 = model->fitTo(*data,Save());
 
  build_pdf(*ws,"crystal_ball", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_sig3 = model->fitTo(*data,Save());

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

void fit_syst_error_bin(TString fname, double bin_min, double bin_max, RooArgSet &c_vars){
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

  build_pdf(*ws, "nominal", c_vars);
  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres_nom = model->fitTo(*data_bin,Save());

  build_pdf(*ws,"bkg_poly", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_bgpol = model->fitTo(*data_bin,Save());

  build_pdf(*ws,"bkg_range", c_vars);
  model = ws->pdf("model");
  RooFitResult* fitres_bgrange = model->fitTo(*data_bin,Range("peakright"),Save());

  build_pdf(*ws,"sig1gauss", c_vars);
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

double get_yield_syst(RooDataSet* data_bin, TString syst_src, RooArgSet &c_vars, double pt_min, double pt_max) {
  //returns the yield's value per bin
  
  //cout << "aaa 0\n";
  //data_bin->Print();
  //cout << "aaa 1\n";

  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  ws->import(*data_bin, Rename("data"));
  //  ws->Import(data_bin);

  RooRealVar* Bmass = ws->var("Bmass");

  //  TString mdl = (syst_src.Contains("range")) ? "nominal" : syst_src;
  TString rng = (syst_src.Contains("range")) ? "peakright" : "all";
  cout << "pdf = " << syst_src.Data() << endl;
  build_pdf(*ws,syst_src.Data(), c_vars);

  RooAbsPdf* model = ws->pdf("model");
  RooFitResult* fitres;
  if((particle == 2) && (MC == 0)){fitres = model->fitTo(*data_bin,Range(rng),Save(),Constrain(c_vars));}
  else{fitres = model->fitTo(*data_bin,Range(rng),Save());}

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

  RooPlot* massframe = Bmass->frame(Title(""));
  data_bin->plotOn(massframe,RooFit::Name("Data"));
  model->paramOn(massframe,Layout(0.60,0.99,0.95));
  model->plotOn(massframe, RooFit::Name("Fit"),Range("all"));

  massframe->Draw();

  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  RooPlot* pull_plot = Bmass->frame(Title(""));

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
    b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_Bu.gif");
    b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_Bu.pdf");
  }
  else if(particle == 1){
    b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_Bs.gif");
    b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_Bs.pdf");
  }
  else if(particle == 2){
    b.SaveAs(Form("./results/B0/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_B0.gif");
    b.SaveAs(Form("./results/B0/Bpt/%lf_%lf", pt_min, pt_max)+syst_src+"_fit_plot_B0.pdf");   
  }

  RooRealVar* n1_var = 0;
  if( (particle == 2) && (MC == 0) ){n1_var = (RooRealVar*) fitres ->floatParsFinal().find("RT_yield");}
  else{n1_var = (RooRealVar*) fitres ->floatParsFinal().find("n_signal");}

  double n1 = n1_var->getVal();

  return n1; 
}
//get_yield_syst ends

void plot_complete_fit(RooWorkspace& w, RooArgSet &c_vars){

  RooAbsPdf*  model = w.pdf("model");
  RooDataSet* data = (RooDataSet*) w.data("data");
  data->Print();

  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* lambda   = w.var("lambda");

  if( (particle == 2) && (MC == 0) ){model->fitTo(*data,Range("all"),Constrain(c_vars),Extended(kTRUE));}
  else{model->fitTo(*data,Range("all"));}

  TFile* f;
  if(MC == 1){
    if(component == 0){f = new TFile("./results/B0/MC/RT_fit.root", "RECREATE");}
    else if(component == 1){f = new TFile("./results/B0/MC/WT_fit.root", "RECREATE");}
  }
  else if(MC == 0){f = new TFile("./results/B0/MC/DATA_fit.root", "RECREATE");}

  RooFitResult* r = model->fitTo(*data,Range("all"),Save());
  r->Print();
  f->cd();
  r->Write();
  f->Close();
 
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
  }
  else if(particle == 1){
    data->plotOn(massframe, RooFit::Name("Data"));
    model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal1"),Range("all"),LineColor(kOrange),LineStyle(kDashed));
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("Bmass (GeV)");
  }  
  else if(particle == 2){
      data->plotOn(massframe, RooFit::Name("Data"));
      model->plotOn(massframe, RooFit::Name("Fit"), Range("all"), LineColor(kMagenta), LineStyle(1), LineWidth(2));
      if(MC == 0){
        model->plotOn(massframe, RooFit::Name("Corr Tag"), RooFit::Components("rt_pdf"), RooFit::Range("all"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
        model->plotOn(massframe, RooFit::Name("Mis Tag"), RooFit::Components("wt_pdf"), RooFit::Range("all"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));  
        model->plotOn(massframe, RooFit::Name("Combinatorial"), RooFit::Components("fit_side"), RooFit::Range("all"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
      }
  }
  model->paramOn(massframe,Layout(0.55,0.95,0.90));

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

  /*
  TLatex* tex11 = new TLatex(0.6,0.8,"302.3 pb^{-1} (pp) 5.02 TeV");
  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  tex11->Draw();
  tex11 = new TLatex(0.6,0.85,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  tex11->Draw();
  
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
  */

  TLegend *leg = new TLegend (0.7, 0.7, 0.9, 0.9);

  if(particle == 0){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/psi X", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/psi pi", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  }
  else if(particle == 1){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  }
  else if(particle == 2){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
    if(MC == 0){
      leg->AddEntry(massframe->findObject("Fit"), "Total sig", "l");
      leg->AddEntry(massframe->findObject("Corr Tag"), "Corr. sig", "l");
      leg->AddEntry(massframe->findObject("Mis Tag"), "Mis-tag sig", "l");
      leg->AddEntry(massframe->findObject("Combinatorial"), "Comb. bkg", "l");
    }
  }
  //leg->Draw("same");

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
    d.SaveAs("./results/Bu/complete_fit_Bu.pdf");
    d.SaveAs("./results/Bu/complete_fit_Bu.gif");
  }
  else if(particle == 1){
    d.SaveAs("./results/Bs/complete_fit_Bs.pdf");
    d.SaveAs("./results/Bs/complete_fit_Bs.gif");
  }
  else if(particle == 2){
    if(MC == 0){
      d.SaveAs("./results/B0/DATA_fit_B0.gif");
      d.SaveAs("./results/B0/DATA_fit_B0.pdf");
    }
    else if(MC == 1){
      if(component == 0){
        d.SaveAs("./results/B0/MC_RT_fit_B0.gif");
        d.SaveAs("./results/B0/MC_RT_fit_B0.pdf");
      }
      else if(component == 1){
        d.SaveAs("./results/B0/MC_WT_fit_B0.gif");
        d.SaveAs("./results/B0/MC_WT_fit_B0.pdf");
      }
    }
  }
}
//plot_complete_fit ends

//SIDEBAND SUBTRACTION//
std::vector<TH1D*> sideband_subtraction(RooWorkspace w, int* n, int n_var){
  
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooAbsPdf* BpModel;
  if(particle ==0 ){BpModel =  w.pdf("signal");}
  else if(particle == 1){BpModel =  w.pdf("signal1");}
  else if(particle == 2){BpModel = w.pdf("rt_pdf");}
  RooAbsPdf* BgModel = w.pdf("fit_side");

  RooRealVar Bmass = *(w.var("Bmass"));

  vector<RooRealVar> variables;

  variables.push_back(*(w.var("Bmass"))); 
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
    variables.push_back(*(w.var("BDT_pt_0_2")));
    variables.push_back(*(w.var("BDT_pt_0_3")));
    variables.push_back(*(w.var("BDT_pt_3_5")));
    variables.push_back(*(w.var("BDT_pt_5_7")));
    variables.push_back(*(w.var("BDT_pt_7_10")));
    variables.push_back(*(w.var("BDT_pt_10_15")));
    variables.push_back(*(w.var("BDT_pt_15_20")));
    variables.push_back(*(w.var("BDT_pt_20_50")));
  }
  if((particle == 1) || (particle == 2)){
    variables.push_back(*(w.var("Btrk2Pt")));
    variables.push_back(*(w.var("Btrk2Eta")));
    variables.push_back(*(w.var("Btrk2PtErr")));
    variables.push_back(*(w.var("Btrk2Dz1")));
    variables.push_back(*(w.var("Btrk2DzError1")));
    variables.push_back(*(w.var("Btrk2Dxy1")));
    variables.push_back(*(w.var("Btrk2DxyError1")));
  }
  if(particle == 1){
    variables.push_back(*(w.var("BDT_pt_1_2")));
    variables.push_back(*(w.var("BDT_pt_2_3")));
    variables.push_back(*(w.var("BDT_pt_3_5")));
    variables.push_back(*(w.var("BDT_pt_5_7")));
    variables.push_back(*(w.var("BDT_pt_7_10")));
    variables.push_back(*(w.var("BDT_pt_10_15")));
    variables.push_back(*(w.var("BDT_pt_15_20")));
    variables.push_back(*(w.var("BDT_pt_20_50")));
  }
  if(particle == 2){
    variables.push_back(*(w.var("BDT_pt_0_2")));
    variables.push_back(*(w.var("BDT_pt_2_3")));
    variables.push_back(*(w.var("BDT_pt_3_5")));
    variables.push_back(*(w.var("BDT_pt_5_7")));
    variables.push_back(*(w.var("BDT_pt_7_10")));
    variables.push_back(*(w.var("BDT_pt_10_15")));
    variables.push_back(*(w.var("BDT_pt_15_20")));
    variables.push_back(*(w.var("BDT_pt_20_50")));
  } 
  
  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central; 

  double left;
  double right;

  if(particle == 0){
    left = 5.2;
    right = 5.4;
  }
  else if(particle == 1){
    left = 5.3;
    right = 5.45;
  }
  else if(particle == 2){
    left = 5.22;
    right = 5.32;
  }

  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);
  Bmass.setRange("peakright",left,Bmass.getMax());
  Bmass.setRange("total", Bmass.getMin(), Bmass.getMax());
  
  if(particle == 0){reduceddata_side = (RooDataSet*)data->reduce(Form("Bmass>%lf", right));}
  else if( (particle == 1) || (particle == 2) ){reduceddata_side =  (RooDataSet*)data->reduce(Form("Bmass>%lf || Bmass<%lf", right, left));}

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

  double factor;
  if(particle == 0){factor = (int_fit_peak->getVal())/(int_fit_side_right->getVal());}
  else if( (particle == 1) || (particle == 2) ){(factor = int_fit_peak->getVal())/(int_fit_side_right->getVal() + int_fit_side_left->getVal());}
  std::cout << std::endl << "Factor: " << factor << std::endl;

  for(int i=0; i<n_var; i++){
    std::cout << "bins: " << n[i] << std::endl;
  } 
  std::vector<TH1D*> histos;

  if(particle == 0){
    histos.push_back(create_histogram(variables[1], "By",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "BsvpvDistance_2D",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "BsvpvDisErr_2D",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[29]));
    histos.push_back(create_histogram(variables[31], "BDT_pt_0_2",factor, reduceddata_side, reduceddata_central, data, n[30]));
    histos.push_back(create_histogram(variables[32], "BDT_pt_0_3",factor, reduceddata_side, reduceddata_central, data, n[31]));
    histos.push_back(create_histogram(variables[33], "BDT_pt_3_5",factor, reduceddata_side, reduceddata_central, data, n[32]));
    histos.push_back(create_histogram(variables[34], "BDT_pt_5_7",factor, reduceddata_side, reduceddata_central, data, n[33]));
    histos.push_back(create_histogram(variables[35], "BDT_pt_7_10",factor, reduceddata_side, reduceddata_central, data, n[34]));
    histos.push_back(create_histogram(variables[36], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[35]));
    histos.push_back(create_histogram(variables[37], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[36]));
    histos.push_back(create_histogram(variables[38], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[37]));
  }else if(particle == 1){
    histos.push_back(create_histogram(variables[1], "By",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "BsvpvDistance_2D",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "BsvpvDisErr_2D",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[29]));
    histos.push_back(create_histogram(variables[31], "Btrk2Pt",factor, reduceddata_side, reduceddata_central, data, n[30]));
    histos.push_back(create_histogram(variables[32], "Btrk2Eta",factor, reduceddata_side, reduceddata_central, data, n[31]));
    histos.push_back(create_histogram(variables[33], "Btrk2PtErr",factor, reduceddata_side, reduceddata_central, data, n[32]));
    histos.push_back(create_histogram(variables[34], "Btrk2Dz1",factor, reduceddata_side, reduceddata_central, data, n[33]));
    histos.push_back(create_histogram(variables[35], "Btrk2DzError1",factor, reduceddata_side, reduceddata_central, data, n[34]));
    histos.push_back(create_histogram(variables[36], "Btrk2Dxy1",factor, reduceddata_side, reduceddata_central, data, n[35]));
    histos.push_back(create_histogram(variables[37], "Btrk2DxyError1",factor, reduceddata_side, reduceddata_central, data, n[36]));
    histos.push_back(create_histogram(variables[38], "BDT_pt_1_2",factor, reduceddata_side, reduceddata_central, data, n[37]));
    histos.push_back(create_histogram(variables[39], "BDT_pt_2_3",factor, reduceddata_side, reduceddata_central, data, n[38]));
    histos.push_back(create_histogram(variables[40], "BDT_pt_3_5",factor, reduceddata_side, reduceddata_central, data, n[39]));
    histos.push_back(create_histogram(variables[41], "BDT_pt_5_7",factor, reduceddata_side, reduceddata_central, data, n[40]));
    histos.push_back(create_histogram(variables[42], "BDT_pt_7_10",factor, reduceddata_side, reduceddata_central, data, n[41]));
    histos.push_back(create_histogram(variables[43], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[42]));
    histos.push_back(create_histogram(variables[44], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[43]));
    histos.push_back(create_histogram(variables[45], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[44]));
  }else if(particle == 2){
    histos.push_back(create_histogram(variables[1], "By",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "BsvpvDistance_2D",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "BsvpvDisErr_2D",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[29]));
    histos.push_back(create_histogram(variables[31], "Btrk2Pt",factor, reduceddata_side, reduceddata_central, data, n[30]));
    histos.push_back(create_histogram(variables[32], "Btrk2Eta",factor, reduceddata_side, reduceddata_central, data, n[31]));
    histos.push_back(create_histogram(variables[33], "Btrk2PtErr",factor, reduceddata_side, reduceddata_central, data, n[32]));
    histos.push_back(create_histogram(variables[34], "Btrk2Dz1",factor, reduceddata_side, reduceddata_central, data, n[33]));
    histos.push_back(create_histogram(variables[35], "Btrk2DzError1",factor, reduceddata_side, reduceddata_central, data, n[34]));
    histos.push_back(create_histogram(variables[36], "Btrk2Dxy1",factor, reduceddata_side, reduceddata_central, data, n[35]));
    histos.push_back(create_histogram(variables[37], "Btrk2DxyError1",factor, reduceddata_side, reduceddata_central, data, n[36]));
    histos.push_back(create_histogram(variables[38], "BDT_pt_0_2",factor, reduceddata_side, reduceddata_central, data, n[37]));
    histos.push_back(create_histogram(variables[39], "BDT_pt_2_3",factor, reduceddata_side, reduceddata_central, data, n[37]));
    histos.push_back(create_histogram(variables[40], "BDT_pt_3_5",factor, reduceddata_side, reduceddata_central, data, n[38]));
    histos.push_back(create_histogram(variables[41], "BDT_pt_5_7",factor, reduceddata_side, reduceddata_central, data, n[40]));
    histos.push_back(create_histogram(variables[42], "BDT_pt_7_10",factor, reduceddata_side, reduceddata_central, data, n[41]));
    histos.push_back(create_histogram(variables[43], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[42]));
    histos.push_back(create_histogram(variables[44], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[43]));
    histos.push_back(create_histogram(variables[45], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[44]));
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
  }else if(particle == 2){
    c.SaveAs("./results/B0/sideband_sub/"+name + "sideband_sub_B0.pdf");
    c.SaveAs("./results/B0/sideband_sub/"+name + "sideband_sub_B0.gif");
  }

  if(background == 0){return dist_peak;}
  else if(background == 1){return dist_side;}

}
//create_histogram ends

void do_splot(RooWorkspace& w, RooArgSet &c_vars){

  RooDataSet* data = (RooDataSet*) w.data("data");   
  RooAbsPdf* model = w.pdf("model");
  //we need the fit and the dataset previously saved in the woorkspace

  RooRealVar* BpYield = 0;
  if( (particle == 2) && (MC == 0) ){BpYield = w.var("RT_yield");}
  else{BpYield = w.var("n_signal");}
  RooRealVar* BgYield = w.var("n_combinatorial");
  //we need the n values previously saved in the woorkspace

  //fit the model to the data
  if( (particle == 2) && (MC == 0) ){model->fitTo(*data,Extended(),Constrain(c_vars));}
  else{model->fitTo(*data,Extended());}

  //sPlot technique requires model parameters (other than the yields) to be fixed
  RooRealVar* lambda = w.var("lambda");
  RooRealVar* mean  = w.var("mean");
  RooRealVar* sigma1 = w.var("sigma1");

  lambda->setConstant();
  mean->setConstant();
  sigma1->setConstant();

  RooRealVar* sigma2;
  RooRealVar* cofs;

  if((particle == 0) || (particle == 2)){
    sigma2 = w.var("sigma2");
    cofs = w.var("cofs");

    sigma2->setConstant();
    cofs->setConstant();
  }

  RooRealVar* alpha1;
  RooRealVar* alpha2;
  RooRealVar* n1;
  RooRealVar* n2;
  RooRealVar* mean_difference;
  RooRealVar* sigma1_swp;
  RooRealVar* alpha1_swp;
  RooRealVar* alpha2_swp;
  RooRealVar* n1_swp;
  RooRealVar* n2_swp;

  if( (particle == 2) && (MC == 0) ){
    alpha1 = w.var("alpha1");
    alpha2 = w.var("alpha2");
    n1 = w.var("n1");
    n2 = w.var("n2");
    mean_difference = w.var("mean_difference");
    sigma1_swp = w.var("sigma1_swp");
    alpha1_swp = w.var("alpha1_swp");
    alpha2_swp = w.var("alpha2_swp");
    n1_swp = w.var("n1_swp");
    n2_swp = w.var("n2_swp");

    alpha1->setConstant();
    alpha2->setConstant();
    n1->setConstant();
    n2->setConstant();
    mean_difference->setConstant();
    sigma1_swp->setConstant();
    alpha1_swp->setConstant();
    alpha2_swp->setConstant();
    n1_swp->setConstant();
    n2_swp->setConstant();
  }

  RooMsgService::instance().setSilentMode(true);

  //add sWeights to dataset based on model and yield variables
  //sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
  SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*BpYield,*BgYield));

  cout << endl <<  "Yield of B+ is "
       << BpYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("RT_yield") << endl;

  cout << "Yield of background is "
       << BgYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;
  
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
  RooAbsPdf* BpModel;

  if( (particle == 2) && (MC == 0) ){BpModel = w.pdf("rt_pdf");}
  else if(particle == 0){BpModel = w.pdf("signal");}
  else if(particle == 1){BpModel = w.pdf("signal1");}
  RooAbsPdf* BgModel = w.pdf("fit_side");

  RooRealVar* Bmass  = w.var("Bmass");
  RooRealVar* variable = w.var(label);

  RooRealVar* BpYield = 0;
  if( ( particle == 2) && (MC == 0) ){BpYield = w.var("RT_yield");}
  else{ BpYield = w.var("n_signal");}
  RooRealVar* BgYield = w.var("n_combinatorial");

  double sigYield = BpYield->getVal();
  double bkgYield  = BgYield->getVal();

  RooDataSet* data = (RooDataSet*)w.data("data");

  cdata->cd(1);
  RooPlot* mframe = Bmass->frame();
  if(particle == 0){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of B+ [GeV]"));
  }else if(particle == 1){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of Bs [GeV]"));
  }else if(particle == 2){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of B0 [GeV]"));
  }

  data->plotOn(mframe);
  model->plotOn(mframe,LineColor(kRed));
  model->plotOn(mframe,Components(*BpModel),LineStyle(kDashed),LineColor(kOrange));
  model->plotOn(mframe,Components(*BgModel),LineStyle(kDashed),LineColor(kBlue));
  model->paramOn(mframe,Layout(0.60,0.99,0.95));
  mframe->SetTitle("Bmass");
  mframe->Draw();

  cdata->cd(2);
  RooPlot* ptframe = variable->frame();
  data->plotOn(ptframe);
  if(particle == 0){
    ptframe->SetTitle(label + " of B+: total sample");
  }else if(particle == 1){
    ptframe->SetTitle(label + " of Bs: total sample");
  }else if(particle == 2){
    ptframe->SetTitle(label + " of B0: total sample");
  }
  ptframe->Draw();

  //get the dataset with sWeights
  RooDataSet* dataW = (RooDataSet*) w.data("dataWithSWeights");
  RooDataSet* dataWBp;
  if( (particle == 2) && (MC == 0) ){dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"RT_yield_sw");}
  else{dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");}
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
  }else if(particle == 2){
    ptframe2Bp->GetXaxis()->SetTitle(label + " of B0");
    ptframe2Bg->GetXaxis()->SetTitle(label + " of B0");
  }

  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));
  dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(n));

  if(particle == 0){
    ptframe2Bp->SetTitle(label+" distribution of B+ for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of B+ for background (splot)");
  }else if(particle == 1){
    ptframe2Bp->SetTitle(label+" distribution of Bs for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of Bs for background (splot)");
  }else if(particle == 2){
    ptframe2Bp->SetTitle(label+" distribution of B0 for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of B0 for background (splot)");
  }

  cdata->cd(3);  ptframe2Bp->Draw();
  cdata->cd(4);  ptframe2Bg->Draw();

  if(particle == 0){
    cdata->SaveAs("./results/Bu/splot/Bmass/"+label+"sPlot_Bu.gif");
    cdata->SaveAs("./results/Bu/splot/Bmass/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    cdata->SaveAs("./results/Bs/splot/Bmass/"+label+"sPlot_Bs.gif");
    cdata->SaveAs("./results/Bs/splot/Bmass/"+label+"sPlot_Bs.pdf");
  }else if(particle == 2){
    cdata->SaveAs("./results/B0/splot/Bmass/"+label+"sPlot_B0.gif");
    cdata->SaveAs("./results/B0/splot/Bmass/"+label+"sPlot_B0.pdf");
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
  }else if(particle == 2){
    prov->SaveAs("./results/B0/splot/sig/"+label+"sPlot_B0.gif");
    prov->SaveAs("./results/B0/splot/sig/"+label+"sPlot_B0.pdf");
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
  }else if(particle == 2){
    prov_bkg->SaveAs("./results/B0/splot/bkg/"+label+"sPlot_B0.gif");
    prov_bkg->SaveAs("./results/B0/splot/bkg/"+label+"sPlot_B0.pdf");
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
  }else if(particle == 2){
    sig_bkg->SaveAs("./results/B0/splot/sig_bkg/"+label+"sPlot_B0.gif");
    sig_bkg->SaveAs("./results/B0/splot/sig_bkg/"+label+"sPlot_B0.pdf");
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

void validate_fit(RooWorkspace* w, RooArgSet &c_vars)
{
  RooRealVar Bmass = *(w->var("Bmass"));
  RooAbsPdf* model  = w->pdf("model");
  
  vector<RooRealVar> params;
  if( (particle == 2) && (MC == 0) ){params.push_back(*(w->var("RT_yield")));}
  else{params.push_back(*(w->var("n_signal")));}

  int params_size = params.size();  

  RooMCStudy* mcstudy;
  if( (particle == 2) && (MC == 0) ){mcstudy = new RooMCStudy(*model, Bmass, Constrain(c_vars), Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));}
  else{mcstudy = new RooMCStudy(*model, Bmass, Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));}

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
    c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.pdf");
    c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.gif");
    c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.pdf");
    c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.gif");
  }else if(particle == 1){
    c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.pdf");
    c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.gif");
    c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.pdf");
    c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.gif");
  }else if(particle == 2){
    c_pull->SaveAs("./results/B0/pulls/pulls_poisson_B0.pdf");
    c_pull->SaveAs("./results/B0/pulls/pulls_poisson_B0.gif");
    c_params->SaveAs("./results/B0/pulls/pulls_params_poisson_B0.pdf");
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
    float mass_min, mass_max;
    float y_min, y_max;
    float pt_min, pt_max;
    float trk1pt_min, trk1pt_max;
    float trk1eta_min, trk1eta_max;
    float trk1pterr_min, trk1pterr_max;
    float chi2cl_min, chi2cl_max;
    float svpvDistance_min, svpvDistance_max;
    float svpvDisErr_min, svpvDisErr_max;
    float svpvDistance2D_min, svpvDistance2D_max;
    float svpvDisErr2D_min, svpvDisErr2D_max;
    float mumumass_min, mumumass_max; 
    float mu1eta_min, mu1eta_max;
    float mu2eta_min, mu2eta_max;
    float mu1pt_min, mu1pt_max;
    float mu2pt_min, mu2pt_max;
    float mu1dxyPV_min, mu1dxyPV_max;
    float mu2dxyPV_min, mu2dxyPV_max;
    float mu1dzPV_min, mu1dzPV_max;
    float mu2dzPV_min, mu2dzPV_max;
    float d0_min, d0_max;
    float d0err_min, d0err_max;
    float dtheta_min, dtheta_max;
    float alpha_min, alpha_max;
    float trk1Dz1_min, trk1Dz1_max;
    float trk1DzError1_min, trk1DzError1_max;
    float trk1Dxy1_min, trk1Dxy1_max;
    float trk1DxyError1_min, trk1DxyError1_max;
    float mumueta_min, mumueta_max;
    float mumuphi_min, mumuphi_max;
    float mumupt_min, mumupt_max;
    float BDT_0_2_min, BDT_0_2_max;
    float BDT_0_3_min, BDT_0_3_max;
    float BDT_3_5_min, BDT_3_5_max;
    float BDT_5_7_min, BDT_5_7_max;
    float BDT_7_10_min, BDT_7_10_max;
    float BDT_10_15_min, BDT_10_15_max;
    float BDT_15_20_min, BDT_15_20_max;
    float BDT_20_50_min, BDT_20_50_max;

    mass_min = 5.;
    mass_max = 6.;
    
    y_min = -2.4;
    y_max = 2.4;

    pt_min = 0.;
    pt_max = 100.;

    trk1pt_min = 0.;
    trk1pt_max = 60.;

    trk1eta_min = -6.;
    trk1eta_max = 6.;

    trk1pterr_min = 0.;  
    trk1pterr_max = 2;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 25.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.4;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 2.5;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.04;

    mumumass_min = 2.97;
    mumumass_max = 3.22; 
 
    mu1eta_min = -6.;
    mu1eta_max = 6.;

    mu2eta_min = -6.;
    mu2eta_max = 6.;

    mu1pt_min = 0.;
    mu1pt_max = 80.;

    mu2pt_min = 0.;
    mu2pt_max = 80.;

    mu1dxyPV_min = -1; 
    mu1dxyPV_max = 1;

    mu2dxyPV_min = -0.5;
    mu2dxyPV_max = 0.5;

    mu1dzPV_min = -20.; 
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.; 
    d0_max = 2;

    d0err_min = 0.;
    d0err_max = 0.0005;

    dtheta_min = 0.;
    dtheta_max = 3.2;

    alpha_min = 0.;
    alpha_max = 3.2;

    trk1Dz1_min = -100.;
    trk1Dz1_max = 100.;

    trk1DzError1_min = -0.5;
    trk1DzError1_max = 0.5;

    trk1Dxy1_min = -50;
    trk1Dxy1_max = 50;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.5;

    mumueta_min = -6.;
    mumueta_max = 6.;

    mumuphi_min = -3.2;
    mumuphi_max = 3.2;

    mumupt_min = 0.;
    mumupt_max = 100.;

    BDT_0_2_min = -0.3;
    BDT_0_2_max = 0.2;

    BDT_0_3_min = -0.9;
    BDT_0_3_max = 0.3;

    BDT_3_5_min = -0.4;
    BDT_3_5_max = 0.4;

    BDT_5_7_min = -0.3;
    BDT_5_7_max = 0.5;

    BDT_7_10_min = -0.2;
    BDT_7_10_max = 0.4;

    BDT_10_15_min = -0.2;
    BDT_10_15_max = 0.4;

    BDT_15_20_min = -0.2;
    BDT_15_20_max = 0.4;

    BDT_20_50_min = -0.2;
    BDT_20_50_max = 0.5;

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
    RooRealVar BDT_pt_0_2("BDT_pt_0_2", "BDT_pt_0_2", BDT_0_2_min, BDT_0_2_max);
    RooRealVar BDT_pt_0_3("BDT_pt_0_3", "BDT_pt_0_3", BDT_0_3_min, BDT_0_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_3_5", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);

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
    w.import(BDT_pt_0_2);
    w.import(BDT_pt_0_3);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
 }
      
 else if(particle == 1){
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
    double BDT_1_2_min, BDT_1_2_max;
    double BDT_2_3_min, BDT_2_3_max;
    double BDT_3_5_min, BDT_3_5_max;
    double BDT_5_7_min, BDT_5_7_max;
    double BDT_7_10_min, BDT_7_10_max;
    double BDT_10_15_min, BDT_10_15_max;
    double BDT_15_20_min, BDT_15_20_max;
    double BDT_20_50_min, BDT_20_50_max;

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
    trk1pterr_max = 0.4;

    trk2pterr_min = 0.;
    trk2pterr_max = 0.4;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 25.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.5;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 1.0;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.05;

    mumumass_min = 2.95;
    mumumass_max = 3.25;

    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 0.;
    mu1pt_max = 80.;

    mu2pt_min = 0.;
    mu2pt_max = 80.;

    mu1dxyPV_min = -0.2;
    mu1dxyPV_max = 0.2;

    mu2dxyPV_min = -0.2;
    mu2dxyPV_max = 0.2;

    mu1dzPV_min = -20.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 1.;

    d0err_min = 0.;
    d0err_max = 0.0005;

    dtheta_min = 0.;
    dtheta_max = 3.2;

    alpha_min = 0.;
    alpha_max = 3.2;

    trk1Dz1_min = -20.;
    trk1Dz1_max = 20.;

    trk2Dz1_min = -20.;
    trk2Dz1_max = 20.;

    trk1DzError1_min = -0.5;
    trk1DzError1_max = 0.5;

    trk2DzError1_min = -0.5;
    trk2DzError1_max = 0.5;

    trk1Dxy1_min = -0.5;
    trk1Dxy1_max = 0.5;

    trk2Dxy1_min = -0.5;
    trk2Dxy1_max = 0.5;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.2;

    trk2DxyError1_min = 0.;
    trk2DxyError1_max = 0.2;

    mumueta_min = -4.5;
    mumueta_max = 4.5;

    mumuphi_min = -3.2;
    mumuphi_max = 3.2;

    mumupt_min = 0.;
    mumupt_max = 100.;

    BDT_1_2_min = -0.6;
    BDT_1_2_max = 0.5;

    BDT_2_3_min = -0.5;
    BDT_2_3_max = 0.5;

    BDT_3_5_min = -0.6;
    BDT_3_5_max = 0.3;

    BDT_5_7_min = -0.6;
    BDT_5_7_max = 0.5;

    BDT_7_10_min = -0.45;
    BDT_7_10_max = 0.45;

    BDT_10_15_min = -0.5;
    BDT_10_15_max = 0.55;

    BDT_15_20_min = -0.5;
    BDT_15_20_max = 0.5;

    BDT_20_50_min = -0.2;
    BDT_20_50_max = 0.4;

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
    RooRealVar BDT_pt_1_2("BDT_pt_1_2", "BDT_pt_1_2", BDT_1_2_min, BDT_1_2_max);
    RooRealVar BDT_pt_2_3("BDT_pt_2_3", "BDT_pt_2_3", BDT_2_3_min, BDT_2_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_3_5", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);

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
    w.import(BDT_pt_1_2);
    w.import(BDT_pt_2_3);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
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
    float BDT_0_2_min, BDT_0_2_max;
    float BDT_2_3_min, BDT_2_3_max;
    float BDT_3_5_min, BDT_3_5_max;
    float BDT_5_7_min, BDT_5_7_max;
    float BDT_7_10_min, BDT_7_10_max;
    float BDT_10_15_min, BDT_10_15_max;
    float BDT_15_20_min, BDT_15_20_max;
    float BDT_20_50_min, BDT_20_50_max;

    if((MC == 1) && (component == 1)){
      mass_min = 4.9;
      mass_max = 5.7;
    }
    else{
      mass_min = 5.0;
      mass_max = 5.6;
    }

    y_min = -2.4;
    y_max = 2.4;

    pt_min = 0.;
    pt_max = 50.;

    trk1pt_min = 0.;
    trk1pt_max = 20.;

    trk2pt_min = 0.;
    trk2pt_max = 20.;

    trk1eta_min = -2.5;
    trk1eta_max = 2.5;

    trk2eta_min = -2.5;
    trk2eta_max = 2.5;

    trk1pterr_min = 0.;
    trk1pterr_max = 0.5;

    trk2pterr_min = 0.;
    trk2pterr_max = 0.5;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 20.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.1;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 1.;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.02;

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

    mu1dxyPV_min = -0.25;
    mu1dxyPV_max = 0.2;

    mu2dxyPV_min = -0.25;
    mu2dxyPV_max = 0.2;

    mu1dzPV_min = -20.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 1.;

    d0err_min = 0.;
    d0err_max = 0.001;

    dtheta_min = 0.;
    dtheta_max = 3.2;

    alpha_min = 0.;
    alpha_max = 3.2;

    trk1Dz1_min = -50.;
    trk1Dz1_max = 50.;

    trk2Dz1_min = -50.;
    trk2Dz1_max = 50.;

    trk1DzError1_min = -0.5;
    trk1DzError1_max = 0.5;

    trk2DzError1_min = -0.5;
    trk2DzError1_max = 0.5;

    trk1Dxy1_min = -20;
    trk1Dxy1_max = 20;

    trk2Dxy1_min = -20;
    trk2Dxy1_max = 20;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.2;

    trk2DxyError1_min = 0.;
    trk2DxyError1_max = 0.2;

    mumueta_min = -6.;
    mumueta_max = 6.;

    mumuphi_min = -3.2;
    mumuphi_max = 3.2;

    mumupt_min = 0.;
    mumupt_max = 45.;

    BDT_0_2_min = -0.8;
    BDT_0_2_max = 0.3;

    BDT_2_3_min = -0.6;
    BDT_2_3_max = 0.4;

    BDT_3_5_min = -0.5;
    BDT_3_5_max = 0.5;

    BDT_5_7_min = -0.6;
    BDT_5_7_max = 0.5;

    BDT_7_10_min = -0.4;
    BDT_7_10_max = 0.6;

    BDT_10_15_min = -0.45;
    BDT_10_15_max = 0.7;

    BDT_15_20_min = -0.6;
    BDT_15_20_max = 0.8;

    BDT_20_50_min = -0.7;
    BDT_20_50_max = 0.9;

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
    RooRealVar BDT_pt_0_2("BDT_pt_0_2", "BDT_pt_0_2", BDT_0_2_min, BDT_0_2_max);
    RooRealVar BDT_pt_2_3("BDT_pt_2_3", "BDT_pt_2_3", BDT_2_3_min, BDT_2_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_3_5", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);

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
    w.import(BDT_pt_0_2);
    w.import(BDT_pt_2_3);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
  }
}
