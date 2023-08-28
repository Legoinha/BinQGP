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

#include <cmath>
#include <math.h>
#include "RooRealProxy.h"
#include "RooAbsReal.h"
// #include "RooDoubleCBFast.h"
// #include "RooDoubleCBFast.cc"
#include "CMS_lumi.C"
#include <algorithm>
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
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBifurGauss.h>
#include "TMath.h"
#include <RooGenericPdf.h>
#include <RooExtendPdf.h>
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
#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <string>
#include <stdio.h>
#include <unordered_map>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "RooCBShape.h"
#include "RooCmdArg.h"

//particle
// 0 = Bu
// 1 = Bs
// 2 = B0
#define particle 0

constexpr bool early = false;
// constexpr bool early = true;

// Use NP Jpsi shape for mass fitting (Err fn + Gaussians)
constexpr bool include_np = true;

const std::vector<TString> BDTdir = {"Bu", "Bs"};
const std::unordered_map<std::string, unsigned> iBDT_bs = {
  {"BDT_pt_7_10", 0},
  {"BDT_pt_10_15", 1},
  {"BDT_pt_15_20", 2},
  {"BDT_pt_20_50", 3}};
const std::unordered_map<std::string, unsigned> iBDT_bp = {
  {"BDT_pt_5_7", 0},
  {"BDT_pt_7_10", 1},
  {"BDT_pt_10_15", 2},
  {"BDT_pt_15_20", 3},
  {"BDT_pt_20_50", 4}};

// only include events with BDT output higher than this lower bound
const std::vector<double> bdt_lower_bound_bp =
  {-0.01,
   -0.12,
   -1.0,
   -1.0,
   -1.0};
const std::vector<double> bdt_lower_bound_bs =
  {-0.01,
   -1.0,
   -1.0,
   -1.0};


// min and max for BDT histograms
const std::vector<double> BDTmin_bs = {-0.1, -0.2, -0.05, -0.1};
const std::vector<double> BDTmax_bs = {0.85, 0.821, 0.92, 0.9};

const std::vector<double> BDTmin_bp = {0.00, -0.00, -0.18, -0.12, -0.1};
const std::vector<double> BDTmax_bp = {0.8, 0.82, 0.74, 0.74, 0.77};
// min and max for BDT histograms

// number of bins
const std::vector<int> BDTnbins_bs = {11, 30, 20, 30};
const std::vector<int> BDTnbins_bp(5, 40);

// initial values for error function / signal ratio
const std::vector<double> ini_f_erf = {0.2, 0.2, 0.2, 0.2, 0.2};

using d_matrix = std::vector<std::vector<double> >;
   // initial values for NP Jpsi fit
const d_matrix ini_jpsi_poly_f = {{0.9, 0.9},
                                  {0.8, 0.7},
                                  {0.8, 0.8},
                                  {0.6, 0.8},
                                  {0.6, 0.6}};
const d_matrix ini_jpsi_p3 = {{0, -7},
                              // {0, -8.4},
                              {0, -1.2},
                              {0, 0.002},
                              {0, -0.15},
                              {0, -0.1}};

const d_matrix ini_jpsi_p2 = {{0, 70},
                              // {0, 42},
                              {0, 14},
                              {0, 0.007},
                              {0, -0.15},
                              {0, -0.1}};
// slope of the continuum bg
const d_matrix ini_jpsi_p1 = {{0, -156},
                              // {0, -248},
                              {0, -37},
                              {0, -0.26},
                              {0, -0.15},
                              {0, -0.12}};

const d_matrix ini_jpisnp_jpsipi_f = {{0, 0.01},
                                      {0, 0.05},
                                      {0, 0.05},
                                      {0, 0.05},
                                      {0, 0.05}};

const d_matrix ini_erf_scale = {{0, 0.05},
                                {0, 0.05},
                                {0, 0.05},
                                {0, 0.05},
                                {0, 0.05}};

const d_matrix ini_p1 = {{0, -23},
                         {0, -23},
                         {0, 0.002},
                         {0, 0},
                         {0, 0}};
const d_matrix ini_p2 = {{0, 9.5},
                         {0, 9.5},
                         {0, 0.003},
                         {0, 0},
                         {0, 0}};
const d_matrix ini_p3 = {{0, -0.9},
                         {0, -0.9},
                         {0, 0},
                         {0, 0},
                         {0, 0}};

using namespace RooStats;
using namespace RooFit;
using namespace std;

std::vector<TH1D*> sideband_subtraction(RooWorkspace& w, std::vector<TString> label, int n_var);
std::vector<TH1D*> splot_method(RooWorkspace& w, std::vector<TString> label, int n_var);

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram_mc(RooWorkspace& w, RooRealVar var, int n, TString weight); 
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n); 
void AddWeights(TTree* t);
void read_MC_DATA_samples(RooWorkspace& w, std::vector<TString> label, TString f_input, TString dataORmc);
void read_samples(RooWorkspace& w, std::vector<TString>, TString fName, TString treeName, TString sample);
void build_pdf (RooWorkspace& w, TString choice, RooArgSet &c_vars, int ipt, int iy);
void fit_jpsinp (RooWorkspace& w, std::string choice, const RooArgSet &c_vars, int ipt, int iy, bool inclusive=false, bool includeSignal=true);
void plot_complete_fit(RooWorkspace& w, RooArgSet &c_vars, TString subname, int iy);
void do_splot(RooWorkspace& w, RooArgSet &c_vars);
TH1D* make_splot(RooWorkspace& w, int n, TString label);
void validate_fit(RooWorkspace* w, RooArgSet &c_vars);
void get_ratio( std::vector<TH1D*>,  std::vector<TH1D*>,  std::vector<TString>, TString, TString);
void constrainVar(TString input_file, TString inVarName, RooArgSet &c_vars, RooArgSet &c_pdfs);
double MC_fit_result(TString input_file, TString inVarName);
TString ystring(int iy);
void save_validation_plot(TCanvas& can, TString name, TString comp, TString ptdir, int iy);
void fix_signal_shape(RooWorkspace& w, bool release=false);
void fix_parameters(RooWorkspace& w, TString pdfName, bool release=false);
void fix_parameters(RooWorkspace& w, RooArgList& parlist, bool release=false);

template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds,
                TString plotName, TString title, Targs... options);

void plot_jpsifit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds,
                  TString plotName, TString title, double nGen, RooRealVar& n_signal);

// change according to what we want to compute: "Bpt"  or  "By"  or  "nMul"

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

# define component 1


#if particle == 0
auto BDTmin = BDTmin_bp;
auto BDTmax = BDTmax_bp;
auto BDTnbins = BDTnbins_bp;
auto indexBDT = iBDT_bp;
auto BDT_lower_bound = bdt_lower_bound_bp;
#elif particle == 1
auto BDTmin = BDTmin_bs;
auto BDTmax = BDTmax_bs;
auto BDTnbins = BDTnbins_bs;
auto indexBDT = iBDT_bs;
auto BDT_lower_bound = bdt_lower_bound_bs;
#endif


const std::vector<int> ptlist = {5, 7, 10, 15, 20, 50};
const std::vector<double> ylist = {0, 1.5, 2.4};
const std::vector<TString> particleList = {"Bu", "Bs"};

   // for B+, ipt is 0--4
   // for Bs, ipt is 1--4

void bmesons_new(int ipt = 3, int iy = 1){

  gSystem->mkdir("./results/Bu" ,true );
  gSystem->mkdir("./results/Bs" ,true );

  gROOT->SetBatch();
  bool inclusive = false;
  if (ipt < 0) {
    inclusive = true;
  }

  int n_var;
  TString input_file_data;
  TString input_file_mc;
  TString input_file_mc_swap;
  // Inclusive pT comparison
  if (ipt < 0) {

    input_file_data = (particle == 0)?
      "/eos/user/h/hmarques/work/data/BPData_nom.root" : "/eos/user/h/hmarques/work/data/BsData_nom.root";
    input_file_mc = (particle == 0)?
      "/eos/user/h/hmarques/work/data/BPMC_nom.root" : "/eos/user/h/hmarques/work/data/BsMC_nom.root";
  } else {
    // Binned pT comparison
    if (particle == 0) {
      input_file_data = Form("/afs/cern.ch/user/t/tsheng/public/forHenrique/trk5/BPData_noBDT_trk5_%i_%i.root",
                             ptlist.at(ipt), ptlist.at(ipt+1));
    } else if(particle == 1) {
      input_file_data = Form("/afs/cern.ch/user/t/tsheng/public/forHenrique/trk5/BsData_noBDT_trk5_%i_%i.root",
                             ptlist.at(ipt), ptlist.at(ipt+1));
    }
    if (particle == 0) {
      input_file_mc = Form("/afs/cern.ch/user/t/tsheng/public/forHenrique/trk5/BPMC_noBDT_trk5_%i_%i.root",
                           ptlist.at(ipt), ptlist.at(ipt+1));
    } else if(particle == 1) {
      input_file_mc = Form("/afs/cern.ch/user/t/tsheng/public/forHenrique/trk5/BsMC_noBDT_trk5_%i_%i.root",
                           ptlist.at(ipt), ptlist.at(ipt+1));
    }
  }


  std::vector<TH1D*> histos_sideband_sub;
  std::vector<TH1D*> histos_mc;
  std::vector<TH1D*> histos_splot;

#if particle == 0
  std::vector<TString> variables = {"BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50", "By", "nMult", "Bpt" /*, "Btrk1Pt", "Btrk1Eta", "BsvpvDisErr", "Btrk1PtErr", "Bchi2cl" , "BsvpvDistance", "Bmumumass", "Bmu1eta","Bmu2eta", "Bmu1pt", "Bmu2pt","Bmu1dxyPV","Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bdtheta", "Balpha"*/};

#elif particle == 1

  std::vector<TString> variables = {"BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50", "By", "Bpt", "nMult", /*"Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "Btrk2Pt", "Btrk2Eta", "Btrk2PtErr" */};
#elif particle == 2

  std::vector<TString> variables = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "BsvpvDistance_2D", "BsvpvDisErr_2D", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "Btrk2Pt", "Btrk2Eta", "Btrk2PtErr", "BDT_pt_0_2", "BDT_pt_2_3", "BDT_pt_3_5", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50","nMult"};
#endif

  if( (particle != 2) && (MC == 1) ){return;} // only fits the MC for B0

  n_var = variables.size();
  cout << "number of variables in TString: "<< n_var << endl;

   RooWorkspace* ws = new RooWorkspace("ws");

   set_up_workspace_variables(*ws);

  /*if( (MC == 1) && (component == 0) ){read_MC_DATA_samples(*ws, variables, input_file_mc, "mc");}
  else if( (MC == 1) && (component == 1)){read_MC_DATA_samples(*ws, variables, input_file_mc_swap, "mc");}*/
  
  if(MC == 0) {
    read_MC_DATA_samples(*ws, variables, input_file_data, "data");
    read_MC_DATA_samples(*ws, variables, input_file_mc, "mc");
  }

  RooArgSet c_vars;

  // use the setting of pT 10-15 for inclusive pT comparison
  if (inclusive) {
    ipt = 2;
    cout << "Using initial values for pT 10-15" << "\n";
  }
  TString signal_shape = (particle == 0 && ipt == 0)? "sig3gauss" : "nominal";
  // TString signal_shape = "nominal";
  cout << "choice:" << signal_shape << "\n";

  build_pdf(*ws, signal_shape, c_vars, ipt, iy);

  if (include_np && (particle == 0)) {
    fit_jpsinp(*ws, "nominal", c_vars, ipt, iy, inclusive);
  }

  TString subname = TString::Format("%i_%i", ptlist.at(ipt), ptlist.at(ipt + 1));
  if (inclusive) {
    subname = "inclusive";
  }
  plot_complete_fit(*ws, c_vars, subname, iy);
  if (early) {return;}

//validate_fit(ws, c_vars);
//  return;


  //SIDEBAND SUBTRACTION (needs to be run after plot_complete_fit)
  histos_sideband_sub = sideband_subtraction(*ws, variables , n_var);
 
  //SPLOT (fixes parameters of the fit -> they need to be unfixed for pT analysis) 
  do_splot(*ws,c_vars); 
  histos_splot = splot_method(*ws, variables, n_var); 

  //MONTE CARLO HISTOGRAMS
  TFile *fin_mc = new TFile(input_file_mc); //use this file to add the weights (to clone original tree) and make data-MC comparisons without weights
  //TFile *fin_mc = new TFile(input_file_reweighted_mc); //use this file to make data-MC comparisons with weights

  TTree* t1_mc;

  if(particle == 0){t1_mc = (TTree*)fin_mc->Get("ntKp");}
  else if(particle == 1){t1_mc = (TTree*)fin_mc->Get("ntphi");}
  else if(particle == 2){t1_mc = (TTree*)fin_mc->Get("ntKstar");}

  std::vector<TString> names;

  for(int i=0; i<n_var; i++){
    TString weight = "weight";                
    histos_mc.push_back(create_histogram_mc(*ws, (*ws->var(variables[i])) , 40, weight));
    names.push_back(TString(variables[i]));}
  
  TString ptdir = "./results/" + particleList.at(particle);
  if (inclusive) {
    ptdir += "/inclusive";
  } else {
    ptdir += Form("/%i_%i", ptlist[ipt], ptlist[ipt + 1]);
  }
gSystem->mkdir(Form("%s/mc_validation_plots", ptdir.Data()) );

  // RATIO BETWEEN DATA (SPLOT) AND MC
  if (weights == 1){ get_ratio(histos_splot, histos_mc, names,"weights.root", ptdir);}
  
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

	cout << ss_comp_mc[i]->GetXaxis()->GetXmax()<< "  max  "<<mc_comp_ss[i]->GetXaxis()->GetXmax()<< endl;
	cout << ss_comp_mc[i]->GetXaxis()->GetXmin()<< "  min  "<<mc_comp_ss[i]->GetXaxis()->GetXmin()<< endl;

    auto rp = new TRatioPlot(ss_comp_mc[i] ,mc_comp_ss[i], "divsym");

    c.SetTicks(0, 1);
    rp->SetH1DrawOpt("E");
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("Data(ss)/MC");
    rp->GetUpperRefYaxis()->SetTitle("normalized entries");
    c.Update();

    TLegend* leg;
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(ss_comp_mc[i]->GetName(), "S. Subtraction", "LE");
    leg->AddEntry(mc_comp_ss[i]->GetName(), "Monte Carlo", "LE");
    leg->SetTextSize(0.03);
    leg->Draw("same");
   
    save_validation_plot(c, names[i], "ss_mc", ptdir, iy);
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
    leg->AddEntry(ss_comp_sp[i]->GetName(), "S. Subtraction", "LE");
    leg->AddEntry(sp_comp_ss[i]->GetName(), "SPlot", "LE");
    leg->SetTextSize(0.03);
    leg->Draw("same");
 
    save_validation_plot(a, names[i], "ss_sp", ptdir, iy);
   }

  //SPlot vs. Monte Carlo
  vector<TH1D*> sp_comp_mc(histos_splot);
  vector<TH1D*> mc_comp_sp(histos_mc);

  for(int i=0; i<n_var; i++){
    TCanvas b;
    mc_comp_sp[i]->SetXTitle(variables[i]);
    mc_comp_sp[i]->SetYTitle("Normalized entries");
    mc_comp_sp[i]->SetStats(0);
    sp_comp_mc[i]->SetStats(0);

    if(particle == 0){mc_comp_ss[i]->SetTitle("");}
    else if (particle == 1){mc_comp_ss[i]->SetTitle("");}
    else if(particle == 2){mc_comp_ss[i]->SetTitle("");}
 
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
    rp->GetLowerRefYaxis()->SetTitle("Data/MC");
    rp->GetUpperRefYaxis()->SetTitle("Normalized entries");
    rp->GetLowerRefGraph()->SetMinimum(-1);
    rp->GetLowerRefGraph()->SetMaximum(3);
    // rp->GetUpperRefXaxis()->SetRange(0., 0.825);
    // rp->GetLowerRefXaxis()->SetRange(0., 0.825);
    //b.Update();

	  TLatex* texB = new TLatex(0.5,0.5,"");
	  if(particle==1){ texB = new TLatex(0.15,0.85, "B^{0}_{s}");}
	  if(particle==0){ texB = new TLatex(0.15,0.85, "B^{+}");}
		texB->SetNDC();
		texB->SetTextFont(62);
		texB->SetTextSize(0.04);
		texB->SetLineWidth(2);
		texB->Draw();    
 
    TLegend* leg;	
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(mc_comp_sp[i]->GetName(), "Monte Carlo", "LE");
    leg->AddEntry(sp_comp_mc[i]->GetName(), "sPlot", "LE");
    leg->SetTextSize(0.03);
    leg->Draw("same");
	
    save_validation_plot(b, names[i], "mc_sp", ptdir, iy);
  }
 
  //Sideband subtraction vs. Monte Carlo vs SPlot
  vector<TH1D*> sp_comp(histos_splot);
  vector<TH1D*> mc_comp(histos_mc);
  vector<TH1D*> ss_comp(histos_sideband_sub);

  for(int i=0; i<n_var; i++){
    TCanvas d;	
    mc_comp[i]->SetXTitle(variables[i]);
    mc_comp[i]->SetYTitle("Normalized entries");
    mc_comp[i]->SetStats(0);
    sp_comp[i]->SetStats(0);
    ss_comp[i]->SetStats(0);

    if(particle == 0){ss_comp[i]->SetTitle("");}
    else if (particle == 1){ss_comp[i]->SetTitle("");}
    else if(particle == 2){ss_comp[i]->SetTitle("");}

    //normalization
    mc_comp[i]->Scale(1/mc_comp[i]->Integral());
    sp_comp[i]->Scale(1/sp_comp[i]->Integral());
    ss_comp[i]->Scale(1/ss_comp[i]->Integral());	
	
    // y axis: maximum and minimum
    double ymax = std::max(mc_comp[i]->GetMaximum(),
                           std::max(sp_comp[i]->GetMaximum(),
                                    ss_comp[i]->GetMaximum()));
    double ymin = std::min(mc_comp[i]->GetMinimum(),
                           std::min(sp_comp[i]->GetMinimum(),
                                    ss_comp[i]->GetMinimum()));
    mc_comp[i]->GetYaxis()->SetRangeUser(0.1 * ymin, 1.1 * ymax);
    ss_comp[i]->GetYaxis()->SetRangeUser(0.1 * ymin, 1.1 * ymax);

    auto rp_sp = new TRatioPlot(sp_comp[i], mc_comp[i], "divsym");
    rp_sp->Draw();
    TGraph* ratio_sp = rp_sp->GetLowerRefGraph();
    d.Clear();

    auto rp_ss = new TRatioPlot(ss_comp[i], mc_comp[i], "divsym");
    d.SetTicks(0, 1);
    rp_ss->SetH1DrawOpt("E");
    rp_ss->Draw();
    rp_ss->GetLowerRefYaxis()->SetTitle("Data/MC");
    rp_ss->GetUpperRefYaxis()->SetTitle("Normalized entries");
    rp_ss->GetLowerPad()->cd();
    ratio_sp->Draw("e same");
    TGraph* ratio_ss = rp_ss->GetLowerRefGraph();
    // Limit the range of y axis for ratio plots
    double ratio_range = 2;
    auto ratio_ss_max = TMath::MaxElement(ratio_ss->GetN(), ratio_ss->GetY());
    auto ratio_sp_max = TMath::MaxElement(ratio_sp->GetN(), ratio_sp->GetY());
    auto ratio_ss_min = TMath::MinElement(ratio_ss->GetN(), ratio_ss->GetY());
    auto ratio_sp_min = TMath::MinElement(ratio_sp->GetN(), ratio_sp->GetY());
    auto minmax = std::vector<double> {ratio_ss_min, ratio_ss_max,
                                       ratio_sp_min, ratio_sp_max};
    std::transform(minmax.cbegin(), minmax.cend(), minmax.begin(), [](double m) {
      return std::abs(m - 1);});
    // make the ratio range symmetric
    if (*std::max_element(minmax.begin(), minmax.end()) < ratio_range) {
      ratio_range = *std::max_element(minmax.begin(), minmax.end());
    }
    ratio_ss->SetMinimum(1 - ratio_range);
    ratio_ss->SetMaximum(1 + ratio_range);

    rp_ss->GetUpperPad()->cd();
    sp_comp[i]->Draw("e same");
    d.Update();

    TLegend* leg;
    leg = new TLegend(0.7, 0.9, 0.9, 1.0);
    leg->AddEntry(sp_comp[i]->GetName(), "SPlot", "LE");
    leg->AddEntry(ss_comp[i]->GetName(), "S. Subtraction", "LE");
    leg->AddEntry(mc_comp[i]->GetName(), "Monte Carlo", "LE");
    leg->SetTextSize(0.03);
    leg->Draw("same");

    save_validation_plot(d, names[i], "ss_mc_sp", ptdir, iy);
  }

//comparisons end
}
//main function ends

double MC_fit_result(TString input_file, TString inVarName){

  TFile* f = new TFile(input_file);
  RooFitResult* fitresult;
  fitresult = (RooFitResult*)f->Get("fitresult_model_data");
  RooRealVar* var = (RooRealVar*)fitresult->floatParsFinal().find(inVarName);
 
  return var->getVal();}

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

//get the ratio between the data (splot method) and the MC and save it in a root file
void get_ratio( std::vector<TH1D*> data, std::vector<TH1D*> mc,
                std::vector<TString> v_name, TString filename,
                TString ptdir){

  TString dir_name;
  dir_name = ptdir + "/mc_validation_plots/weights/";
  gSystem->mkdir( dir_name.Data() );
  if(particle == 2){dir_name = "./results/B0/mc_validation_plots/weights/";}

  gSystem->Exec("mkdir -p " + dir_name);
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
    gSystem->Exec("mkdir -p " + dir_name );
    c.SaveAs(dir_name+v_name.at(i) + "_weights.pdf");   
    //output: a root file and plots 
  }
  f_wei->Close();
  return;
}
//get_ratio ends

//read_data
void read_MC_DATA_samples(RooWorkspace& w, std::vector<TString> label, TString f_input, TString dataORmc){
  TString treeName;
  if(particle == 0) {
    treeName = "ntKp";
    cout << "reading J/psi inclusive file" << "\n";
    std::vector<TString> jpsi_vars = {"Bgen", "By", "Bpt", "BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15","BDT_pt_15_20", "BDT_pt_20_50"};
    read_samples(w, jpsi_vars, "/eos/user/h/hmarques/work/data/jpsinp_nom.root", treeName, "jpsinp");
  } 
  else if (particle == 1) { treeName = "ntphi";}
  else {treeName = "ntKstar";}
  read_samples(w, label, f_input, treeName, dataORmc.Data());
}

// work horse to read data/MC/jpsi
void read_samples(RooWorkspace& w, std::vector<TString> label, TString fName, TString treeName, TString sample){
  TFile* fin = new TFile(fName);
  TTree* t1;

  t1 = (TTree*) fin->Get(treeName);
  int n_var = label.size();
  cout << "size: " << n_var  << "\n";

  RooArgList arg_list ("arg_list");
  arg_list.add(*(w.var("Bmass")));  // read the fitting variable

  // read additional variables
  for(auto lab : label){
    cout << lab << "\n";
    arg_list.add(*(w.var(lab)));
    cout << "added " << lab << "\n";
  }

  RooDataSet* ds = new RooDataSet(sample, sample, t1, arg_list);
  cout << "input filename = " << fName << "; entries before FID region: " << ds->sumEntries() << endl;
  ds = (RooDataSet*)ds->reduce("(Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10)");  //FID REGION
  cout << "input filename = " << fName << "; entries after FID region: " << ds->sumEntries() << endl;

  //should be optimized: use the global arrays and vectors to make it easier changing the values
  #if particle == 0
  ds = (RooDataSet*)ds->reduce("(BDT_pt_5_7 > 0.08 && Bpt >= 5 && Bpt < 7) || (BDT_pt_7_10 > 0.07 && Bpt >= 7 && Bpt < 10) || (BDT_pt_10_15 > 0 && Bpt >= 10 && Bpt < 15) || (BDT_pt_15_20 > 0.02 && Bpt >= 15 && Bpt < 20) || (BDT_pt_20_50 > 0.04 && Bpt >= 20 && Bpt < 50) || (Bpt >= 50 && Bpt < 60) ");
  #elif particle == 1
  ds = (RooDataSet*)ds->reduce("(BDT_pt_7_10 > 0.06 && Bpt >= 7 && Bpt < 10) || (BDT_pt_10_15 > -0.04 && Bpt >= 10 && Bpt < 15) || (BDT_pt_15_20 > 0.05 && Bpt >= 15 && Bpt < 20) || (BDT_pt_20_50 > -1 && Bpt >= 20 && Bpt < 50) ");
  #endif
  cout << "input filename = " << fName << "; entries before FID region: " << ds->sumEntries() << endl;

  w.import(*ds, Rename(sample));
}

//build_pdf
void build_pdf(RooWorkspace& w, TString choice, RooArgSet &c_vars, int ipt=3, int iy=1){

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*)w.data("data");

  TString input_file_RT ;
  TString input_file_WT ;
  TFile* file_RT ;
  TFile* file_WT ;
  if( (particle == 2) && (MC==0)){
  input_file_RT = "./results/B0/MC/RT_fit.root" ;
  input_file_WT = "./results/B0/MC/WT_fit.root" ;
  file_RT = new TFile(input_file_RT);
  file_WT = new TFile(input_file_WT);}
  RooArgSet c_pdfs_RT;
  RooArgSet c_pdfs_WT; 

  //  RooDataSet* reduceddata_central;
  double left;
  double right;
  double mass_peak;

  if(particle == 0){
    left = 5.2;
    right = 5.4;
    mass_peak = 5.279;
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

  // Variable initialization
  RooRealVar* f_swap = 0;
  RooRealVar* mean = 0;
  RooRealVar* mean_swp = 0;
  RooRealVar* sigma1 = 0;
  // RooRealVar* sigma2 = 0;
  RooProduct* sigma2;
  RooRealVar* ratio_sigma12 = 0;
  RooRealVar* ratio_sigma13 = 0;
  // RooRealVar* sigma3 = 0;
  RooProduct* sigma3 = 0;
  RooRealVar* sigma_cb1 = 0;
  RooRealVar* sigma_cb2 = 0;
  RooRealVar* sigma1_swp = 0;
  RooRealVar* sigma2_swp = 0;
  RooRealVar* alpha1 = 0;
  RooRealVar* alpha2 = 0;
  RooRealVar* alpha1_swp = 0;
  RooRealVar* alpha2_swp = 0;
  RooRealVar* n1 = 0;
  RooRealVar* n2 = 0;
  RooRealVar* n3 = 0;
  RooRealVar* n1_swp = 0;
  RooRealVar* n2_swp = 0;
  RooRealVar* cofs = 0;
  RooRealVar* cofs1 = 0;
  RooRealVar* m_nonprompt_scale = 0;
  RooRealVar* m_nonprompt_shift = 0;
  RooRealVar* lambda = 0;
  RooRealVar* slope = 0;
  RooProduct* m_jpsipi_width = 0;
  RooRealVar* m_jpsipi_mean1 = 0;
  RooRealVar* m_jpsipi_sigma1l = 0;
  RooRealVar* m_jpsipi_sigma1r = 0;
  RooRealVar* m_jpsipi_mean2 = 0;
  RooRealVar* m_jpsipi_mean3 = 0;
  RooRealVar* m_jpsipi_sigma2 = 0;
  RooRealVar* m_jpsipi_sigma3 = 0;
  RooRealVar* m_jpsipi_fraction2 = 0;
  RooRealVar* m_jpsipi_fraction3 = 0;
  RooProduct* sigma1_fix = 0;
  RooProduct* sigma2_fix = 0;
  RooProduct* sigma1_swp_fix = 0;
  RooRealVar* scale_factor = new RooRealVar("scale_factor", "scale_factor", 1., 0., 2.);
  RooRealVar* mean_difference = 0;

  RooRealVar* p1 = 0;
  RooRealVar* p2 = 0;
  RooRealVar* p3 = 0;

  if( (particle == 2) && (MC == 0)){

if( choice != "scale_factor"){ 
    cout << "Applying constraints" << endl;
    constrainVar(input_file_RT, "sigma_cb1", c_vars, c_pdfs_RT);
    constrainVar(input_file_RT, "sigma_cb2", c_vars, c_pdfs_RT);
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

 cout << "Initialisating variables with constraints" << endl;
    // mis-tag fraction
    double n_rt = MC_fit_result(input_file_RT, "n_signal");
    double n_wt = MC_fit_result(input_file_WT, "n_signal_swp");
    double fraction = n_wt / (n_wt + n_rt);
    cout << "fraction = " << fraction << endl;
    f_swap = new RooRealVar("f_swap","f_swap", fraction, 0., 1.);       
    f_swap->setConstant();
    cout << "f_swap = " << f_swap->getVal() << endl; 

// RT Component (CB + CB)
    mean = new RooRealVar("mean", "mean", MC_fit_result(input_file_RT, "mean"), MC_fit_result(input_file_RT, "mean")-0.1, MC_fit_result(input_file_RT, "mean")+0.1);
    sigma_cb1 = new RooRealVar("sigma_cb1", "sigma_cb1", MC_fit_result(input_file_RT, "sigma_cb1"), 0.005, 0.5);
    sigma_cb2 = new RooRealVar("sigma_cb2", "sigma_cb2", MC_fit_result(input_file_RT, "sigma_cb2"), 0.005 ,0.5);
    alpha1 = new RooRealVar("alpha1", "alpha1", MC_fit_result(input_file_RT, "alpha1"), 0., 20.);
    alpha2 = new RooRealVar("alpha2", "alpha2", MC_fit_result(input_file_RT, "alpha2"), 0., 20.);
    n1 = new RooRealVar("n1", "n1", MC_fit_result(input_file_RT, "n1"), 0., 300.);
    n2 = new RooRealVar("n2", "n2", MC_fit_result(input_file_RT, "n2"), 0., 300.);
    cofs = new RooRealVar("cofs", "cofs", MC_fit_result(input_file_RT, "cofs"), 0., 1.);

// WT Component (double CB)
    mean_swp = new RooRealVar("mean_swp", "mean_swp", MC_fit_result(input_file_WT, "mean_swp"), MC_fit_result(input_file_WT, "mean_swp")-0.1, MC_fit_result(input_file_WT, "mean_swp")+0.1);
    sigma1_swp = new RooRealVar("sigma1_swp", "sigma1_swp", MC_fit_result(input_file_WT, "sigma1_swp"), 0.005, 0.5);
    alpha1_swp = new RooRealVar("alpha1_swp", "alpha1_swp", MC_fit_result(input_file_WT, "alpha1_swp"), 0., 20.);
    alpha2_swp = new RooRealVar("alpha2_swp", "alpha2_swp", MC_fit_result(input_file_WT, "alpha2_swp"), 0., 20.);
    n1_swp = new RooRealVar("n1_swp", "n1_swp", MC_fit_result(input_file_WT, "n1_swp"), 0., 300.);
    n2_swp = new RooRealVar("n2_swp", "n2_swp", MC_fit_result(input_file_WT, "n2_swp"), 0., 300.);

    /*elif(choice == "scale_factor"){
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
      }*/
  }else{
    cout << "Initialising variables without constraints" << endl;
    mean = new RooRealVar("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.2);
    sigma_cb1 = new RooRealVar("sigma_cb1", "sigma_cb1", 0.012,0.010,0.030);
    sigma_cb2 = new RooRealVar("sigma_cb2", "sigma_cb2", 0.012,0.010,0.030);
    alpha1 = new RooRealVar("alpha1", "alpha1", 5., 0., 20.); 
    alpha2 = new RooRealVar("alpha2", "alpha2", 5., 0., 20.); 
    n1 = new RooRealVar("n1", "n1", 10., 0., 300.);
    n2 = new RooRealVar("n2", "n2", 10., 0., 300.);
    cofs = new RooRealVar("cofs", "cofs", 0.3, 0., 1.);
    mean_swp = new RooRealVar("mean_swp","mean_swp",mass_peak,mass_peak-0.1,mass_peak+0.2);   
    sigma1_swp = new RooRealVar("sigma1_swp","sigma1_swp",0.02,0.005,0.5);
    alpha1_swp = new RooRealVar("alpha1_swp", "alpha1_swp", 5., 0., 20.);
    alpha2_swp = new RooRealVar("alpha2_swp", "alpha2_swp", 5., 0., 20.);
    n1_swp = new RooRealVar("n1_swp", "n1_swp", 10., 0., 300.);
    n2_swp = new RooRealVar("n2_swp", "n2_swp", 10., 0., 300.);
    }

   sigma1 = new RooRealVar("sigma1","sigma1",0.02,0.005,0.025);
   ratio_sigma12 = new RooRealVar("ratio_sigma12","ratio_sigma12", 2, 0.01, 10);
   ratio_sigma13 = new RooRealVar("ratio_sigma13","ratio_sigma13", 2, 0.01, 10);
   // sigma2 = new RooRealVar("sigma2","sigma2",0.01,0.005,0.5);
   sigma2 = new RooProduct("sigma2", "sigma2", RooArgList(*sigma1, *ratio_sigma12));
   f_swap = new RooRealVar("f_swap","f_swap",0.,0.,1.);
   cofs1 = new RooRealVar("cofs1", "cofs1", 0.3, 0., 1.);
   lambda = new RooRealVar("lambda","lambda",-2.,-10.,1.0);
   slope = new RooRealVar("slope","slope", -500, -1000, 1);
   n3 = new RooRealVar("n3", "n3", 100., 0., 400.);
   // sigma3 = new RooRealVar("sigma3","sigma3",0.012,0.010,0.030);
   sigma3 = new RooProduct("sigma3", "sigma3", RooArgList(*sigma1, *ratio_sigma13));
   p1 = new RooRealVar("p1", "p1", ini_p1[ipt][iy], -150., 150.);
   p2 = new RooRealVar("p2", "p2", ini_p2[ipt][iy], -15., 15.);
   p3 = new RooRealVar("p3", "p3", ini_p3[ipt][iy], -2., 2.);

  m_nonprompt_scale = new RooRealVar("m_nonprompt_scale", "m_nonprompt_scale",
                                     ini_erf_scale[ipt][iy], 0, 0.1);
  m_nonprompt_shift = new RooRealVar("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  RooRealVar jpsipi_to_signal_width_ratio("jpsipi_to_signal_width_ratio",
                                          "jpsipi_to_signal_width_ratio",
                                          1 / sigma1->getVal());
  m_jpsipi_width = new RooProduct("m_jpsipi_width","m_jpsipi_width",
                                  RooArgList(jpsipi_to_signal_width_ratio, *sigma1));
  m_jpsipi_mean1 = new RooRealVar("m_jpsipi_mean1","m_jpsipi_mean1",5.34693, 5.3, 5.5);
  m_jpsipi_sigma1l = new RooRealVar("m_jpsipi_sigma1l","m_jpsipi_sigma1l",0.0290762,0.010,0.150);
  m_jpsipi_sigma1r = new RooRealVar("m_jpsipi_sigma1r","m_jpsipi_sigma1r",0.0652519,0.010,0.200);
  m_jpsipi_mean2 = new RooRealVar("m_jpsipi_mean2","m_jpsipi_mean2",5.46876, 5.3, 5.6);
  m_jpsipi_mean3 = new RooRealVar("m_jpsipi_mean3","m_jpsipi_mean3",5.48073);
  m_jpsipi_sigma2 = new RooRealVar("m_jpsipi_sigma2","m_jpsipi_sigma2",0.0994712,0.020,0.500);
  m_jpsipi_sigma3 = new RooRealVar("m_jpsipi_sigma3","m_jpsipi_sigma3",0.330152,0.020,0.500);
  m_jpsipi_fraction2 = new RooRealVar("m_jpsipi_fraction2","m_jpsipi_fraction2",0.234646,0.0,1.0);
  m_jpsipi_fraction3 = new RooRealVar("m_jpsipi_fraction3","m_jpsipi_fraction3",0.114338,0.0,1.0);

  // scaled widths
  RooRealVar m_jpsipi_sigma2l("m_jpsipi_sigma2l","m_jpsipi_sigma2l",0.0994712,0.020,0.500);
  RooRealVar m_jpsipi_sigma2r("m_jpsipi_sigma2r","m_jpsipi_sigma2r",0.0994712,0.020,0.500);
  RooProduct jpsipi_sigma1l("jpsipi_sigma1l", "jpsipi_sigma1l",
                            RooArgList(*m_jpsipi_width, *m_jpsipi_sigma1l));
  RooProduct jpsipi_sigma1r("jpsipi_sigma1r", "jpsipi_sigma1r",
                            RooArgList(*m_jpsipi_width, *m_jpsipi_sigma1r));
  RooProduct jpsipi_sigma2("jpsipi_sigma2", "jpsipi_sigma2",
                            RooArgList(*m_jpsipi_width, *m_jpsipi_sigma2));
  RooProduct jpsipi_sigma2l("jpsipi_sigma2l", "jpsipi_sigma2l",
                           RooArgList(*m_jpsipi_width, m_jpsipi_sigma2l));
  RooProduct jpsipi_sigma2r("jpsipi_sigma2r", "jpsipi_sigma2r",
                           RooArgList(*m_jpsipi_width, m_jpsipi_sigma2r));

//SIGNAL PDF
cout << "Defining PDF" << endl;
  //GAUSSIANS
  RooGaussian* signal1 = new RooGaussian("signal1","signal_gauss1",Bmass,*mean,*sigma1);
  RooGaussian* signal2 = new RooGaussian("signal2","signal_gauss2",Bmass,*mean,*sigma2); 
	// double gaussian
	RooAddPdf* signal = new RooAddPdf("signal", "signal", RooArgList(*signal1,*signal2),*cofs);
       // single gaussian
       // RooGaussian* signal = new RooGaussian(*signal1, "signal");
      // triple gaussian
  		RooGaussian* signal3 = new RooGaussian("signal3","signal3",Bmass, *mean, *sigma3);
  		RooAddPdf* signal_triple = new RooAddPdf("signal_triple","signal_triple",RooArgList(*signal1,*signal2,*signal3),RooArgList(*cofs,*cofs1));
	
  //CRYSTAL BALL
  RooCBShape* CB1 = new RooCBShape("CB1", "CB1",Bmass, *mean, *sigma_cb1,*alpha1,*n1);
  RooCBShape* CB2 = new RooCBShape("CB2", "CB2",Bmass, *mean, *sigma_cb2,*alpha2,*n2);
  	// double crystal ball
  	RooAddPdf* sum_CB = new RooAddPdf("sum_CB","sum_CB",RooArgList(*CB1,*CB2),*cofs);
 
  //GAUSSIANS AND CRYSTAL BALL COMBINATIONS
  RooAddPdf* gauss_CB = new RooAddPdf("gauss_CB","gauss_CB",RooArgList(*signal1,*CB1),*cofs);
  RooAddPdf* two_gauss_CB = new RooAddPdf("two_gauss_CB","two_gauss_CB",RooArgList(*signal,*CB1),*cofs);

  //FOR B0, WT COMPONENT             
  // RooDoubleCBFast* double_CB_swp = new RooDoubleCBFast("double_CB_swp", "double_CB_swp", Bmass, *mean_swp, *sigma1_swp, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);

//BACKGROUND//

  //error function (for JPsi X peaking background)
  RooGenericPdf erfn ("erfn", "erfn", "TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)",
                      RooArgList(Bmass, *m_nonprompt_scale, *m_nonprompt_shift));

  // non-prompt jpsi
  RooRealVar np_p0("np_p0", "np_p0", 100, 0, 2000);
  RooRealVar np_p1("np_p1", "np_p1", ini_jpsi_p1[ipt][iy], -1000., 1000.);
  RooRealVar np_p2("np_p2", "np_p2", ini_jpsi_p2[ipt][iy], -1000., 1000.);
  RooRealVar np_p3("np_p3", "np_p3", ini_jpsi_p3[ipt][iy], -10, 10);
  RooRealVar np_mean1("np_mean1", "np_mean1", 5.1, 5., 5.3);
  RooRealVar np_sigma1("np_sigma1", "np_sigma1", 0.02, 0.005, 0.05);
  RooRealVar np_mean2("np_mean2", "np_mean2", 5.55, 5, 6);
  RooRealVar np_sigma2("np_sigma2", "np_sigma2", 0.2, 0.05, 0.5);

  // For B->Jpsi pi, bifur Gaussian + Gaussian
  RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1", "m_jpsipi_gaussian1", Bmass,
                                   *m_jpsipi_mean1, jpsipi_sigma1l, jpsipi_sigma1r);
  // RooGaussian* m_jpsipi_gaussian2 =
  //   new RooGaussian("m_jpsipi_gaussian2", "m_jpsipi_gaussian2", Bmass,
  //                   *m_jpsipi_mean1, jpsipi_sigma2);
  RooBifurGauss m_jpsipi_gaussian2("m_jpsipi_gaussian2", "m_jpsipi_gaussian2", Bmass,
                                   *m_jpsipi_mean1, jpsipi_sigma2l, jpsipi_sigma2r);
  RooGaussian* m_jpsipi_gaussian3 = new RooGaussian("m_jpsipi_gaussian3","m_jpsipi_gaussian3",Bmass,*m_jpsipi_mean3,*m_jpsipi_sigma3);
  // RooAddPdf* jpsipi = new RooAddPdf("jpsipi","jpsipi",RooArgList(*m_jpsipi_gaussian3,*m_jpsipi_gaussian2,*m_jpsipi_gaussian1),RooArgList(*m_jpsipi_fraction3,*m_jpsipi_fraction2));
  RooAddPdf* jpsipi = new RooAddPdf("jpsipi", "jpsipi", RooArgList(m_jpsipi_gaussian2, m_jpsipi_gaussian1), RooArgList(*m_jpsipi_fraction2));
  // RooAbsPdf* jpsipi = new RooBifurGauss(*m_jpsipi_gaussian1, "jpsipi");


  RooPolynomial* poly_jpsi = 0;
  if (ipt <= 2) {
    poly_jpsi = new RooPolynomial("poly_jpsi", "poly_jpsi", Bmass, RooArgList(np_p1, np_p2, np_p3));
  } else {
    poly_jpsi = new RooPolynomial("poly_jpsi", "poly_jpsi", Bmass, RooArgList(np_p1));
  }

  RooRealVar jpsinp_poly_fraction("jpsinp_poly_fraction", "fraction",
                                  ini_jpsi_poly_f[ipt][iy], 0.01, 1);
  // fraction of B+ -> J/psi pi+ in the NP
  RooRealVar jpsinp_jpsipi_fraction("jpsinp_jpsipi_fraction", "fraction",
                                    ini_jpisnp_jpsipi_f[ipt][iy], 0.001, 1);
  RooAbsPdf* erf = 0;
  erf = new RooGenericPdf(erfn, "erf");
  RooAddPdf model_jpsinp_cont("m_jpsinp_cont", "model for jpsi nonprompt bg",
                             RooArgList(*poly_jpsi, *erf),
                             RooArgList(jpsinp_poly_fraction));
  // erf = new RooAddPdf(model_jpsinp, "erf");
  w.import(model_jpsinp_cont);
  w.import(*jpsipi);

  // polynomials for continuum background
  RooRealVar p2nd2("p2nd2", "", -0.01, -20., 20.);
  RooRealVar p2ratio("p2ratio", "", -2 * 5.7, -50, 50);
  // RooRealVar p2nd1("p2nd1", "", -2 * 5.7 * p2nd2.getVal(), -1000., 1000.);
  RooProduct p2nd1("p2nd1", "", RooArgList(p2nd2, p2ratio));

  // RooRealVar cont_mean1("cont_mean1", "cont_mean1", 5.7, 5.5, 5.8);
  RooRealVar cont_mean1("cont_mean1", "cont_mean1", 5.7, 5.2, 5.8);
  RooRealVar cont_sigma1("cont_sigma1", "cont_sigma1", 0.4, 0.05, 2);
  RooGaussian cont_g1("cont_g1","cont_g1", Bmass, cont_mean1, cont_sigma1);
  RooRealVar cont_g1_fraction("cont_g1_fraction", "fraction", 0.2, 0., 1);

  RooRealVar cont_mean2("cont_mean2", "cont_mean2", 5.5, 5.3, 5.6);
  RooRealVar cont_sigma2("cont_sigma2", "cont_sigma2", 0.1, 0.05, 1);
  RooGaussian cont_g2("cont_g2","cont_g2", Bmass, cont_mean2, cont_sigma2);
  RooRealVar cont_g2_fraction("cont_g2_fraction", "fraction", 0.2, 0., 1);

  RooRealVar cont_p1_1("cont_p1_1", "cont_p1_1", 0, -10, 10);
  RooPolynomial cont_p1("cont_p1", "cont_p1", Bmass, RooArgList(cont_p1_1), 1);

  RooAbsPdf* fit_side = 0;
  bool use_polynomial_for_background = (ipt <= 1);
  const int poly_order = 3 - ipt;

  switch (poly_order) {
  case 3: {
    fit_side = new RooPolynomial("fit_side", "fit_side", Bmass, RooArgList(*p1, *p2, *p3), 1);
    break;
  } case 2: {
      fit_side = new RooPolynomial("fit_side", "fit_side", Bmass, RooArgList(*p1, *p2, *p3), 1);
      break;
  } case 1: {
      // fit_side = new RooPolynomial("fit_side", "fit_side", Bmass, RooArgList(p2nd1, p2nd2), 1);
      fit_side = new RooPolynomial("fit_side", "fit_side", Bmass, RooArgList(*p1, *p2, *p3), 1);
    break;
  } default:
  //exponential (for combinatorial background)
  fit_side = new RooExponential("fit_side","fit_side",Bmass,*lambda);
    break;
  }
  fit_side->Print();


  // 1st order polynomial (combinatorial background - pdf systematics)
  RooPolynomial* poly_bkg = new RooPolynomial("poly_bkg","poly_bkg",Bmass,*slope);


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
 
cout << "Definig B0 model" << endl;
  // B0 FULL MODEL
  // RT component
  RooCBShape* cb1_rt = new RooCBShape("cb1_rt","cb1_rt",Bmass,*mean,*sigma_cb1,*alpha1,*n1);
  RooCBShape* cb2_rt = new RooCBShape("cb2_rt","cb2_rt",Bmass,*mean,*sigma_cb2,*alpha2,*n2);
  RooAddPdf* sum_cb_rt = new RooAddPdf("sum_cb_rt","sum_cb_rt",RooArgList(*cb1_rt,*cb2_rt),*cofs);
  RooRealVar* RT_yield = new RooRealVar("RT_yield", "RT_yield", RT_yield_initial, 0.,  data->sumEntries());
  RooArgList constr_rt_list = RooArgList(c_pdfs_RT);
  constr_rt_list.add(*sum_cb_rt);
  RooProdPdf* rt_pdf = new RooProdPdf("rt_pdf", "rt_pdf", constr_rt_list);   

  // WT component
  RooFormulaVar* mass_swp;
  // RooDoubleCBFast* double_CB_wt;
  if((particle == 2)&&(MC==0)){ 
    mean_difference = new RooRealVar("mean_difference", "mean_difference", MC_fit_result(input_file_WT, "mean_swp") - MC_fit_result(input_file_RT, "mean"), -2, 2, "GeV");
    RooProduct* mean_difference_fix = 0;
    if(choice == "scale_factor"){
      mean_difference->setConstant();
      mean_difference_fix = new RooProduct("mean_difference_fix", "mean_difference_fix", RooArgList(*scale_factor,*mean_difference));
      mass_swp = new RooFormulaVar("mass_swp", "mass_swp", "@0+@1", RooArgList(*mean,*mean_difference_fix));
    }else if (choice != "scale_factor"){mass_swp = new RooFormulaVar("mass_swp", "mass_swp", "@0+@1", RooArgList(*mean,*mean_difference));}

  // double_CB_wt = new RooDoubleCBFast("double_CB_wt", "double_CB_wt", Bmass, *mass_swp, *sigma1_swp, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);
  }
  RooProduct* WT_yield = new RooProduct("WT_yield","WT_yield",RooArgList(*f_swap,*RT_yield));
  
  // Scale factor
  RooCBShape* cb1_rt_sf;
  RooCBShape* cb2_rt_sf; 
  RooAddPdf* sum_cb_rt_sf; 
  // RooDoubleCBFast* double_CB_wt_sf;
  // if((particle == 2) && (MC == 0) && (choice == "scale_factor")){
  //   cb1_rt_sf = new RooCBShape("cb1_rt_sf","cb1_rt_sf",Bmass,*mean,*sigma1_fix,*alpha1,*n1);
  //   cb2_rt_sf = new RooCBShape("cb2_rt_sf","cb2_rt_sf",Bmass,*mean,*sigma2_fix,*alpha2,*n2);
  //   sum_cb_rt_sf = new RooAddPdf("sum_cb_rt_sf","sum_cb_rt_sf",RooArgList(*cb1_rt_sf,*cb2_rt_sf),*cofs);
  //   double_CB_wt_sf = new RooDoubleCBFast("double_CB_wt_sf", "double_CB_wt_sf", Bmass, *mass_swp, *sigma1_swp_fix, *alpha1_swp, *n1_swp, *alpha2_swp, *n2_swp);}

  // NORMALISATIONS
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar* n_signal = new RooRealVar("n_signal","n_signal",n_signal_initial,0.,(data->sumEntries())*2);
  RooRealVar* n_signal_swp = new RooRealVar("n_signal_swp","n_signal_swp",n_signal_initial,0.,(data->sumEntries())*2);
  RooRealVar* n_combinatorial = new RooRealVar("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
  RooRealVar* f_erf = new RooRealVar("f_erf","f_erf", ini_f_erf[ipt], 0.02, 5);
  RooProduct* n_erf = new RooProduct("n_erf","n_erf",RooArgList(*n_signal,*f_erf));
  // RooRealVar* f_jpsipi = new RooRealVar("f_jpsipi","f_jpsipi",0.03996, 0.038, 0.040);
  // RooRealVar* f_jpsipi = new RooRealVar("f_jpsipi","f_jpsipi",0.03996, 0.01, 1);
  // RooProduct* n_jpsipi = new RooProduct("n_jpsipi","n_jpsipi",RooArgList(*n_signal,*f_jpsipi)); 
  double n_nonprompt_initial = 0.1 * n_combinatorial_initial;

  // don't fit the peaking backgrounds for lower pT
  // if (ipt < 2) {
    // f_erf->setConstant(kTRUE);
  // }

  //mean difference constraint
  RooFitResult* fr;
  RooFitResult* fr_swp;
  RooRealVar* constrained_mean;
  RooRealVar* constrained_mean_swp;
  double mean_diff_val ;
  double mean_diff_error ;
  RooGaussian* mean_constr ;
  RooArgList constr_wt_list;
  RooProdPdf* wt_pdf ;
  if((particle == 2) && (MC == 0)){
  fr = (RooFitResult*)file_RT->Get("fitresult_model_data");
  fr_swp = (RooFitResult*)file_WT->Get("fitresult_model_data");
  constrained_mean =  (RooRealVar*)fr->floatParsFinal().find("mean");
  constrained_mean_swp =  (RooRealVar*)fr_swp->floatParsFinal().find("mean_swp");
  mean_diff_val = (constrained_mean_swp->getVal()) - (constrained_mean->getVal());
  mean_diff_error = sqrt( pow(constrained_mean->getError(),2) + pow(constrained_mean_swp->getError(),2) );
  mean_constr = new RooGaussian("mean_constr", "mean_constr", *mean_difference, RooConst(mean_diff_val), RooConst(mean_diff_error));
  c_vars.add(*mean_difference);
  constr_wt_list = RooArgList(c_pdfs_WT);
  constr_wt_list.add(*mean_constr);
  // constr_wt_list.add(*double_CB_wt);
  wt_pdf = new RooProdPdf("wt_pdf", "wt_pdf", constr_wt_list);
  }


  RooAddPdf model_np("model_nonprompt", "erf and jpsipi",
                     RooArgList(*jpsipi, *erf),
                     RooArgList(jpsinp_jpsipi_fraction));
  RooRealVar n_nonprompt("n_nonprompt", "n_nonprompt",
                         n_nonprompt_initial, 0., data->sumEntries());
  RooRealVar jpsipi_to_signal_ratio("jpsipi_to_signal_ratio", "jpsipi_to_signal_ratio",
                                    0.05, 0, 1);
  RooProduct* n_jpsipi = new RooProduct("n_jpsipi_by_signal",
                                        "number of jpsi pi with fixed ratio to n_signal",
                                        RooArgList(*n_signal, jpsipi_to_signal_ratio));

  if(particle == 0){//B+
    if(choice == "nominal"){
      RooAddPdf model("model", "model",
                      RooArgList(*signal, *jpsipi, *fit_side, *erf),
                      RooArgList(*n_signal, *n_jpsipi, *n_combinatorial, n_nonprompt));
      RooAddPdf model_nonp("model_nonp","model",
                           RooArgList(*signal,*fit_side),
                           RooArgList(*n_signal,*n_combinatorial));
      if (include_np) {
        w.import(model, RecycleConflictNodes());
      } else {
        model_nonp.SetName("model");
        w.import(model_nonp, RecycleConflictNodes());
      }
      w.import(*lambda);
      w.import(*f_erf);
    } else if (choice == "sig3gauss") {
      RooAddPdf model("model", "model",
                      RooArgList(*signal_triple, *jpsipi, *fit_side, *erf),
                      RooArgList(*n_signal, *n_jpsipi, *n_combinatorial, n_nonprompt));
      signal_triple->SetName("signal");
      sigma1->setMax(0.015);
      w.import(model, RecycleConflictNodes());
    w.import(*lambda);
    w.import(*f_erf);
    } else if(choice == "bkg_poly"){
      RooAddPdf model("model","model",RooArgList(*signal,*poly_bkg,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    } else if(choice == "bkg_range"){
      RooAddPdf model("model","model",RooArgList(*signal,*fit_side,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_jpsipi));
      w.import(model);
    }else if(choice == "sig1gauss"){
      RooAddPdf model("model","model",RooArgList(*signal3,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "gauss_crystal"){
      RooAddPdf model("model","model",RooArgList(*gauss_CB,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "crystal_ball"){
      RooAddPdf model("model","model",RooArgList(*CB1,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }else if(choice == "double_crystal"){
      RooAddPdf model("model","model",RooArgList(*sum_CB,*fit_side,*erf,*jpsipi),RooArgList(*n_signal,*n_combinatorial,*n_erf,*n_jpsipi));
      w.import(model);
    }
  }

  else if(particle == 1){//Bs 
    if(choice == "nominal"){
      RooAddPdf model("model","model",RooArgList(*signal,*fit_side),RooArgList(*n_signal,*n_combinatorial)); 
      w.import(model);
      w.import(*lambda);
      w.import(*f_erf);
    }else if(choice == "bkg_poly"){
      RooAddPdf model("model","model",RooArgList(*signal,*poly_bkg),RooArgList(*n_signal,*n_combinatorial));
      w.import(model);
    }else if(choice == "gauss"){
      RooAddPdf model("model","model",RooArgList(*signal1,*fit_side),RooArgList(*n_signal,*n_combinatorial));
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
        // RooAddPdf model("model","model",RooArgList(*sum_cb_rt_sf,*double_CB_wt_sf,*fit_side),RooArgList(*RT_yield,*WT_yield,*n_combinatorial));
        // w.import(model);
      }
    }
    else if(MC == 1){
      if(component == 0){ //RT
        RooAddPdf model("model","model",RooArgList(*sum_CB),RooArgList(*n_signal));
        w.import(model);
      }
      else if(component == 1){ //WT
        // RooAddPdf model("model", "model",RooArgList(*double_CB_swp),RooArgList(*n_signal_swp));
        // w.import(model);
      }
    }
  }
} 
//build_pdf ends



/**
   Fit the non-prompt J/psi MC sample
   The shapes of erf, signal and J/psi pi peak are determined
   by their individual components.
   The relative yields and widths between signal and J/psi pi are fixed
   during successive fits
 */
void fit_jpsinp(RooWorkspace& w, std::string choice, const RooArgSet &c_vars,
                int ipt, int iy, bool inclusive, bool includeSignal) {
  int pti = ptlist[ipt];
  int ptf = ptlist[ipt + 1];

  RooAbsPdf*  model_cont = w.pdf("m_jpsinp_cont");
  RooDataSet* fullds = (RooDataSet*) w.data("jpsinp");

  RooDataSet* ds = fullds;
  // get B -> jpsi pi dataset before applying pT cut
  RooDataSet* ds_jpsipi = (RooDataSet*) ds->reduce("Bgen == 23335");
  // Apply pT cuts for all the other datasets
  if (!inclusive) {
    ds = (RooDataSet*) fullds->
      reduce(TString::Format("Bpt > %f && Bpt < %f", (double) pti, (double) ptf));
  }
  // Get rid of B+ at gen level
  RooDataSet* ds_cont = (RooDataSet*) ds->reduce("Bgen != 23333 && Bgen != 23335");
  RooDataSet* ds_sig = (RooDataSet*) ds->reduce("Bgen == 23333");

  std::vector<double> n_signal_initial = {3e3, 5e3, 5e3, 5000, 200};
  std::vector<double> n_cont_initial = {3e3, 1e4, 1.5e4, 1e3, 1e3};
  std::vector<double> n_erf_initial = {800, 1e3, 6e3, 1e3, 1e3};
  RooAbsPdf* signal = w.pdf("signal");
  RooAbsPdf* jpsipi = w.pdf("jpsipi");
  RooAbsPdf* erf = w.pdf("erf");
  RooAbsPdf* poly = w.pdf("poly_jpsi");
  RooRealVar np_sig_fraction("np_sig_fraction", "np_sig_fraction", 0.1, 0, 1);
  RooRealVar n_signal("n_signal_np", "n_signal_np",
                      n_signal_initial[ipt], 0., (ds->sumEntries())*2);
  RooRealVar n_cont("n_cont_np", "n_cont_np",
                    n_cont_initial[ipt], 0., (ds->sumEntries())*2);
  RooRealVar n_erf("n_nonprompt", "n_nonprompt",
                    n_erf_initial[ipt], 0., (ds->sumEntries())*2);

  TString ystr = "";
  // get relative yields of Jpsi pi to signal
  RooRealVar n_jpsipi_ext("n_jpsipi_ext", "n_jpsipi_ext",
                          0.1 * n_erf.getVal() , 0., (ds->sumEntries())*2);
  RooExtendPdf jpsipi_ext("jpsipi_ext", "extended jpsipi", *jpsipi, n_jpsipi_ext);
  w.import(n_jpsipi_ext);
  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* signal_width = w.var("sigma1");
  RooRealVar* m_jpsipi_width = w.var("m_jpsipi_width");
  RooRealVar* m_jpsipi_mean1 = w.var("m_jpsipi_mean1");
  // Fit jpsipi to MC truth samples
  // The width of jpsi pi is a constant relative to that of the signal's
  signal_width->setConstant();
  Bmass.setRange("bjpsipi", 5.22, 5.8);
  auto jpsipi_result = jpsipi_ext.fitTo(*ds_jpsipi, Range("bjpsipi"), Save(), Extended());
  auto jpsipi_par_list = jpsipi_result->floatParsFinal();
  TString jpsipiPlot = "./results/Bu/" +
    TString::Format("%i_%i/np_gen_jpsipi_pt_all%s.pdf",
                    pti, ptf, ystr.Data());
  if (inclusive) {
    cout << "make inclusive plot" << "\n";

    gSystem->MakeDirectory("./results/Bu/inclusive" );
    jpsipiPlot = "./results/Bu/inclusive/np_gen_jpsipi_all.pdf";
  }
  plot_mcfit(w, &jpsipi_ext, ds_jpsipi, jpsipiPlot,
             "NP gen-matched B^{+} #rightarrow J/#psi #pi^{+}",
             RooFit::Name("MCFit"), NormRange("bjpsipi"), LineColor(kRed),
             LineStyle(1), LineWidth(2));
  // fix the shape
  fix_parameters(w, jpsipi_par_list);
  double jpsipi_width_allpt = signal_width->getVal();
  signal_width->setConstant(false);
  m_jpsipi_mean1->setConstant(false);

  // jpsi pi fit with pT selection
  ds_jpsipi = (RooDataSet*) ds->reduce("Bgen == 23335");
  Bmass.setRange("bjpsipi_narrow", 5.25, 5.7);
  jpsipi_result = jpsipi_ext.fitTo(*ds_jpsipi, Range("bjpsipi_narrow"), Save(), Extended());
  jpsipi_par_list = jpsipi_result->floatParsFinal();
  jpsipiPlot = "./results/Bu/" +
    TString::Format("%i_%i/np_gen_jpsipi_pt%i-%i%s.pdf",
                    pti, ptf, pti, ptf, ystr.Data());
  if (inclusive) {
    jpsipiPlot = "./results/Bu/inclusive/np_gen_jpsipi.pdf";
  }
  plot_mcfit(w, &jpsipi_ext, ds_jpsipi, jpsipiPlot,
             "NP gen-matched B^{+} #rightarrow J/#psi #pi^{+}",
             RooFit::Name("MCFit"), NormRange("bjpsipi"), LineColor(kRed),
             LineStyle(1), LineWidth(2));
  double jpsipi_width_val = signal_width->getVal() / jpsipi_width_allpt;

  // NP background fit with no B+ or J/psi pi component
  auto cont_result = model_cont->fitTo(*ds_cont, Save());
  auto cont_par_list = cont_result->floatParsFinal();
  TString jpsi_fit_plot = "./results/Bu/" +
    TString::Format("%i_%i/np_fit_pt%i-%i%s.pdf", pti, ptf, pti, ptf, ystr.Data());
  if (inclusive) {
    jpsi_fit_plot = "./results/Bu/inclusive/np_fit.pdf";
  }
  plot_jpsifit(w, model_cont, ds_cont, jpsi_fit_plot, "Non-prompt J/#psi", 0, n_signal);
  fix_parameters(w, cont_par_list);

  if (!includeSignal) {
    return;
  }
  // Fit pure signal
  Bmass.setRange("bmc", 5.18, 5.38);
  RooExtendPdf signal_ext("signal_ext", "extended signal pdf", *signal, n_signal);
  auto signal_result = signal_ext.fitTo(*ds_sig, Range("bmc"), Save(), Extended());
  auto signal_par_list = signal_result->floatParsFinal();

  TString signalPlot = "./results/Bu/" +
    TString::Format("%i_%i/np_gen_signal_pt%i-%i%s.pdf",
                    pti, ptf, pti, ptf, ystr.Data());
  if (inclusive) {
    signalPlot = "./results/Bu/inclusive/np_gen_signal.pdf";
  }
  plot_mcfit(w, &signal_ext, ds_sig, signalPlot, "NP gen-matched signal",
             RooFit::Name("MCFit"), Range("bmc"), LineColor(kRed),
             LineStyle(1), LineWidth(2));


  // Fix the ratio of jpsipi to signal
  RooRealVar* jpsipi_to_signal_ratio = w.var("jpsipi_to_signal_ratio");
  RooRealVar* jpsinp_poly_fraction = w.var("jpsinp_poly_fraction");
  RooRealVar* jpsipi_fraction = w.var("jpsinp_jpsipi_fraction");
  jpsipi_to_signal_ratio->setVal(n_jpsipi_ext.getVal() / n_signal.getVal());
  jpsipi_to_signal_ratio->setConstant();

  RooProduct* n_jpsipi = new RooProduct("n_jpsipi_by_signal",
                                        "number of jpsi pi with fixed ratio to n_signal",
                                        RooArgList(n_signal, *jpsipi_to_signal_ratio));
  cout << n_jpsipi->getVal() << "\n";

  // Fix the width of jpsipi
  RooRealVar* jpsipi_to_signal_width_ratio = w.var("jpsipi_to_signal_width_ratio");
  jpsipi_to_signal_width_ratio->setVal(jpsipi_width_val / signal_width->getVal());

  // include signal pdf
  RooAddPdf* model_inclusive = new
    RooAddPdf("np_signal", "NP with B+",
              RooArgList(*signal, *jpsipi, *erf, *poly),
              RooArgList(n_signal, *n_jpsipi, n_erf, n_cont));
  // fit the J/psi using fixed signal and jpsipi shape
  fix_signal_shape(w);
  // fix the erf shape
  fix_parameters(w, "erf");
  // fix the jpsi pi mean
  m_jpsipi_mean1->setConstant();
  model_inclusive->fitTo(*ds, Save(), Extended(), NumCPU(4));
  TString jpsi_plot_with_sig = "./results/Bu/" +
    TString::Format("%i_%i/np_fit_signal_pt%i-%i%s.pdf",
                    pti, ptf, pti, ptf, ystr.Data());
  if (inclusive) {
    jpsi_plot_with_sig = "./results/Bu/inclusive/np_fit_signal.pdf";
  }
  plot_jpsifit(w, model_inclusive, ds, jpsi_plot_with_sig, "Non-prompt J/#psi with signal",
               ds_sig->sumEntries(), n_signal);

  // float the signal parameters after the fit
  fix_signal_shape(w, true);

  cout << "MC fit for NP J/psi complete" << "\n";
}

void plot_complete_fit(RooWorkspace& w, RooArgSet &c_vars, TString subname, int iy=1){
cout <<"ploting complete fit"<< endl;
  RooAbsPdf*  model = w.pdf("model");
  RooAbsPdf*  signal = w.pdf("signal");
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooDataSet* mc = (RooDataSet*) w.data("mc");

  TString ystr = "";

  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* lambda   = w.var("lambda");
  RooRealVar* n_signal   = w.var("n_signal");
  RooAddPdf modelmc("modelmc","modelmc",RooArgList(*signal), RooArgList(*n_signal));

  if( (particle == 2) && (MC == 0) ){
    model->fitTo(*data,Range("all"),Constrain(c_vars),Extended(kTRUE));
  } else{
    // determine gaussian width with MC

    if (particle == 0) {
      Bmass.setRange("bmc", 5.18, 5.38);
    } else if (particle == 1) {
      Bmass.setRange("bmc", 5.26, 5.46);
    }
    // modelmc.fitTo(*mc, Range("bmc"), Extended(kTRUE));
    fix_signal_shape(w, true);
    auto signal_result = signal->fitTo(*mc, Range("bmc"), Save());
    auto signal_par_list = signal_result->floatParsFinal();
    cout << "MC fit complete" << "\n";

    TString signalPlot;
    if (particle == 0) {
      signalPlot = "./results/Bu/" + subname + "/mc_fit_Bu.pdf";
    } else if (particle == 1) {
      signalPlot = "./results/Bs/" + subname + "/mc_fit_Bs.pdf";
    }
    plot_mcfit(w, signal, mc, signalPlot, "MC signal",
               RooFit::Name("MCFit"), NormRange("bmc"), LineColor(kRed),
               LineStyle(1), LineWidth(2));
    fix_signal_shape(w);
    // float the mean of signal for data
    RooRealVar* signal_mean = w.var("mean");
    signal_mean->setConstant(false);

  }

  TFile* f;
  if(MC == 1){
    if(component == 0){f = new TFile("./results/B0/MC/RT_fit.root", "RECREATE");}
    else if(component == 1){f = new TFile("./results/B0/MC/WT_fit.root", "RECREATE");}
  }
  else if(MC == 0){f = new TFile("./results/" + BDTdir[particle] + "/" + subname + "/DATA_fit.root", "RECREATE");}

  RooFitResult* r = model->fitTo(*data,Range("all"),Save(), NumCPU(4));
  r->Print();
  f->cd();
  r->Write();
  f->Close();

  RooPlot* massframe = Bmass.frame();
  if(particle == 0){
    data->plotOn(massframe, RooFit::Name("Data"), MarkerSize(0.9));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),NormRange("all"),LineColor(kOrange-3),LineStyle(kDashed),LineWidth(3), FillStyle(3002),FillColor(kOrange-3),VLines(),DrawOption("LF")); 
    model->plotOn(massframe, RooFit::Name("Fit"),NormRange("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),NormRange("all"),LineColor(kBlue),LineStyle(kDashed));

    model->plotOn(massframe, RooFit::Name("B->J/psi X"),Components("erf"),NormRange("all"),LineColor(kGreen+3),LineStyle(1),LineWidth(3),FillStyle(3005),FillColor(kGreen+3),VLines(),DrawOption("LF"));
    model->plotOn(massframe, RooFit::Name("B->J/psi pi"),Components("jpsipi"),NormRange("all"),LineColor(kPink+10),LineStyle(kDashed));
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("m(#mu^{+}#mu^{-}K^{+}) [GeV]");
    model->paramOn(massframe,  Layout(0.35, 0.6, 0.9),
                   Format("NEU", AutoPrecision(1)));
  }
  else if(particle == 1){
    data->plotOn(massframe, RooFit::Name("Data"), MarkerSize(0.9));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),NormRange("all"),LineColor(kOrange-3),LineStyle(kDashed),LineWidth(3), FillStyle(3002),FillColor(kOrange-3),VLines(),DrawOption("LF")); 
    model->plotOn(massframe, RooFit::Name("Fit"),NormRange("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),NormRange("all"),LineColor(kBlue),LineStyle(kDashed));
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("m(#mu^{+}#mu^{-}K^{+}K^{-}) [GeV]");
  }  
  else if(particle == 2){
    data->plotOn(massframe, RooFit::Name("Data"), MarkerSize(0.9));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),NormRange("all"),LineColor(kOrange-3),LineStyle(kDashed),LineWidth(3), FillStyle(3002),FillColor(kOrange-3),VLines(),DrawOption("LF")); 
    model->plotOn(massframe, RooFit::Name("Fit"),NormRange("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),NormRange("all"),LineColor(kBlue),LineStyle(kDashed));
      if(MC == 0){
        model->plotOn(massframe, RooFit::Name("Corr Tag"), RooFit::Components("rt_pdf"), RooFit::NormRange("all"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
        model->plotOn(massframe, RooFit::Name("Mis Tag"), RooFit::Components("wt_pdf"), RooFit::NormRange("all"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));  
        model->plotOn(massframe, RooFit::Name("Combinatorial"), RooFit::Components("fit_side"), RooFit::NormRange("all"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
        model->plotOn(massframe, RooFit::Name("B->J/#psi X"),Components("erf"),NormRange("all"),LineColor(kCyan),LineStyle(kDashed));
        massframe->SetXTitle("m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]");
      }
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
  
  massframe->SetTitle(" ");
  massframe->Draw();

  
  TLatex* tex11 = new TLatex(0.66,0.92,"302.3 pb^{-1} (pp) 5.02 TeV");
  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  tex11->Draw();
  tex11 = new TLatex(0.1,0.92,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  tex11->Draw();
 /* 
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
  TLegend *leg = new TLegend (0.7, 0.45, 0.9, 0.85);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(massframe->findObject("Data"), "Data", "p)");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "f");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");

  if(particle == 0){
    leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/#psi X", "l");
    leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/#psi #pi^{+}", "l");
  }
  else if(particle == 2){
    if(MC == 0){
      leg->AddEntry(massframe->findObject("Total Signal"), "Total sig", "l");
      leg->AddEntry(massframe->findObject("Corr Tag"), "Corr. sig", "l");
      leg->AddEntry(massframe->findObject("Mis Tag"), "Mis-tag sig", "l");
      leg->AddEntry(massframe->findObject("Combinatorial"), "Comb. bkg", "l");
    }
  }
  leg->SetBorderSize(0);
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
    d.SaveAs("./results/Bu/" + subname + "/complete_fit_Bu" + ystr + ".pdf");
  }
  else if(particle == 1){
    d.SaveAs("./results/Bs/" + subname + "/complete_fit_Bs" + ystr + ".pdf");
  }
  else if(particle == 2){
    if(MC == 0){
      d.SaveAs("./results/B0/DATA_fit_B0.pdf");
    }
    else if(MC == 1){
      if(component == 0){
        d.SaveAs("./results/B0/MC_RT_fit_B0.pdf");
      }
      else if(component == 1){
        d.SaveAs("./results/B0/MC_WT_fit_B0.pdf");
      }
    }
  }
}
//plot_complete_fit ends

//SIDEBAND SUBTRACTION//
std::vector<TH1D*> sideband_subtraction(RooWorkspace& w, std::vector<TString> label, int n_var) {
  
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooAbsPdf* BgModel = w.pdf("fit_side");
  RooRealVar Bmass = *(w.var("Bmass"));
  
  vector<RooRealVar> variables;
  variables.push_back(*(w.var("Bmass"))); 

  for(int i=0; i<n_var; i++){variables.push_back(*(w.var(label[i])));}

  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central; 

  double left;
  double right;
  double max, min;
  double mean = w.var("mean")->getValV();
  double sigma = 0.015;

  if(particle == 0){
    left = mean - 3 * sigma;
    right = mean + 3 * sigma;
    min = mean - 9 * sigma;
    max = mean + 9 * sigma;
  }
  else if(particle == 1){
    left = mean - 3 * sigma;
    right = mean + 3 * sigma;
    min = mean - 9 * sigma;
    max = mean + 9 * sigma;
  }
  else if(particle == 2){
    left = 5.22;
    right = 5.32;
  }

  Bmass.setRange("right", right, max);
  Bmass.setRange("left", min, left);
  Bmass.setRange("peak", left, right);

  if(particle == 0){reduceddata_side = (RooDataSet*)data->reduce(Form("Bmass>%lf", right))->reduce(Form("Bmass < %lf", max));}   //only events w bigger mass than the peak ?partlialy recosntructed  background ? 
  else if( (particle == 1) || (particle == 2) ){
    reduceddata_side =  (RooDataSet*)data->reduce(Form("Bmass>%lf || Bmass<%lf", right, left))->reduce(Form("Bmass > %lf && Bmass < %lf", min, max));
  }

  reduceddata_central = (RooDataSet*)data->reduce(Form("Bmass>%lf",left));
  reduceddata_central = (RooDataSet*)reduceddata_central->reduce(Form("Bmass<%lf",right));


  //Integrating the background distribution 
  RooAbsReal* int_fit_side_right = BgModel->createIntegral(Bmass, Bmass, "right");
  RooAbsReal* int_fit_side_left = BgModel->createIntegral(Bmass, Bmass, "left");
  RooAbsReal* int_fit_peak = BgModel->createIntegral(Bmass, Bmass, "peak");
 
  cout << "Integral left band = " << int_fit_side_left->getVal() << endl; 
  cout << "Integral right band = " << int_fit_side_right->getVal() << endl;
  cout << "Integral peak = " << int_fit_peak->getVal() << endl; 

  double factor;
  if (particle == 0) {
    factor = (int_fit_peak->getVal()) / (int_fit_side_right->getVal());
  } else if ( (particle == 1) || (particle == 2) ) {
    factor = (int_fit_peak->getVal()) /
      (int_fit_side_right->getVal() + int_fit_side_left->getVal());
  }
  std::cout << std::endl << "Factor: " << factor << std::endl;

  /*for(int i=0; i<n_var; i++){
    std::cout << "bins: " << n[i] << std::endl;
  }*/
 
  std::vector<TH1D*> histos;

  for(int i=1; i<(n_var+1);i++){histos.push_back(create_histogram(variables[i], std::string(label[i-1]),factor, reduceddata_side, reduceddata_central, data, 40));}
  return histos; 

}
//sideband_subtraction ends

//create_histogram_mc starts
TH1D* create_histogram_mc(RooWorkspace& w, RooRealVar var, int n, TString weight){

  RooDataSet* mcSAMPLE = (RooDataSet*) w.data("mc");

  TH1D* h = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());

  if(std::string(var.GetName()) == "BsvpvDisErr"){          h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));}
  else if(std::string(var.GetName()) == "Bpt"){             h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 5, 60));}
  else if(std::string(var.GetName()) == "BsvpvDistance"){   h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 1));}
  else if(std::string(var.GetName()) == "Bmumumass"){       h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 3, 3.2));}
  else if(std::string(var.GetName()) == "Bchi2cl"){         h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0.05, 1));}
  else if(std::string(var.GetName()) == "Balpha"){          h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.1));}
  else if(std::string(var.GetName()) == "Bd0"){             h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.5));}
  else if(std::string(var.GetName()) == "Bd0Err"){          h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0., 0.0001));}
  else if(std::string(var.GetName()) == "Bdtheta"){         h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0.05, 0.1));}
  else if(std::string(var.GetName()) == "Bmu1dxyPV"){       h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Bmu1dzPV"){        h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));}
  else if(std::string(var.GetName()) == "Bmu1pt"){          h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 16));}
  else if(std::string(var.GetName()) == "Bmu2dxyPV"){       h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Bmu2dzPV"){        h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));}
  else if(std::string(var.GetName()) == "Bmu2pt"){          h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 16));}
  else if(std::string(var.GetName()) == "Bmumueta"){        h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -3.5, 3.5));}
  else if(std::string(var.GetName()) == "Bmumupt"){         h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 30));}
  else if(std::string(var.GetName()) == "BsvpvDisErr_2D"){  h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0.002, 0.012));}
  else if(std::string(var.GetName()) == "BsvpvDistance_2D"){h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 0.5));}
  else if(std::string(var.GetName()) == "Btrk1Dxy1"){       h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Btrk1DxyError1"){  h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));}
  else if(std::string(var.GetName()) == "Btrk1Dz1"){        h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, -0.25, 0.25));}
  else if(std::string(var.GetName()) == "Btrk1Pt"){         h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0 , 8));}
  else if(std::string(var.GetName()) == "Btrk1DzError1"){   h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0 ,0.05));}
  else if(std::string(var.GetName()) == "nMult"){           h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0, 100 ));}
  else if(std::string(var.GetName()) == "Btrk1PtErr"){      h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0, 0.1 ));}
  else if(std::string(var.GetName()) == "Btrk2PtErr"){      h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0, 0.1));}
  else if(std::string(var.GetName()) == "Btrk2Pt"){         h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n,0 , 8));}
  else if (indexBDT.count(var.GetName())) {
    auto iBDT = indexBDT.at(var.GetName());
    h = (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var,Binning(BDTnbins[iBDT],BDTmin[iBDT],BDTmax[iBDT]));
  }
  else{ h= (TH1D*) mcSAMPLE->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));}

  h->SetTitle("");
  //h->SetMarkerStyle(29);
  h->SetMarkerColor(kGreen);
  h->SetMarkerSize(1);
  h->SetLineColor(kGreen);
  return h;

}
//create_histogram_mc ends

//create_histogram  TO BE USED IN SIDEBAND
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n){

  TH1D* dist_side;
  TH1D* hist_dist_peak;
cout << endl;

  if(std::string(var.GetName()) == "BsvpvDisErr"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));
	dist_side      = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));}
  else if(std::string(var.GetName()) == "Bpt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 5, 60));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 5, 60));}
  else if(std::string(var.GetName()) == "BsvpvDistance"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 1));}
  else if(std::string(var.GetName()) == "Bmumumass"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 3, 3.2));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 3, 3.2));}
  else if(std::string(var.GetName()) == "Bchi2cl"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.05, 1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0.05, 1));}
  else if(std::string(var.GetName()) == "Balpha"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.1));}
  else if(std::string(var.GetName()) == "Bd0"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.5));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.5));}
  else if(std::string(var.GetName()) == "Bd0Err"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.0001));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0., 0.0001));}
  else if(std::string(var.GetName()) == "Bdtheta"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.0, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0.05, 0.1));}
  else if(std::string(var.GetName()) == "Bmu1dxyPV"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Bmu1dzPV"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));}
  else if(std::string(var.GetName()) == "Bmu1pt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0 , 16));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 16));}
  else if(std::string(var.GetName()) == "Bmu2dxyPV"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.1 , 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Bmu2dzPV"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.2, 0.2));}
  else if(std::string(var.GetName()) == "Bmu2pt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0 , 16));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 16));}
  else if(std::string(var.GetName()) == "Bmumueta"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -3.5 , 3.5));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -3.5, 3.5));}
  else if(std::string(var.GetName()) == "Bmumupt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 30));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 30));}
  else if(std::string(var.GetName()) == "BsvpvDisErr_2D"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0.002, 0.012));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0.002, 0.012));}
  else if(std::string(var.GetName()) == "BsvpvDistance_2D"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.5));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 0.5));}
  else if(std::string(var.GetName()) == "Btrk1Dxy1"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.1, 0.1));}
  else if(std::string(var.GetName()) == "Btrk1DxyError1"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, 0, 0.03));}
  else if(std::string(var.GetName()) == "Btrk1Dz1"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, -0.25, 0.25));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, -0.25, 0.25));}
  else if(std::string(var.GetName()) == "Btrk1Pt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n,0 , 8));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0 , 8));}
  else if(std::string(var.GetName()) == "Btrk1DzError1"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.05));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0 ,0.05));}
  else if(std::string(var.GetName()) == "nMult"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 100));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0, 100 ));}
  else if(std::string(var.GetName()) == "Btrk1PtErr"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0, 0.1 ));}
  else if(std::string(var.GetName()) == "Btrk2PtErr"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, 0, 0.1));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0, 0.1));}
  else if(std::string(var.GetName()) == "Btrk2Pt"){
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n,0 , 8));
	dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n,0 , 8));}
  else if (indexBDT.count(var.GetName())) {
    auto iBDT = indexBDT.at(var.GetName());
    hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var,
                                                      Binning(BDTnbins[iBDT],
                                                              BDTmin[iBDT],
                                                              BDTmax[iBDT]));
    dist_side = (TH1D*) reduced->createHistogram(var.GetName(), var,
                                                 Binning(BDTnbins[iBDT],
                                                         BDTmin[iBDT],
                                                         BDTmax[iBDT]));
  }
  else{
	hist_dist_peak = (TH1D*) central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
	dist_side      = (TH1D*) reduced->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));}

cout << endl;

  dist_side->SetMarkerColor(kRed);
  dist_side->SetLineColor(kRed);
  dist_side->SetNameTitle("dist_side", "");

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

  hist_dist_peak->Draw();  //signal + background
  dist_side->Draw("same"); //background
  dist_peak->Draw("same");
  
  dist_peak->SetXTitle(var.GetName());
  dist_side->SetXTitle(var.GetName());
  hist_dist_peak->SetXTitle(var.GetName());


hist_dist_peak->GetYaxis()->SetRangeUser(0, 1.3*hist_dist_peak->GetMaximum());
  TLatex* tex = new TLatex(0.6,0.8,"302.3 pb^{-1} (pp) 5.02 TeV");
  /*tex->SetNDC(kTRUE);
 *   tex->SetLineWidth(2);
 *     tex->SetTextSize(0.04);
 *       tex->Draw();
 *         tex = new TLatex(0.68,0.85,"CMS Preliminary");
 *           tex->SetNDC(kTRUE);
 *             tex->SetTextFont(42);
 *               tex->SetTextSize(0.04);
 *                 tex->SetLineWidth(2);
 *                   tex->Draw();*/

TLegend *leg = new TLegend (0.7, 0.9, 0.9, 1.0);
  leg->AddEntry("dist_peak", "Signal", "LE");
  leg->AddEntry("dist_side", "Background", "LE");
  leg->AddEntry("hist_dist_peak", "Total", "LE");
  leg->Draw("same");

  if(particle == 0){
    gSystem->Exec("mkdir -p ./results/Bu/sideband_sub/");
    c.SaveAs("./results/Bu/sideband_sub/"+ name + "sideband_sub_Bu.pdf");
    }else if(particle == 1){
    gSystem->Exec("mkdir -p ./results/Bs/sideband_sub/");
    c.SaveAs("./results/Bs/sideband_sub/"+ name + "sideband_sub_Bs.pdf");
    }else if(particle == 2){
    c.SaveAs("./results/B0/sideband_sub/"+ name + "sideband_sub_B0.pdf");
    }
   
    if(background == 0){return dist_peak;}
    else if(background == 1){return dist_side;}
    }
//create_histogram  TO BE USED IN SIDEBAND

//do_splot
void do_splot(RooWorkspace& w, RooArgSet &c_vars){

  RooDataSet* data = (RooDataSet*) w.data("data");   
  RooAbsPdf* model = w.pdf("model");
  //we need the fit and the dataset previously saved in the woorkspace

  RooRealVar* BpYield = 0;
  if( (particle == 2) && (MC == 0) ){BpYield = w.var("RT_yield");}
  else{BpYield = w.var("n_signal");}
  RooRealVar* combYield = w.var("n_combinatorial");
  RooRealVar* npYield = w.var("n_nonprompt");
  //we need the n values previously saved in the woorkspace

  //fit the model to the data
  if( (particle == 2) && (MC == 0) ){model->fitTo(*data,Extended(),Constrain(c_vars));}
  else{model->fitTo(*data,Extended());}

  //sPlot technique requires model parameters (other than the yields) to be fixed
  RooRealVar* alpha1;
  RooRealVar* alpha2;
  RooRealVar* n1;
  RooRealVar* n2;
  RooRealVar* mean_difference;
  RooRealVar* sigma1_swp;
  RooRealVar* sigma_cb1;
  RooRealVar* sigma_cb2;
  RooRealVar* alpha1_swp;
  RooRealVar* alpha2_swp;
  RooRealVar* n1_swp;
  RooRealVar* n2_swp;
  RooRealVar* sigma1;
  RooRealVar* sigma2;

  RooRealVar* mean = w.var("mean");
  RooRealVar* cofs = w.var("cofs");
  RooRealVar* lambda = w.var("lambda") ;
  RooRealVar* slope = w.var("slope") ;
  mean->setConstant();
  cofs->setConstant();
  lambda->setConstant();
  // slope->setConstant();
  RooRealVar* f_erf = w.var("f_erf") ;
  f_erf->setConstant();

  fix_signal_shape(w);
  fix_parameters(w, "fit_side");
  if(particle != 2){  //both Bu and Bs nominal model are a double gaussian
  sigma1 = w.var("sigma1");
    sigma1->setConstant();
  }
  else{
  sigma1_swp = w.var("sigma1_swp");
  sigma_cb1 = w.var("sigma_cb1");
  sigma_cb2 = w.var("sigma_cb2");
  mean_difference = w.var("mean_difference");
  n1 = w.var("n1");
  n2 = w.var("n2");
  n1_swp = w.var("n1_swp");
  n2_swp = w.var("n2_swp");
  alpha1 = w.var("alpha1");
  alpha2 = w.var("alpha2");  
  alpha1_swp = w.var("alpha1_swp");  
  alpha2_swp = w.var("alpha2_swp"); 
    alpha1->setConstant();
    alpha1_swp->setConstant(); 
    alpha2->setConstant();
    alpha2_swp->setConstant(); 
    sigma1_swp->setConstant();
    sigma_cb1->setConstant();
    sigma_cb2->setConstant();
    n1->setConstant();
    n1_swp->setConstant();
    n2->setConstant();
    n2_swp->setConstant();
    mean_difference->setConstant();}

  // cout << "do check fit" << "\n";
  // model->fitTo(*data, Extended());
  // cout << "finished check fit" << "\n";
  cout << "before splot" << "\n";
  RooMsgService::instance().setSilentMode(true);

  //add sWeights to dataset based on model and yield variables
  //sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
  // RooArgList pdfComp(*BpYield, *combYield);
  RooArgList pdfComp(*BpYield, *combYield);
  // For B+, include non-prompt J/psi component
  if (particle == 0) {pdfComp.add(*npYield);}

  SPlot* sData = new SPlot("sData","An sPlot",*data, model, pdfComp);





    // Get the sWeights
    const RooArgList& sWeights = sData->GetSWeightVars();

    // Loop over the dataset and print the sWeights for each event
    for (int i = 0; i < 100; ++i) {
        const RooArgSet* event = data->get(i);
        const RooRealVar* sWeightVar = dynamic_cast<const RooRealVar*>(sWeights.at(0)); // Assuming there's only one sWeight variable

        if (sWeightVar) {
            double sWeight = event->getRealValue(sWeightVar->GetName());
            cout << "xxxx Event " << i << ": sWeight = " << sWeight << endl;
        } else {
            cerr << "Error: Unable to access sWeight variable." << endl;
        }
    }






  cout << endl <<  "Yield of B+ is "
       << BpYield->getVal() << ".  From sWeights it is "
       // << sData->GetYieldFromSWeight("RT_yield") << endl;
       << sData->GetYieldFromSWeight("n_signal") << endl;

  cout << "Yield of background is "
       << combYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;
  
  w.import(*data, Rename("dataWithSWeights"));
  //the reweighted data is saved in the woorkspace 
}
//do_splot ends

//make_splot
TH1D* make_splot(RooWorkspace& w, int n, TString label){

  //saves the plots of signal distributions, background distributions and signal+background distributions
  //in the end returns the histogram of signal

  TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
  cdata->Divide(2,2);

  RooAbsPdf* model  = w.pdf("model");
  RooAbsPdf* BpModel;

  if( (particle == 2) && (MC == 0) ){BpModel = w.pdf("rt_pdf");}
  else if(particle == 0){BpModel = w.pdf("signal");}
  else if(particle == 1){BpModel = w.pdf("signal");}
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
  model->plotOn(mframe,NormRange("all"), LineColor(kRed));
  model->plotOn(mframe,NormRange("all"), Components(*BpModel), LineStyle(kDashed), LineColor(kOrange));
  model->plotOn(mframe,NormRange("all"), Components(*BgModel), LineStyle(kDashed), LineColor(kBlue));

  model->paramOn(mframe,Layout(.060,0.99,0.95));

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

  if(label == "BsvpvDisErr"){
	 ptframe2Bp = variable->frame(0.0,0.03,40);
	 ptframe2Bg = variable->frame(0.0,0.03,40);}

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
    gSystem->Exec("mkdir -p ./results/Bu/splot/Bmass/");
    cdata->SaveAs("./results/Bu/splot/Bmass/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    gSystem->Exec("mkdir -p ./results/Bs/splot/Bmass/");
    cdata->SaveAs("./results/Bs/splot/Bmass/"+label+"sPlot_Bs.pdf");
  }else if(particle == 2){
    cdata->SaveAs("./results/B0/splot/Bmass/"+label+"sPlot_B0.pdf");
  }


  TH1D* histo_Bp_sig ;
  TH1D* histo_Bp_bkg;

  if(label == "BsvpvDisErr"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.03));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.03));}
  else if(label == "BsvpvDistance"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 1));}
  else if(label == "Bpt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 5, 60));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 5, 60));}
  else if(label == "Bmumumass"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 3, 3.2));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 3, 3.2));}
  else if(label == "Bchi2cl"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0.05, 1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0.05, 1));}
  else if(label == "Balpha"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0., 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0., 0.1));}
  else if(label == "Bd0"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0., 0.5));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0., 0.5));}
  else if(label == "Bd0Err"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0., 0.0001));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0., 0.0001));}
  else if(label == "Bdtheta"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0., 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0., 0.1));}
  else if(label == "Bmu1dxyPV"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.1, 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.1, 0.1));}
  else if(label == "Bmu1dzPV"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.2, 0.2));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.2, 0.2));}
  else if(label == "Bmu1pt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 16));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 16));}
  else if(label == "Bmu2dxyPV"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.1, 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.1, 0.1));}
  else if(label == "Bmu2dzPV"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.2, 0.2));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.2, 0.2));}
  else if(label == "Bmu2pt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 16));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 16));}
  else if(label == "Bmumueta"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -3.5, 3.5));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -3.5, 3.5));}
  else if(label == "Bmumupt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 30));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 30));}
  else if(label == "BsvpvDisErr_2D"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0.002, 0.012));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0.002, 0.012 ));}
  else if(label == "BsvpvDistance_2D"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.5));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.5 ));}
  else if(label == "Btrk1Dxy1"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.1, 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.1, 0.1 ));}
  else if(label == "Btrk1DxyError1"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.03));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.03 ));}
  else if(label == "Btrk1Dz1"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, -0.25, 0.25));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, -0.25, 0.25 ));}
  else if(label == "Btrk1Pt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 8));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 8));}
  else if(label == "Btrk1DzError1"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.05));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.05));}
  else if(label == "nMult"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 100));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 100));}
  else if(label == "Btrk1PtErr"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.1));}
  else if(label == "Btrk2PtErr"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 0.1));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 0.1));}
  else if(label == "Btrk2Pt"){
         histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,Binning(40, 0, 8));
 	 histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,Binning(40, 0, 8));}
  else if (indexBDT.count(label.Data())) {
    auto iBDT = indexBDT.at(label.Data());
    histo_Bp_sig = (TH1D*) dataWBp->createHistogram(label,*variable,
                                                    Binning(BDTnbins[iBDT],
                                                            BDTmin[iBDT],
                                                            BDTmax[iBDT]));
    histo_Bp_bkg = (TH1D*) dataWBg->createHistogram(label,*variable,
                                                    Binning(BDTnbins[iBDT],
                                                            BDTmin[iBDT],
                                                            BDTmax[iBDT]));
  }
  else{
    cout << "Create histogram " << label << " using variable range" << "\n";
   	histo_Bp_sig = (TH1D*)dataWBp->createHistogram(label,n,variable->getMin(),variable->getMax());
   	histo_Bp_bkg = (TH1D*)dataWBg->createHistogram(label,n,variable->getMin(),variable->getMax());}
  
  for (int i=1; i<=n; i++) {
    if (histo_Bp_sig->GetBinContent(i)==0) histo_Bp_sig->SetBinError(i,0.);
    if (histo_Bp_bkg->GetBinContent(i)==0) histo_Bp_bkg->SetBinError(i,0.);
     histo_Bp_sig->SetBinContent(i,histo_Bp_sig->GetBinContent(i)/sigYield);
     histo_Bp_sig->SetBinError(i,histo_Bp_sig->GetBinError(i)/sigYield);
     histo_Bp_bkg->SetBinContent(i,histo_Bp_bkg->GetBinContent(i)/bkgYield);
     histo_Bp_bkg->SetBinError(i,histo_Bp_bkg->GetBinError(i)/bkgYield);}

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
    prov->SaveAs("./results/Bu/splot/sig/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    prov->SaveAs("./results/Bs/splot/sig/"+label+"sPlot_Bs.pdf");
  }else if(particle == 2){
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
    prov_bkg->SaveAs("./results/Bu/splot/bkg/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    prov_bkg->SaveAs("./results/Bs/splot/bkg/"+label+"sPlot_Bs.pdf");
  }else if(particle == 2){
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
    gSystem->Exec("mkdir -p ./results/Bu/splot/sig_bkg");
    sig_bkg->SaveAs("./results/Bu/splot/sig_bkg/"+label+"sPlot_Bu.pdf");
  }else if(particle == 1){
    gSystem->Exec("mkdir -p ./results/Bs/splot/sig_bkg");
    sig_bkg->SaveAs("./results/Bs/splot/sig_bkg/"+label+"sPlot_Bs.pdf");
  }else if(particle == 2){
    sig_bkg->SaveAs("./results/B0/splot/sig_bkg/"+label+"sPlot_B0.pdf");
  }

  //cleanup
  delete cdata;
  delete prov;
  delete prov_bkg;
  delete sig_bkg;
  
  if (background == 0){return histo_Bp_sig;}
  else if(background == 1){return histo_Bp_bkg;}
} 
//make_splot ends

//SPLOT_METHOD//
std::vector<TH1D*> splot_method(RooWorkspace& w, std::vector<TString> label, int n_var){

  std::vector<TH1D*> histos;
  for(int i = 0;i<n_var;i++){histos.push_back(make_splot(w,40,label[i]));}

  return histos;
}


//splot_method endss

void validate_fit(RooWorkspace* w, RooArgSet &c_vars){
  RooRealVar Bmass = *(w->var("Bmass"));
  RooAbsPdf* model  = w->pdf("model");
  
  vector<RooRealVar> params;
  if( (particle == 2) && (MC == 0) ){params.push_back(*(w->var("RT_yield")));}
  else{params.push_back(*(w->var("n_signal")));}

  int params_size = params.size();  

  RooMCStudy* mcstudy;
  if( (particle == 2) && (MC == 0) ){mcstudy = new RooMCStudy(*model, Bmass, Constrain(c_vars), Binned(kTRUE), Silence(), Extended(kTRUE), FitOptions(Save(kTRUE), PrintEvalErrors(0)));}
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
    c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.pdf");
  }else if(particle == 1){
    c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.pdf");
    c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.pdf");
  }else if(particle == 2){
    c_pull->SaveAs("./results/B0/pulls/pulls_poisson_B0.pdf");
    c_params->SaveAs("./results/B0/pulls/pulls_params_poisson_B0.pdf");
  }  
}



void AddWeights(TTree* t){
  TString input_file;
  input_file = particle ? "./results/Bs/mc_validation_plots/weights/weights.root": "./results/Bu/mc_validation_plots/weights/weights.root";
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
  Float_t nMult;
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
  t->SetBranchAddress("nMult", &nMult);

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

void set_up_workspace_variables(RooWorkspace& w){

  if(particle == 0){
    float mass_min, mass_max;
    float y_min, y_max;
    float pt_min, pt_max;
    float nMult_min, nMult_max;
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
    //float BDT_1_2_min, BDT_1_2_max;
    //float BDT_2_3_min, BDT_2_3_max;
    float BDT_3_5_min, BDT_3_5_max;
    float BDT_5_7_min, BDT_5_7_max;
    float BDT_7_10_min, BDT_7_10_max;
    float BDT_10_15_min, BDT_10_15_max;
    float BDT_15_20_min, BDT_15_20_max;
    float BDT_20_50_min, BDT_20_50_max;
   
    mass_min = 5.0;
    mass_max = 6.0;
    
    y_min = -2.6;
    y_max = 2.6;   //2.6

    pt_min = 0.;  //0 -> 100
    pt_max = 100.;

    nMult_min = 0;
    nMult_max = 130;

    trk1pt_min = 0.;
    trk1pt_max = 30.;

    trk1eta_min = -3;
    trk1eta_max = 3;

    trk1pterr_min = 0.0;  
    trk1pterr_max = 0.5; 

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 40.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = .15;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 10;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.08;

    mumumass_min = 2.9;
    mumumass_max = 3.5; 
 
    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 0.;
    mu1pt_max = 80.;

    mu2pt_min = 0.;
    mu2pt_max = 80.;

    mu1dxyPV_min = -3; 
    mu1dxyPV_max = 3;

    mu2dxyPV_min = -0.5;
    mu2dxyPV_max = 0.5;

    mu1dzPV_min = -20.; 
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.; 
    d0_max = 5;

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

    mumueta_min = -10.;
    mumueta_max = 10.;

    mumuphi_min = -3.5;
    mumuphi_max = 3.5;

    mumupt_min = 0.;
    mumupt_max = 100.;
/*
    BDT_1_2_min = -0.25;
    BDT_1_2_max = 0.2;

    BDT_2_3_min = -0.85;
    BDT_2_3_max = 0.25;
*/
    BDT_3_5_min = -0.2;
    BDT_3_5_max = 0.85;

    BDT_5_7_min = -1;
    BDT_5_7_max = 1;

    BDT_7_10_min = -1;
    BDT_7_10_max = 1;

    BDT_10_15_min = -1;
    BDT_10_15_max = 1;

    BDT_15_20_min = -1;
    BDT_15_20_max = 1;

    BDT_20_50_min = -1;
    BDT_20_50_max = 1;

    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);  
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar nMult("nMult", "nMult", nMult_min, nMult_max);
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
    //RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
    RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
    RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    //RooRealVar BDT_pt_1_2("BDT_pt_1_2", "BDT_pt_1_2", BDT_1_2_min, BDT_1_2_max);
    //RooRealVar BDT_pt_2_3("BDT_pt_2_3", "BDT_pt_2_3", BDT_2_3_min, BDT_2_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_5_7", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);

    RooRealVar Bgen("Bgen", "Bgen", 0, 30000);

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
    //w.import(Btrk1DzError1);
    w.import(Btrk1Dxy1);
    w.import(Btrk1DxyError1);
    w.import(Bmumueta);
    w.import(Bmumuphi);
    w.import(Bmumupt);
    //w.import(BDT_pt_1_2);
    //w.import(BDT_pt_2_3);
    // w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    w.import(nMult);
    w.import(Bgen);
}
      
 else if(particle == 1){
    double mass_min, mass_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double nMult_min, nMult_max;
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
    //double trk2Dz1_min, trk2Dz1_max;
    double trk1DzError1_min, trk1DzError1_max;
    //double trk2DzError1_min, trk2DzError1_max;
    double trk1Dxy1_min, trk1Dxy1_max;
    //double trk2Dxy1_min, trk2Dxy1_max;
    double trk1DxyError1_min, trk1DxyError1_max;
    //double trk2DxyError1_min, trk2DxyError1_max;
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
    // double BDT_30_50_min, BDT_30_50_max;

    mass_min = 5.0;
    mass_max = 6.0;

    y_min = -2.8;
    y_max = 2.8;

    pt_min = 0;
    pt_max = 100.;  //100
      
    nMult_min = 0.;
    nMult_max = 130;   

    trk1pt_min = 0.;
    trk1pt_max = 10.;

    trk2pt_min = 0.;
    trk2pt_max = 10.;

    trk1eta_min = -2.6;
    trk1eta_max = 2.6;

    trk2eta_min = -2.6;
    trk2eta_max = 2.6;

    trk1pterr_min = 0.;
    trk1pterr_max = 0.3;

    trk2pterr_min = 0.;
    trk2pterr_max = 0.15;

    chi2cl_min = 0.0;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 20.;

    svpvDisErr_min = 0.;
    svpvDisErr_max = 0.1;

    svpvDistance2D_min = 0.;
    svpvDistance2D_max = 0.8;

    svpvDisErr2D_min = 0.;
    svpvDisErr2D_max = 0.05;

    mumumass_min = 2.9;
    mumumass_max = 3.3;

    mu1eta_min = -2.5;
    mu1eta_max = 2.5;

    mu2eta_min = -2.5;
    mu2eta_max = 2.5;

    mu1pt_min = 0.;
    mu1pt_max = 20.;

    mu2pt_min = 0.;
    mu2pt_max = 20.;

    mu1dxyPV_min = -0.15;
    mu1dxyPV_max = 0.15;

    mu2dxyPV_min = -0.15;
    mu2dxyPV_max = 0.15;

    mu1dzPV_min = -20.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 0.6;

    d0err_min = 0.;
    d0err_max = 0.0005;

    dtheta_min = 0.;
    dtheta_max = 3.4;

    alpha_min = 0.;
    alpha_max = 3.4;

    trk1Dz1_min = -20.;
    trk1Dz1_max = 20.;

    //trk2Dz1_min = -25.;
    //trk2Dz1_max = 20.;

    trk1DzError1_min = 0.0;
    trk1DzError1_max = 0.5;

    //trk2DzError1_min = -0.5;
    //trk2DzError1_max = 0.5;

    trk1Dxy1_min = -0.5;
    trk1Dxy1_max = 0.5;

    //trk2Dxy1_min = -1.;
    //trk2Dxy1_max = 1.;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.3;

    //trk2DxyError1_min = 0.;
    //trk2DxyError1_max = 0.4;

    mumueta_min = -4.5;
    mumueta_max = 4.5;

    mumuphi_min = -3.4;
    mumuphi_max = 3.4;

    mumupt_min = 0.;
    mumupt_max = 100.;

    BDT_1_2_min = -0.75;
    BDT_1_2_max = 0.05;

    BDT_2_3_min = -0.55;
    BDT_2_3_max = 0.35;

    BDT_3_5_min = -0.55;
    BDT_3_5_max = 0.1;

    BDT_5_7_min = -1;
    BDT_5_7_max = 1;

    BDT_7_10_min = -1;
    BDT_7_10_max = 1;

    BDT_10_15_min = -1;
    BDT_10_15_max = 1;

    BDT_15_20_min = -1;
    BDT_15_20_max = 1;

    BDT_20_50_min = -1;
    BDT_20_50_max = 1;

    // BDT_30_50_min = -0.7;
    // BDT_30_50_max = 0.85;

    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar nMult("nMult", "nMult", nMult_min, nMult_max);
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
    //RooRealVar Btrk2Dz1("Btrk2Dz1","Btrk2Dz1",trk2Dz1_min,trk2Dz1_max);
    //RooRealVar Btrk2DzError1("Btrk2DzError1","Btrk2DzError1",trk2DzError1_min,trk2DzError1_max);
    //RooRealVar Btrk2Dxy1("Btrk2Dxy1","Btrk2Dxy1",trk2Dxy1_min,trk2Dxy1_max);
    //RooRealVar Btrk2DxyError1("Btrk2DxyError1","Btrk2DxyError1",trk2DxyError1_min,trk2DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    RooRealVar BDT_pt_1_2("BDT_pt_1_2", "BDT_pt_1_2", BDT_1_2_min, BDT_1_2_max);
    RooRealVar BDT_pt_2_3("BDT_pt_2_3", "BDT_pt_2_3", BDT_2_3_min, BDT_2_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_5_7", BDT_5_7_min, BDT_5_7_max);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", BDT_7_10_min, BDT_7_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);
    // RooRealVar BDT_pt_30_50("BDT_pt_30_50", "BDT_pt_30_50", BDT_30_50_min, BDT_30_50_max);

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
    //w.import(Btrk2Dz1);
    //w.import(Btrk2DzError1);
    //w.import(Btrk2Dxy1);
    //w.import(Btrk2DxyError1);
    w.import(BDT_pt_1_2);
    w.import(BDT_pt_2_3);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    // w.import(BDT_pt_30_50);
    w.import(nMult);
  }

  else if(particle == 2){
    double mass_min, mass_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double nMult_min, nMult_max;
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
    //double trk2Dz1_min, trk2Dz1_max;
    double trk1DzError1_min, trk1DzError1_max;
    //double trk2DzError1_min, trk2DzError1_max;
    double trk1Dxy1_min, trk1Dxy1_max;
    //double trk2Dxy1_min, trk2Dxy1_max;
    double trk1DxyError1_min, trk1DxyError1_max;
    //double trk2DxyError1_min, trk2DxyError1_max;
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

    mass_min = 5;
    mass_max = 5.6;

    y_min = -2.8;
    y_max = 2.8;

    pt_min = 0.;
    pt_max = 100.;
    
    nMult_min = 0;
    nMult_max = 130;
    
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

    //trk2Dz1_min = -22.;
    //trk2Dz1_max = 20.;

    trk1DzError1_min = -0.5;
    trk1DzError1_max = 0.5;

    //trk2DzError1_min = 0.;
    //trk2DzError1_max = 0.5;

    trk1Dxy1_min = -20;
    trk1Dxy1_max = 20;

    //trk2Dxy1_min = -0.25;
    //trk2Dxy1_max = 0.25;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.2;

    //trk2DxyError1_min = 0.;
    //trk2DxyError1_max = 0.1;

    mumueta_min = -4.;
    mumueta_max = 4.;

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

    BDT_5_7_min = -0.5;
    BDT_5_7_max = 0.8;

    BDT_7_10_min = -0.4;
    BDT_7_10_max = 0.8;

    BDT_10_15_min = -0.45;
    BDT_10_15_max = 0.8;

    BDT_15_20_min = -0.6;
    BDT_15_20_max = 0.8;

    BDT_20_50_min = -0.7;
    BDT_20_50_max = 0.9;

    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar nMult("nMult", "nMult", nMult_min, nMult_max);
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
    RooRealVar Btrk1Dz1("Btrk1Dz1", "Btrk1Dz1", trk1Dz1_min, trk1Dz1_max);
    RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min, trk1DzError1_max);
    RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyError1_min,trk1DxyError1_max);
    RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1", trk1Dxy1_min, trk1Dxy1_max);
    //RooRealVar Btrk2Dz1("Btrk2Dz1","Btrk2Dz1",trk2Dz1_min,trk2Dz1_max);
    //RooRealVar Btrk2DzError1("Btrk2DzError1","Btrk2DzError1",trk2DzError1_min,trk2DzError1_max);
    //RooRealVar Btrk2Dxy1("Btrk2Dxy1","Btrk2Dxy1",trk2Dxy1_min,trk2Dxy1_max);
    //RooRealVar Btrk2DxyError1("Btrk2DxyError1","Btrk2DxyError1",trk2DxyError1_min,trk2DxyError1_max);
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    RooRealVar BDT_pt_0_2("BDT_pt_0_2", "BDT_pt_0_2", BDT_0_2_min, BDT_0_2_max);
    RooRealVar BDT_pt_2_3("BDT_pt_2_3", "BDT_pt_2_3", BDT_2_3_min, BDT_2_3_max);
    RooRealVar BDT_pt_3_5("BDT_pt_3_5", "BDT_pt_3_5", BDT_3_5_min, BDT_3_5_max);
    RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_5_7", BDT_5_7_min, BDT_5_7_max);
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
    //w.import(Btrk2Dz1);
    //w.import(Btrk2DzError1);
    //w.import(Btrk2Dxy1);
    //w.import(Btrk2DxyError1);
    w.import(BDT_pt_0_2);
    w.import(BDT_pt_2_3);
    w.import(BDT_pt_3_5);
    w.import(BDT_pt_5_7);
    w.import(BDT_pt_7_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    w.import(nMult);
  }
}

// get ystring
TString ystring(int iy) {
  TString str = TString::Format("y%.0f-%.0f", 10 * ylist[iy], 10 * ylist[iy + 1]);
  return str;
}

// save the plots MC/SS, MC/
void save_validation_plot(TCanvas& can, TString name, TString comp, TString ptdir, int iy) {
   TString pdfstr;

   TString ystr = "";
   pdfstr.Form("%s/mc_validation_plots/%s/%s_mc_validation_%s%s.pdf",
               ptdir.Data(), comp.Data(), name.Data(),
               particleList.at(particle).Data(), ystr.Data());

    gSystem->mkdir(Form("%s/mc_validation_plots/", ptdir.Data()) ); //create next folder
    gSystem->mkdir(Form("%s/mc_validation_plots/%s", ptdir.Data(), comp.Data()) ); //create next folder

   can.SaveAs(pdfstr);
 }

/** Fix or release the parameters for signal shape */
void fix_signal_shape(RooWorkspace& w, bool release) {
  fix_parameters(w, "signal", release);
  // never fix the overall width
  RooRealVar* sigma1 = w.var("sigma1");
  sigma1->setConstant(false);
}

/** Fix or release the parameters for a given PDF */
void fix_parameters(RooWorkspace& w, TString pdfName, bool release) {
  RooAbsPdf* pdf = w.pdf(pdfName);
  RooAbsData* ds = w.data("data");
  RooArgSet* par_set = pdf->getParameters(*ds);
  auto itr = par_set->createIterator();
  bool toFix = ! release;
  std::string fix_or_float = (toFix)? "fix " : "float ";
  std::cout << fix_or_float << "parameters:";
  for (auto i = 0; i < par_set->getSize(); ++i) {
    RooRealVar* var = (RooRealVar*) itr->Next();
    var->setConstant(true);
    TString name = var->GetName();
    std::cout << name << ", ";
    var->setConstant(toFix);
  }
  std::cout << "\n";
}

/** Fix or release the parameters for a given arg list */
void fix_parameters(RooWorkspace& w, RooArgList& parlist, bool release) {
  for (auto i = 0; i < parlist.getSize(); ++i) {
    TString name = parlist[i].GetName();
    RooRealVar* var = w.var(name);
    bool toFix = ! release;
    std::string fix_or_float = (toFix)? "fix " : "float ";
    std::cout << fix_or_float << "parameters:";
    var->setConstant(toFix);
  }
}

/** plot the mc fit */
template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds,
                TString plotName, TString title, Targs... options) {
  TString label = "";
  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(title));
  ds->plotOn(massframe, MarkerSize(0.9));
  model->plotOn(massframe, options...);

  model->paramOn(massframe, Layout(0.5, 0.8, 0.9),
                 Label(label),
                 Format("NEU", AutoPrecision(1)) ) ;

  TCanvas can_mc;
  can_mc.SetTitle("");
  massframe->Draw();
  can_mc.SaveAs(plotName);
}

void plot_jpsifit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds,
                  TString plotName, TString title, double nGen, RooRealVar& n_signal) {
  bool with_sig = (nGen > 1);
  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(title));
  ds->plotOn(massframe, RooFit::Name("NP"), MarkerSize(0.9));
  model->plotOn(massframe, RooFit::Name("NP Fit"), NormRange("bmass"),
                LineColor(kRed), LineStyle(1), LineWidth(2));
  if (with_sig) {
    model->plotOn(massframe, RooFit::Name("signal"),
                  Components("signal"), NormRange("bmass"),
                  LineColor(kOrange-3), LineStyle(1), LineWidth(3), FillStyle(3002),
                  FillColor(kOrange-3), VLines(), DrawOption("LF"));
    model->plotOn(massframe, RooFit::Name("B->J/#psi #pi"),
                  Components("jpsipi"), NormRange("bmass"),
                  LineColor(kPink+10), LineStyle(kDashed));
  }
  model->plotOn(massframe, RooFit::Name("peaking"),
                Components("erf"), NormRange("bmass"),
                LineColor(kGreen+3), LineStyle(1), LineWidth(3), FillStyle(3005),
                FillColor(kGreen+3), VLines(), DrawOption("LF"));
  model->plotOn(massframe, RooFit::Name("poly"),
                Components("poly_jpsi"), NormRange("bmass"),
                LineColor(kBlue), LineStyle(kDashed));
  model->paramOn(massframe,  Layout(0.35, 0.6, 0.9),
                 Format("NEU", AutoPrecision(1)));

  TCanvas can_np;
  can_np.SetTitle("");
  massframe->Draw();

  TLatex txt;
  TLegend *leg = new TLegend (0.7, 0.35, 0.85, 0.60);
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(massframe->findObject("NP"), "NP MC", "p)");
  leg->AddEntry(massframe->findObject("peaking"), "bg as erf", "f");
  leg->AddEntry(massframe->findObject("poly"), "combinatorial", "l");
  leg->AddEntry(massframe->findObject("NP Fit"),"Fit","l");
  // compare yields with gen particles
  if (with_sig) {
    leg->AddEntry(massframe->findObject("signal"), "signal", "f");
    leg->AddEntry(massframe->findObject("B->J/#psi #pi"), "B->J/#psi #pi", "f");
    txt.DrawLatexNDC(0.42, 0.4, TString::Format("gen: %.0f", nGen));
    txt.DrawLatexNDC(0.42, 0.32, TString::Format("fit: %.0f #pm %.0f",
                                                n_signal.getVal(),
                                                n_signal.getError()));
  }
  leg->Draw();
  can_np.SaveAs(plotName);
}
