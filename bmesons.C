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

std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, int* n, int n_var);
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


// DATA_CUT
// 1 = apply cuts, recd ..strict variable range when reading data -- to be used for mc validation
// 0 = read full data
// note: when reading tratio should assign weight=1 for events out of range

#define DATA_CUT 1

//particle
// 0 = Bu
// 1 = Bs

#define particle 1

//weights
// 1 = calculates ratio between MC and sPlot 
// 0 = does not calculate weights 

#define weights 0

//add_weights
// 1 = adds weights to tree
// 0 = does not add weights to tree

# define add_weights 0


void bmesons(){
  
  int n_var;
  
  TString input_file_data = particle ? "/lstore/cms/nuno/ppdata2017/000814/BsData.root" : "/lstore/cms/nuno/ppdata2017/000814/BPData.root";
  TString input_file_mc   = particle ? "/lstore/cms/nuno/ppdata2017/000814/BsMC.root" : "/lstore/cms/nuno/ppdata2017/000814/BPMC.root";
  TString input_file_reweighted_mc   = particle ? "./results/Bs/mc_validation_plots/weights/tree_with_weight.root" : "./results/Bu/mc_validation_plots/weights/tree_with_weight.root";
  
  std::vector<TH1D*> histos_sideband_sub;
  std::vector<TH1D*> histos_mc;
  std::vector<TH1D*> histos_splot;

#if particle == 0
    int n_bins[] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20} ;

    TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt"} ;    

#elif particle == 1
  int n_bins[] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20} ;

  TString variables[] = {"By", "Bpt", "Btrk1Pt", "Btrk1Eta", "Btrk1PtErr", "Bchi2cl", "BsvpvDistance", "BsvpvDisErr", "Bmumumass", "Bmu1eta", "Bmu2eta", "Bmu1pt", "Bmu2pt", "Bmu1dxyPV", "Bmu2dxyPV", "Bmu1dzPV", "Bmu2dzPV", "Bd0", "Bd0Err", "Bdtheta", "Balpha", "Btrk1Dz1", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bmumueta", "Bmumuphi", "Bmumupt", "BDT_pt_5_10","BDT_pt_10_15", "BDT_pt_15_20", "BDT_pt_20_50"} ;

#endif

  int n_n_bins = sizeof(n_bins)/sizeof(n_bins[0]);
  //std::cout << n_n_bins << std::endl;
  int n_variables = sizeof(variables)/sizeof(variables[0]);
  //std::cout << n_variables << std::endl;

  if(n_n_bins != n_variables){
    std::cout << "Error: number of bins does not correspond to number of variables." << std::endl;
    return;
  }

  n_var = n_variables;
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,input_file_data);
	
  build_pdf(*ws);

  if (DATA_CUT == 0){
    plot_complete_fit(*ws);
    validate_fit(ws);
  }
  
  //if(!DATA_CUT){fit_syst_error(input_file_data);}

  //sideband_sub histograms
  histos_sideband_sub = sideband_subtraction(ws, n_bins, n_var);

  do_splot(*ws);
  histos_splot = splot_method(*ws,n_bins,variables, n_var);

  //monte carlo histograms
  //TFile *fin_mc = new TFile(input_file_mc); //use this file to add the weights (to clone original tree) and make data-MC comparisons without weights
  TFile *fin_mc = new TFile(input_file_reweighted_mc); //use this file to make data-MC comparisons with weights

  TTree* t1_mc = particle ? (TTree*)fin_mc->Get("ntphi") : (TTree*)fin_mc->Get("ntKp");

  std::vector<TString> names;

  for(int i=0; i<n_var; i++){
    //std::cout<< "Var names: "<< histos_sideband_sub[i]->GetName()<<std::endl;
    TString weight = "pthat * weight * weight_Balpha";
    histos_mc.push_back(create_histogram_mc((*ws->var(variables[i])), t1_mc, n_bins[i], weight));
    names.push_back(TString(variables[i]));
  }

  //get the ratio between the data (splot method) and the MC 
 
  if (weights == 1){
    get_ratio(histos_splot, histos_mc,names,"weights.root");
  }
    
  if (add_weights == 1){
    AddWeights(t1_mc);    
  }  
  
  //if(!DATA_CUT){pT_analysis(*ws,n_bins[0], "pT.root", input_file_data);}
 
  //COMPARISONS//
  
  //Sideband Subtraction vs. Monte Carlo
  //clone
  vector<TH1D*> mc_comp_ss(histos_mc);

  vector<TH1D*> ss_comp_mc(histos_sideband_sub);

  if(DATA_CUT == 1){
    for(int i=0; i<n_var; i++) {
      TCanvas c;
 
      mc_comp_ss[i]->SetXTitle(TString(ss_comp_mc[i]->GetName()));

      mc_comp_ss[i]->SetStats(0);
      ss_comp_mc[i]->SetStats(0);

      if(particle == 0){
          mc_comp_ss[i]->SetTitle("B^{+}");
      }
      else if (particle == 1){
          mc_comp_ss[i]->SetTitle("B^{0}_{s}");
      }

      //normalization
      mc_comp_ss[i]->Scale(1/mc_comp_ss[i]->Integral());
      ss_comp_mc[i]->Scale(1/ss_comp_mc[i]->Integral());

      mc_comp_ss[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_ss[i]->GetMinimum(),1.1*mc_comp_ss[i]->GetMaximum());
      mc_comp_ss[i]->Draw();
      ss_comp_mc[i]->Draw("same");

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
	c.SaveAs("./results/Bu/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bu.gif");
      }else if(particle == 1){
	c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.pdf");
	c.SaveAs("./results/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validation_Bs.gif");
      }
      //leg->Delete();
      
    }
  }

  //SPlot vs. Sideband subtraction

  //clone
  vector<TH1D*> sp_comp_ss(histos_splot);
  vector<TH1D*> ss_comp_sp(histos_sideband_sub);

  if(DATA_CUT == 1){
    for(int i=0; i<n_var; i++)
      {
	TCanvas a;
	ss_comp_sp[i]->SetYTitle("normalized entries");
	sp_comp_ss[i]->SetXTitle(TString(ss_comp_sp[i]->GetName()));
	ss_comp_sp[i]->SetStats(0);
	sp_comp_ss[i]->SetStats(0);

        /*if(particle == 0){
          ss_comp_sp[i]->SetTitle("B^{+}");
        }
        else if (particle == 1){
          ss_comp_sp[i]->SetTitle("B^{0}_{s}");
        }*/
	
	//normalization
	ss_comp_sp[i]->Scale(1/ss_comp_sp[i]->Integral());
	sp_comp_ss[i]->Scale(1/sp_comp_ss[i]->Integral());
	
	
	ss_comp_sp[i]->GetYaxis()->SetRangeUser(0.1*ss_comp_sp[i]->GetMinimum(),1.1*ss_comp_sp[i]->GetMaximum());
	ss_comp_sp[i]->Draw();
	sp_comp_ss[i]->Draw("same");
	
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
	  a.SaveAs("./results/Bu/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bu.gif");
	}else if(particle == 1){
	  a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.pdf");
	  a.SaveAs("./results/Bs/mc_validation_plots/ss_sp/" + names[i]+"_mc_validation_Bs.gif");
	}
	//leg->Delete();
      }
  }     
  
  //SPlot vs. Monte Carlo

  //clone
  vector<TH1D*> sp_comp_mc(histos_splot);
  vector<TH1D*> mc_comp_sp(histos_mc);

  if(DATA_CUT == 1){
    for(int i=0; i<n_var; i++)
      {
	TCanvas a;
	mc_comp_sp[i]->SetXTitle(TString(histos_sideband_sub[i]->GetName()));
	mc_comp_sp[i]->SetYTitle("normalized entries");
	sp_comp_mc[i]->SetXTitle(TString(histos_sideband_sub[i]->GetName()));
	
	mc_comp_sp[i]->SetStats(0);
	sp_comp_mc[i]->SetStats(0);

        if(particle == 0){
          mc_comp_sp[i]->SetTitle("B^{+}");
        }
        else if (particle == 1){
          mc_comp_sp[i]->SetTitle("B^{0}_{s}");
        }
	
	//normalization
	mc_comp_sp[i]->Scale(1/mc_comp_sp[i]->Integral());
	sp_comp_mc[i]->Scale(1/sp_comp_mc[i]->Integral());
	
	mc_comp_sp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp_sp[i]->GetMinimum(),1.1*mc_comp_sp[i]->GetMaximum());
	mc_comp_sp[i]->Draw();
	sp_comp_mc[i]->Draw("same");
	
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
	  a.SaveAs("./results/Bu/mc_validation_plots/mc_sp/" + names[i]+"_mc_validation_Bu.gif");
	}else if(particle == 1){
	  a.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
	  a.SaveAs("./results/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validation_Bs.gif");
	}
	
	//leg->Delete();
      }
  }
  //Sideband subtraction vs. Monte Carlo vs SPlot
  
  //clone
  vector<TH1D*> sp_comp(histos_splot);
  vector<TH1D*> mc_comp(histos_mc);
  vector<TH1D*> ss_comp(histos_sideband_sub);

  if(DATA_CUT == 1){
    for(int i=0; i<n_var; i++)
      {
	TCanvas a;
	
	mc_comp[i]->SetXTitle(TString(ss_comp[i]->GetName()));
	mc_comp[i]->SetYTitle("normalized entries");
	sp_comp[i]->SetXTitle(TString(ss_comp[i]->GetName()));
	mc_comp[i]->SetStats(0);
	sp_comp[i]->SetStats(0);
	ss_comp[i]->SetStats(0);

        if(particle == 0){
          mc_comp[i]->SetTitle("B^{+}");
        }
        else if (particle == 1){
          mc_comp[i]->SetTitle("B^{0}_{s}");
        }
	
	//normalization
	mc_comp[i]->Scale(1/mc_comp[i]->Integral());
	sp_comp[i]->Scale(1/sp_comp[i]->Integral());
	ss_comp[i]->Scale(1/ss_comp[i]->Integral());
	
	
	//y axis: maximum and minimum 
	if ( ( mc_comp[i]->GetMaximum() > sp_comp[i]->GetMaximum() ) && ( mc_comp[i]->GetMaximum() > ss_comp[i]->GetMaximum() ) ){
	  mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*mc_comp[i]->GetMaximum());
	}
	else if ( (sp_comp[i]->GetMaximum() > ss_comp[i]->GetMaximum() ) && ( sp_comp[i]->GetMaximum() > mc_comp[i]->GetMaximum() ) ){
	  mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*sp_comp[i]->GetMaximum());
	}
	else {
	  mc_comp[i]->GetYaxis()->SetRangeUser(0.1*mc_comp[i]->GetMinimum(), 1.1*ss_comp[i]->GetMaximum());
	}
	
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
	  a.SaveAs("./results/Bu/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bu.gif");
	}else if(particle == 1){
	  a.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.pdf");
	  a.SaveAs("./results/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validation_Bs.gif");
	}
	
	//leg->Delete();
      }
  }
  
  //comparisons end

}

//main function ends

void pT_analysis(RooWorkspace& w, int n, TString ptfile, TString datafile){

  TString dir_name = particle ? "./results/Bs/Bpt/" : "./results/Bu/Bpt/";

  TFile* f_wei = new TFile(dir_name + "/"+ ptfile, "recreate"); 

  RooAbsPdf*  model = w.pdf("model");
  RooRealVar* Bpt  = w.var("Bpt");
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooRealVar Bmass = *(w.var("Bmass"));

#if particle == 0
  const int n_pt_bins = 7;
  double pt_bins [n_pt_bins + 1] = {5,7,10,15,20,30,50,100};  
#elif particle == 1
  const int n_pt_bins = 4;
  double pt_bins[n_pt_bins + 1] = {5,10,15,20,50};
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

  
  //plots the signal+background and signal distributions in linear and log scales
  TCanvas* a = new TCanvas("pT","pT", 800, 600);
  a->Divide(2,2);

  //signal+bkg distribution

  //linear scale
  a->cd(1);
  RooPlot* ptframe = Bpt->frame();
  data->plotOn(ptframe);
  if(particle == 0){
    ptframe->SetTitle("pT of B+: total sample");
  }else if(particle == 1){
    ptframe->SetTitle("pT of Bs: total sample");
  }
  ptframe->Draw();

  //log scale
  a->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  data->plotOn(ptframe);
  
  if(particle == 0){
    ptframe->SetTitle("pT of B+: total sample");
  }else if(particle == 1){
    ptframe->SetTitle("pT of Bs: total sample");
  }
  ptframe->SetMinimum(1);
  ptframe->Draw();
  
  //signal distribution
  RooDataSet* dataWBp = (RooDataSet*) w.data("dataWBp");

  //linear scale
  a->cd(3);
  RooPlot* ptframe2Bp = Bpt->frame();
  ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(Bpt->getMax()-Bpt->getMin())/n));
  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));

  if(particle == 0){
    ptframe2Bp->SetTitle("Bpt distribution of B+ for signal (splot)");
    ptframe2Bp->GetXaxis()->SetTitle("Bpt of B+");
  }else if(particle == 1){
    ptframe2Bp->SetTitle("Bpt distribution of Bs for signal (splot)");
    ptframe2Bp->GetXaxis()->SetTitle("Bpt of Bs");
  }

  ptframe2Bp->Draw();

  //log scale
  a->cd(4);
  gPad->SetLogx();
  gPad->SetLogy();
  //ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(Bpt->getMax()-Bpt->getMin())/n));
  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));

  if(particle == 0){
    ptframe2Bp->SetTitle("Bpt distribution of B+ for signal (splot)");
    ptframe2Bp->GetXaxis()->SetTitle("Bpt of B+");
  }else if(particle == 1){
    ptframe2Bp->SetTitle("Bpt distribution of Bs for signal (splot)");
    ptframe2Bp->GetXaxis()->SetTitle("Bpt of Bs");
  }

  ptframe2Bp->SetMinimum(1);
  ptframe2Bp->Draw();

  if(particle == 0){
    a->SaveAs("./results/Bu/Bpt/pTdistributions_Bu.gif");
    a->SaveAs("./results/Bu/Bpt/pTdistributions_Bu.pdf");
  }else if(particle == 1){
    a->SaveAs("./results/Bs/Bpt/pTdistributions_Bs.gif");
    a->SaveAs("./results/Bs/Bpt/pTdistributions_Bs.pdf");
  }
  
 
  //applies the splot method and evaluates the weighted average pT per bin

  const int n_pdf_syst=4;
  //number of pdfs
  TString syst_src[n_pdf_syst]={"nominal","bkg_poly","signal1gauss","bkg_range"};
  //array of pdfs to evaluate the systematic uncertainty of the N signal
  double yield_syst[n_pt_bins][n_pdf_syst]; 
  //value of the systematic uncertainty

  
  for(int i=0;i<n_pt_bins;i++){
    //select data subset corresponding to pT bin
    data_pt = (RooDataSet*) data->reduce(Form("Bpt>%lf",pt_bins[i]));
    data_pt = (RooDataSet*) data_pt->reduce(Form("Bpt<%lf",pt_bins[i+1]));
    w.import(*data_pt, Rename(Form("data_pt_%d",i)));
   
    //perform fit and save result
    fit_pt = model->fitTo(*data_pt, Minos(true), Save());

    //plots the fit result
    RooPlot* massframe = Bmass.frame(Title(""));
    data_pt->plotOn(massframe);
    model->paramOn(massframe,Layout(0.60,0.90,0.75));
    model->plotOn(massframe, Range("all"));

    TCanvas b;
    massframe->Draw();
    if(particle == 0){
      b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf_fit_plot_Bu.gif", pt_bins[i], pt_bins[i+1]));
      b.SaveAs(Form("./results/Bu/Bpt/%lf_%lf_fit_plot_Bu.pdf", pt_bins[i], pt_bins[i+1]));
    }else if(particle == 1){
      b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf_fit_plot_Bs.gif", pt_bins[i], pt_bins[i+1]));
      b.SaveAs(Form("./results/Bs/Bpt/%lf_%lf_fit_plot_Bs.pdf", pt_bins[i], pt_bins[i+1]));
    }

    //get yield and its errors

    //floatParsFinal returns the list of floating parameters after fit
    //cout << "Value of floating parameters" << endl;
    //fit_pt->floatParsFinal().Print("s");

    //signal yield
    n_sig_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_signal");
    //combinatorial background yield
    n_comb_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_combinatorial");

    yield[i]= n_sig_pt->getVal();
    yield_err_low [i] = n_sig_pt->getError(); 
    yield_err_high[i] = n_sig_pt->getError(); 

    cout << "test asym error:" << n_sig_pt->getErrorLo() << " " <<  n_sig_pt->getAsymErrorLo() << " symmetric: " <<  n_sig_pt->getError() <<  endl;

    //Loop through the array of pdfs 
    //for(int k=0;k<n_pdf_syst;k++){
    //dúvida
    for(int k = 1; k<n_pdf_syst; k++){
      double val = 0.; 
      double val_nominal = 0.;
      val = get_yield_syst(data_pt, syst_src[k]);
      val_nominal = get_yield_syst(data_pt, syst_src[0]);
      cout<<"syst nominal: "<<syst_src[0]<<endl;
      //yield_syst_rel[i][k] = (val - val_nominal)/val_nominal;
      yield_syst[i][k] = (val - val_nominal);
      cout << "bin:" << i << " range bin min" << pt_bins[i] << " src:" << k << " syst:" << syst_src[k] << " yield syst:" << yield_syst[i][k]<< " rel:"<<  ((val - val_nominal)/val_nominal)*100 << "\%"<<endl;
    }  

    for(int k = 1; k<n_pdf_syst; k++){ 
      yield_err_syst[i] += pow(yield_syst[i][k],2);
    }

    m_yield_err_syst[i] = sqrt(yield_err_syst[i]);
 
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
    } else {
      pt_low [i]= 0.5*(pt_bins[i+1]-pt_bins[i]);
      pt_high[i]= 0.5*(pt_bins[i+1]-pt_bins[i]);
  
    }

    //normalize yield to bin width
    double bin_width = pt_bins[i+1]-pt_bins[i];
    yield[i] = yield[i]/bin_width;
    yield_err_low[i] = yield_err_low[i]/bin_width;
    yield_err_high[i] = yield_err_high[i]/bin_width;
    m_yield_err_syst[i] = m_yield_err_syst[i]/bin_width;

    //cout<<"pt: "<< pt_bins[i]<<"-" << pt_bins[i+1] << " mean_weight:"<<mean_w<< "  Nsig:" <<yield[i]<< "-"<<yield_err_low[i] <<"+" << yield_err_high[i]<<endl;
    // cout <<"pt: "<< pt_bins[i]<<"-" << pt_bins[i+1] <<" systematic uncertainty" <<m_yield_err_syst[i]<<endl;

  }

  //plot yield vs average pT

  TCanvas c;
  TMultiGraph* mg = new TMultiGraph();

  TGraphAsymmErrors* gr = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
  gr->SetTitle("");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(1);
  gr->SetLineColor(1);
  gr->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  gr->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
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
  grs->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  grs->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
  grs->Write();
  f_wei->Close();
  //grs->Draw("same");

  
   mg->Add(gr);
   mg->Add(grs);
   mg->Draw("AP");
   mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
   mg->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
   
  
 
  if(particle == 0){
    c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.pdf");
    c.SaveAs("./results/Bu/Bpt/raw_yield_pt_Bu.gif");}
  else if(particle == 1){
    c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.pdf");
    c.SaveAs("./results/Bs/Bpt/raw_yield_pt_Bs.gif");}

//save both vectors (systematic and statistical errors) in a text file


  TCanvas l;
  //log scale
  l.SetLogx();
  l.SetLogy();
  TGraphAsymmErrors* grlog = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
  grlog->SetTitle("");
  grlog->SetMarkerColor(4);
  grlog->SetMarkerStyle(21);
  grlog->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  grlog->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
  grlog->Draw("AP");

  if(particle == 0){
    l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.pdf");
    l.SaveAs("./results/Bu/Bpt/raw_yield_pt_logscale_Bu.gif");}
  else if(particle == 1){
    l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.pdf");
    l.SaveAs("./results/Bs/Bpt/raw_yield_pt_logscale_Bs.gif");}


  //evaluates the N systematics for each bin
  //for(int i=0;i<n_pt_bins;i++){
  // fit_syst_error_bin(datafile, pt_bins[i], pt_bins[i+1]);
  // }


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
  //TNtupleD* _nt = (TNtupleD*)fin_data->Get("ntKp");
  TTree* t1_data = particle ? (TTree*)fin_data->Get("ntphi") : (TTree*)fin_data->Get("ntKp"); 
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
  
  
  /*if(particle == 0){
    arg_list.add(*(w.var("Btrk1Dz1")));
    arg_list.add(*(w.var("Btrk1DzError1")));
    arg_list.add(*(w.var("Btrk1Dxy1")));
    arg_list.add(*(w.var("Btrk1DxyError1")));
    arg_list.add(*(w.var("Bd0")));
    arg_list.add(*(w.var("Bd0Err")));
    arg_list.add(*(w.var("BDT_pt_5_7")));
    arg_list.add(*(w.var("BDT_pt_7_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_30")));
    arg_list.add(*(w.var("BDT_pt_30_50")));
    arg_list.add(*(w.var("BDT_pt_50_100")));
    arg_list.add(*(w.var("BDT_total")));
  }
  */
  
  if(particle == 1){
    arg_list.add(*(w.var("BDT_pt_5_10")));
    arg_list.add(*(w.var("BDT_pt_10_15")));
    arg_list.add(*(w.var("BDT_pt_15_20")));
    arg_list.add(*(w.var("BDT_pt_20_50")));
    //arg_list.add(*(w.var("BDT_total")));	 
  }
  
 
  RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);

  w.import(*data, Rename("data"));

}

void build_pdf(RooWorkspace& w, std::string choice) {

  cout << "choice " << choice << endl;

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*) w.data("data");

  cout << "buid pdf 0 " << choice << endl;

  //  RooDataSet* reduceddata_central;

  double left = particle ? 5.3 : 5.15; 
  double right = particle ? 5.45 : 5.4; 
  double mass_peak = particle ? 5.366 : 5.265; 

  cout << "buid pdf 1 " << choice << endl;
  
  //  reduceddata_central = (RooDataSet*) data->reduce(Form("Bmass>%lf",left));

  cout << "buid pdf 2 " << choice << endl;

  //  reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Bmass<%lf",right));

  cout << "buid pdf 3 " << choice << endl;

  //SIGNAL//

  RooRealVar mean("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.1);
  RooRealVar sigma1("sigma1","sigma1",0.0252,0.005,0.030);//0.03
  RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);
  RooRealVar sigma2("sigma2","sigma2",0.01052,0.005,0.020);
  RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);
  RooRealVar cofs("cofs", "cofs", 0.317, 0., 1.);
  RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);
  //sigma1.setConstant();
  //sigma2.setConstant();
  //cofs.setConstant();

  //systematics->signal
  RooRealVar sigma3("sigma3","sigma3",0.012,0.010,0.030);
  RooGaussian signal3("signal3","signal3",Bmass, mean, sigma3);
  
  
  
  //BACKGROUND//

  //error function
  RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
  RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  
  m_nonprompt_shift.setConstant(kTRUE);
  m_nonprompt_scale.setConstant(kTRUE);

  RooGenericPdf erf("erf","erf","TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Bmass,m_nonprompt_scale,m_nonprompt_shift));
 
  //exponential
  RooRealVar lambda("lambda","lambda",-2.,-5.,1.0);
  RooExponential fit_side("fit_side", "fit_side_exp", Bmass, lambda);

  //jpsi_pi component
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_sigma1l("m_jpsipi_sigma1l","m_jpsipi_sigma1l",2.90762e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma1r("m_jpsipi_sigma1r","m_jpsipi_sigma1r",6.52519e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma2("m_jpsipi_sigma2","m_jpsipi_sigma2",9.94712e-02,0.020,0.500);

  RooRealVar m_jpsipi_sigma3("m_jpsipi_sigma3","m_jpsipi_sigma3",3.30152e-01,0.020,0.500);
  RooRealVar m_jpsipi_fraction2("m_jpsipi_fraction2","m_jpsipi_fraction2",2.34646e-01,0.0,1.0);
  RooRealVar m_jpsipi_fraction3("m_jpsipi_fraction3","m_jpsipi_fraction3",1.14338e-01,0.0,1.0);

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



  //n values
  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());

  RooRealVar f_erf("f_erf","f_erf",0.2,0.0,1.);//2.50259e-01,0.0,1.
  RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));
  
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); 
  f_jpsipi.setConstant(kTRUE);
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));

  //systematics->bkg
  RooRealVar slope("slope","slope",0,-10,10);
  RooPolynomial poly_bkg("poly_bkg", "poly_bkg", Bmass, slope);

  if(particle == 0){ //B+
    if (choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }else if (choice=="bkg_poly"){
      RooAddPdf model("model", "model", RooArgList(signal,poly_bkg,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    } else if (choice=="bkg_range"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side, jpsipi),RooArgList(n_signal,n_combinatorial, n_jpsipi));
      w.import(model);
    }else if (choice == "signal1gauss"){
      RooAddPdf model("model", "model", RooArgList(signal3,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }
  }

  else if(particle == 1){
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

  build_pdf(*ws);
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

  build_pdf(*ws);
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
    model->paramOn(massframe,Layout(0.,0.,0.));
    massframe->getAttText()->SetTextSize(0.028);
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("Bmass (GeV)");
  }else if(particle == 1){
    data->plotOn(massframe, RooFit::Name("Data"));
    model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
    model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
    model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),Range("all"),LineColor(kOrange),LineStyle(kDashed));
    //model->paramOn(massframe,Layout(0.60,0.90,0.75));
    model->paramOn(massframe,Layout(0.,0.,0.));
    massframe->getAttText()->SetTextSize(0.028);
    massframe->GetYaxis()->SetTitleOffset(1.3);
    massframe->SetXTitle("Bmass (GeV)");
  }

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
  //p2->SetFrameBorderMode(0);
  //p2->SetBorderSize(1); 
  
  p2->Draw();

  p1->cd();

  massframe->Draw();
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

  TLegend *leg = new TLegend (0.7, 0.5, 0.9, 0.7);

  if(particle == 0){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/psi X", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/psi pi", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
    leg->Draw("same");
  }else if(particle == 1){
    leg->SetTextSize(0.03);
    leg->AddEntry(massframe->findObject("Data"), "Data", "l");
    leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
    leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
    leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
    leg->Draw("same");
  }
  
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
 
  //only saves the fit when no cuts are applied
  if(DATA_CUT == 0){
    if(particle == 0){
      d.SaveAs("./results/Bu/complete_fit_Bu.pdf");
      d.SaveAs("./results/Bu/complete_fit_Bu.gif");
    }else if(particle == 1){
      d.SaveAs("./results/Bs/complete_fit_Bs.pdf");
      d.SaveAs("./results/Bs/complete_fit_Bs.gif");
    }
  }
}

//plot_complete_fit ends

//SIDEBAND SUBTRACTION//
std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, int* n, int n_var){
  
  RooDataSet* data = (RooDataSet*) w->data("data");

  RooAbsPdf* fit_side = w->pdf("fit_side");

  vector<RooRealVar> variables;
  
  variables.push_back(*(w->var("Bmass")));
  variables.push_back(*(w->var("By")));
  variables.push_back(*(w->var("Bpt")));
  variables.push_back(*(w->var("Btrk1Pt")));
  variables.push_back(*(w->var("Btrk1Eta")));
  variables.push_back(*(w->var("Btrk1PtErr")));
  variables.push_back(*(w->var("Bchi2cl")));
  variables.push_back(*(w->var("BsvpvDistance")));
  variables.push_back(*(w->var("BsvpvDisErr")));
  variables.push_back(*(w->var("Bmumumass")));
  variables.push_back(*(w->var("Bmu1eta")));
  variables.push_back(*(w->var("Bmu2eta")));
  variables.push_back(*(w->var("Bmu1pt")));
  variables.push_back(*(w->var("Bmu2pt")));
  variables.push_back(*(w->var("Bmu1dxyPV")));
  variables.push_back(*(w->var("Bmu2dxyPV")));
  variables.push_back(*(w->var("Bmu1dzPV")));
  variables.push_back(*(w->var("Bmu2dzPV")));
  variables.push_back(*(w->var("Bd0")));
  variables.push_back(*(w->var("Bd0Err")));
  variables.push_back(*(w->var("Bdtheta")));
  variables.push_back(*(w->var("Balpha")));
  variables.push_back(*(w->var("Btrk1Dz1")));
  variables.push_back(*(w->var("Btrk1DzError1")));
  variables.push_back(*(w->var("Btrk1Dxy1")));
  variables.push_back(*(w->var("Btrk1DxyError1")));
  variables.push_back(*(w->var("Bmumueta")));
  variables.push_back(*(w->var("Bmumuphi")));
  variables.push_back(*(w->var("Bmumupt")));
  
  /*
  if(particle == 0){variables.push_back(*(w->var("Btrk1Dz1")));
    variables.push_back(*(w->var("Btrk1DzError1")));
    variables.push_back(*(w->var("Btrk1Dxy1")));
    variables.push_back(*(w->var("Btrk1DxyError1")));
    variables.push_back(*(w->var("Bd0")));
    variables.push_back(*(w->var("Bd0Err")));
    variables.push_back(*(w->var("BDT_pt_5_7")));
    variables.push_back(*(w->var("BDT_pt_7_10")));
    variables.push_back(*(w->var("BDT_pt_10_15")));
    variables.push_back(*(w->var("BDT_pt_15_20")));
    variables.push_back(*(w->var("BDT_pt_20_30")));
    variables.push_back(*(w->var("BDT_pt_30_50")));
    variables.push_back(*(w->var("BDT_pt_50_100")));
    variables.push_back(*(w->var("BDT_total")));
   }
  */
  
  if(particle == 1){variables.push_back(*(w->var("BDT_pt_5_10")));
    variables.push_back(*(w->var("BDT_pt_10_15")));
    variables.push_back(*(w->var("BDT_pt_15_20")));
    variables.push_back(*(w->var("BDT_pt_20_50")));
    //variables.push_back(*(w->var("BDT_total"))); 
    } 
  
  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central; 

  double left = particle ? 5.3  : 5.15;
  double right = particle ? 5.45 : 5.4;


  reduceddata_side = particle ? (RooDataSet*)data->reduce(Form("Bmass>%lf || Bmass<%lf", right, left)) : (RooDataSet*)data->reduce(Form("Bmass>%lf", right));
  reduceddata_central = (RooDataSet*)data->reduce(Form("Bmass>%lf",left));
  reduceddata_central = (RooDataSet*)reduceddata_central->reduce(Form("Bmass<%lf",right));

  //Integrating the background distribution

  RooAbsReal* int_fit_side_right = fit_side->createIntegral(variables[0], "right");
  RooAbsReal* int_fit_side_left = fit_side->createIntegral(variables[0], "left");
  RooAbsReal* int_fit_peak = fit_side->createIntegral(variables[0], "peak");

  std::cout<< std::endl << "Integral right band: " << int_fit_side_right->getVal() << std::endl;
  if(particle == 1){
    cout << "Integral left band: " << int_fit_side_left->getVal() << endl;
  }

  double factor = particle ? (int_fit_peak->getVal())/(int_fit_side_right->getVal() + int_fit_side_left->getVal()) : (int_fit_peak->getVal())/(int_fit_side_right->getVal()); 

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
    histos.push_back(create_histogram(variables[9], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[27]));
    
  }else if(particle == 1){
    histos.push_back(create_histogram(variables[1], "By",factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "Bpt",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1Pt",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Eta",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk1PtErr",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "BsvpvDisErr",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Bmu1dxyPV",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "Bmu2dxyPV",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "Bmu1dzPV",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Bmu2dzPV",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Bd0Err",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bdtheta",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[20]));
    histos.push_back(create_histogram(variables[22], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[21]));
    histos.push_back(create_histogram(variables[23], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[22]));
    histos.push_back(create_histogram(variables[24], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[23]));
    histos.push_back(create_histogram(variables[25], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[24]));
    histos.push_back(create_histogram(variables[26], "Bmumueta",factor, reduceddata_side, reduceddata_central, data, n[25]));
    histos.push_back(create_histogram(variables[27], "Bmumuphi",factor, reduceddata_side, reduceddata_central, data, n[26]));
    histos.push_back(create_histogram(variables[28], "Bmumupt",factor, reduceddata_side, reduceddata_central, data, n[27]));
    histos.push_back(create_histogram(variables[29], "BDT_pt_5_10",factor, reduceddata_side, reduceddata_central, data, n[28]));
    histos.push_back(create_histogram(variables[30], "BDT_pt_10_15",factor, reduceddata_side, reduceddata_central, data, n[29]));    
    histos.push_back(create_histogram(variables[31], "BDT_pt_15_20",factor, reduceddata_side, reduceddata_central, data, n[30]));    
    histos.push_back(create_histogram(variables[32], "BDT_pt_20_50",factor, reduceddata_side, reduceddata_central, data, n[31])); 
    //histos.push_back(create_histogram(variables[20], "BDT_total",factor, reduceddata_side, reduceddata_central, data, n[29]));    
     
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
  dist_peak->SetNameTitle(var.GetName(), "");

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
  leg->AddEntry(var.GetName(), "Signal", "l");
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

  return dist_peak;

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
  histo_Bp_sig->GetXaxis()->SetTitle(label );

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

  return histo_Bp_sig;

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
    h[i]->Fit("gaus","","",-20,20);//
    h[i]->GetFunction("gaus")->SetLineColor(4);
    h[i]->GetFunction("gaus")->SetLineWidth(5);
    h[i]->GetXaxis()->SetTitle("Pull");
    h[i]->GetYaxis()->SetTitle("Toy MCs");
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
  }
  
}

void AddWeights(TTree* t){
  TString input_file;
  input_file = particle ? "./results/Bs/mc_validation_plots/weights/weights.root": "./results/Bu/mc_validation_plots/weights/weights.root";
  TFile* f_wei = new TFile(input_file, "read");
 
  TH1D* h_pt = (TH1D*)f_wei->Get(TString("weights_Bpt"));
  TH1D* h_y = (TH1D*)f_wei->Get(TString("weights_By"));
  TH1D* h_trk1Pt = (TH1D*)f_wei->Get(TString("weights_Btrk1Pt"));
  TH1D* h_trk1Eta = (TH1D*)f_wei->Get(TString("weights_Btrk1Eta"));
  TH1D* h_trk1PtErr = (TH1D*)f_wei->Get(TString("weights_Btrk1PtErr"));
  TH1D* h_chi2cl = (TH1D*)f_wei->Get(TString("weights_Bchi2cl"));
  TH1D* h_svpvDistance = (TH1D*)f_wei->Get(TString("weights_BsvpvDistance"));
  TH1D* h_svpvDisErr = (TH1D*)f_wei->Get(TString("weights_BsvpvDisErr"));
  TH1D* h_mumumass = (TH1D*)f_wei->Get(TString("weights_Bmumumass"));
  TH1D* h_mu1eta = (TH1D*)f_wei->Get(TString("weights_Bmu1eta"));
  TH1D* h_mu2eta = (TH1D*)f_wei->Get(TString("weights_Bmu2eta"));
  TH1D* h_mu1pt = (TH1D*)f_wei->Get(TString("weights_Bmu1pt"));
  TH1D* h_mu2pt = (TH1D*)f_wei->Get(TString("weights_Bmu2pt"));
  TH1D* h_mu1dxyPV = (TH1D*)f_wei->Get(TString("weights_Bmu1dxyPV"));
  TH1D* h_mu2dxyPV = (TH1D*)f_wei->Get(TString("weights_Bmu2dxyPV"));
  TH1D* h_mu1dzPV = (TH1D*)f_wei->Get(TString("weights_Bmu1dzPV"));
  TH1D* h_mu2dzPV = (TH1D*)f_wei->Get(TString("weights_Bmu2dzPV"));
  TH1D* h_d0 = (TH1D*)f_wei->Get(TString("weights_Bd0"));
  TH1D* h_d0Err = (TH1D*)f_wei->Get(TString("weights_Bd0Err"));
  TH1D* h_dtheta = (TH1D*)f_wei->Get(TString("weights_Bdtheta"));
  TH1D* h_alpha = (TH1D*)f_wei->Get(TString("weights_Balpha"));
  TH1D* h_trk1Dz1 = (TH1D*)f_wei->Get(TString("weights_Btrk1Dz1"));
  TH1D* h_trk1DzError1 = (TH1D*)f_wei->Get(TString("weights_Btrk1DzError1"));
  TH1D* h_trk1Dxy1 = (TH1D*)f_wei->Get(TString("weights_Btrk1Dxy1"));
  TH1D* h_trk1DxyError1 = (TH1D*)f_wei->Get(TString("weights_Btrk1DxyError1"));
  TH1D* h_mumueta = (TH1D*)f_wei->Get(TString("weights_Bmumueta"));
  TH1D* h_mumuphi = (TH1D*)f_wei->Get(TString("weights_Bmumuphi"));
  TH1D* h_mumupt = (TH1D*)f_wei->Get(TString("weights_Bmumupt"));

  double pt_min = h_pt->GetXaxis()->GetXmin();
  double y_min = h_y->GetXaxis()->GetXmin();
  double trk1Pt_min = h_trk1Pt->GetXaxis()->GetXmin();
  double trk1Eta_min = h_trk1Eta->GetXaxis()->GetXmin();
  double trk1PtErr_min = h_trk1PtErr->GetXaxis()->GetXmin();
  double chi2cl_min = h_chi2cl->GetXaxis()->GetXmin();
  double svpvDistance_min = h_svpvDistance->GetXaxis()->GetXmin();
  double svpvDisErr_min = h_svpvDisErr->GetXaxis()->GetXmin();
  double mumumass_min = h_mumumass->GetXaxis()->GetXmin();
  double mu1eta_min = h_mu1eta->GetXaxis()->GetXmin();
  double mu2eta_min = h_mu2eta->GetXaxis()->GetXmin();
  double mu1pt_min = h_mu1pt->GetXaxis()->GetXmin();
  double mu2pt_min = h_mu2pt->GetXaxis()->GetXmin();
  double mu1dxyPV_min = h_mu1dxyPV->GetXaxis()->GetXmin();
  double mu2dxyPV_min = h_mu2dxyPV->GetXaxis()->GetXmin();
  double mu1dzPV_min = h_mu1dzPV->GetXaxis()->GetXmin();
  double mu2dzPV_min = h_mu2dzPV->GetXaxis()->GetXmin();
  double d0_min = h_d0->GetXaxis()->GetXmin();
  double d0Err_min = h_d0Err->GetXaxis()->GetXmin();
  double dtheta_min = h_dtheta->GetXaxis()->GetXmin();
  double alpha_min = h_alpha->GetXaxis()->GetXmin();
  double trk1Dz1_min = h_trk1Dz1->GetXaxis()->GetXmin();
  double trk1DzError1_min = h_trk1DzError1->GetXaxis()->GetXmin();
  double trk1Dxy1_min = h_trk1Dxy1->GetXaxis()->GetXmin();
  double trk1DxyError1_min = h_trk1DxyError1->GetXaxis()->GetXmin();
  double mumueta_min = h_mumueta->GetXaxis()->GetXmin();
  double mumuphi_min = h_mumuphi->GetXaxis()->GetXmin();
  double mumupt_min = h_mumupt->GetXaxis()->GetXmin();

  double pt_max = h_pt->GetXaxis()->GetXmax();
  double y_max = h_y->GetXaxis()->GetXmax();
  double trk1Pt_max = h_trk1Pt->GetXaxis()->GetXmax();
  double trk1Eta_max = h_trk1Eta->GetXaxis()->GetXmax();
  double trk1PtErr_max = h_trk1PtErr->GetXaxis()->GetXmax();
  double chi2cl_max = h_chi2cl->GetXaxis()->GetXmax();
  double svpvDistance_max = h_svpvDistance->GetXaxis()->GetXmax();
  double svpvDisErr_max = h_svpvDisErr->GetXaxis()->GetXmax();
  double mumumass_max = h_mumumass->GetXaxis()->GetXmax();
  double mu1eta_max = h_mu1eta->GetXaxis()->GetXmax();
  double mu2eta_max = h_mu2eta->GetXaxis()->GetXmax();
  double mu1pt_max = h_mu1pt->GetXaxis()->GetXmax();
  double mu2pt_max = h_mu2pt->GetXaxis()->GetXmax();
  double mu1dxyPV_max = h_mu1dxyPV->GetXaxis()->GetXmax();
  double mu2dxyPV_max = h_mu2dxyPV->GetXaxis()->GetXmax();
  double mu1dzPV_max = h_mu1dzPV->GetXaxis()->GetXmax();
  double mu2dzPV_max = h_mu2dzPV->GetXaxis()->GetXmax();
  double d0_max = h_d0->GetXaxis()->GetXmax();
  double d0Err_max = h_d0Err->GetXaxis()->GetXmax();
  double dtheta_max = h_dtheta->GetXaxis()->GetXmax();
  double alpha_max = h_alpha->GetXaxis()->GetXmax();
  double trk1Dz1_max = h_trk1Dz1->GetXaxis()->GetXmax();
  double trk1DzError1_max = h_trk1DzError1->GetXaxis()->GetXmax();
  double trk1Dxy1_max = h_trk1Dxy1->GetXaxis()->GetXmax();
  double trk1DxyError1_max = h_trk1DxyError1->GetXaxis()->GetXmax();
  double mumueta_max = h_mumueta->GetXaxis()->GetXmax();
  double mumuphi_max = h_mumuphi->GetXaxis()->GetXmax();
  double mumupt_max = h_mumupt->GetXaxis()->GetXmax();
 
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

  Float_t weight_Bpt;
  Float_t weight_By;
  Float_t weight_Btrk1Pt;
  Float_t weight_Btrk1Eta;
  Float_t weight_Btrk1PtErr;
  Float_t weight_Bchi2cl;
  Float_t weight_BsvpvDistance;
  Float_t weight_BsvpvDisErr;
  Float_t weight_Bmumumass;
  Float_t weight_Bmu1eta;
  Float_t weight_Bmu2eta;
  Float_t weight_Bmu1pt;
  Float_t weight_Bmu2pt;
  Float_t weight_Bmu1dxyPV;
  Float_t weight_Bmu2dxyPV;
  Float_t weight_Bmu1dzPV;
  Float_t weight_Bmu2dzPV;
  Float_t weight_Bd0;
  Float_t weight_Bd0Err;
  Float_t weight_Bdtheta;
  Float_t weight_Balpha;
  Float_t weight_Btrk1Dz1;
  Float_t weight_Btrk1DzError1;
  Float_t weight_Btrk1Dxy1;
  Float_t weight_Btrk1DxyError1;
  Float_t weight_Bmumueta;
  Float_t weight_Bmumuphi;
  Float_t weight_Bmumupt;

  t->SetBranchAddress("Bpt", &Bpt);
  t->SetBranchAddress("By", &By);
  t->SetBranchAddress("Btrk1Pt", &Btrk1Pt);
  t->SetBranchAddress("Btrk1Eta", &Btrk1Eta);
  t->SetBranchAddress("Btrk1PtErr", &Btrk1PtErr);
  t->SetBranchAddress("Bchi2cl", &Bchi2cl);
  t->SetBranchAddress("BsvpvDistance", &BsvpvDistance);
  t->SetBranchAddress("BsvpvDisErr", &BsvpvDisErr);
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

  tw->Branch("weight_Bpt", &weight_Bpt);
  tw->Branch("weight_By", &weight_By);
  tw->Branch("weight_Btrk1Pt", &weight_Btrk1Pt);
  tw->Branch("weight_Btrk1Eta", &weight_Btrk1Eta);
  tw->Branch("weight_Btrk1PtErr", &weight_Btrk1PtErr);
  tw->Branch("weight_Bchi2cl", &weight_Bchi2cl);
  tw->Branch("weight_BsvpvDistance", &weight_BsvpvDistance);
  tw->Branch("weight_BsvpvDisErr", &weight_BsvpvDisErr);
  tw->Branch("weight_Bmumumass", &weight_Bmumumass);
  tw->Branch("weight_Bmu1eta", &weight_Bmu1eta);
  tw->Branch("weight_Bmu2eta", &weight_Bmu2eta);
  tw->Branch("weight_Bmu1pt", &weight_Bmu1pt);
  tw->Branch("weight_Bmu2pt", &weight_Bmu2pt);
  tw->Branch("weight_Bmu1dxyPV", &weight_Bmu1dxyPV);
  tw->Branch("weight_Bmu2dxyPV", &weight_Bmu2dxyPV);
  tw->Branch("weight_Bmu1dzPV", &weight_Bmu1dzPV);
  tw->Branch("weight_Bmu2dzPV", &weight_Bmu2dzPV);
  tw->Branch("weight_Bd0", &weight_Bd0);
  tw->Branch("weight_Bd0Err", &weight_Bd0Err);
  tw->Branch("weight_Bdtheta", &weight_Bdtheta);
  tw->Branch("weight_Balpha", &weight_Balpha);
  tw->Branch("weight_Btrk1Dz1", &weight_Btrk1Dz1);
  tw->Branch("weight_Btrk1DzError1", &weight_Btrk1DzError1);
  tw->Branch("weight_Btrk1Dxy1", &weight_Btrk1Dxy1);
  tw->Branch("weight_Btrk1DxyError1", &weight_Btrk1DxyError1);
  tw->Branch("weight_Bmumueta", &weight_Bmumueta);
  tw->Branch("weight_Bmumuphi", &weight_Bmumuphi);
  tw->Branch("weight_Bmumupt", &weight_Bmumupt);

  cout << "Starting Cycle" << endl;

  for (int i = 0; i < (t->GetEntries()); i++){
     t->GetEntry(i);

     if(Bpt>=pt_min && Bpt<=pt_max){
       weight_Bpt = h_pt->GetBinContent(h_pt->FindBin(Bpt));
     }
     if(By>=y_min && By<=y_max){
       weight_By = h_y->GetBinContent(h_y->FindBin(By));
     }
     if(Btrk1Pt>=trk1Pt_min && Btrk1Pt<=trk1Pt_max){
       weight_Btrk1Pt = h_trk1Pt->GetBinContent(h_trk1Pt->FindBin(Btrk1Pt));
     }
     if(Btrk1Eta>=trk1Eta_min && Btrk1Eta<=trk1Eta_max){
       weight_Btrk1Eta = h_trk1Eta->GetBinContent(h_trk1Eta->FindBin(Btrk1Eta));
     }
     if(Btrk1PtErr>=trk1PtErr_min && Btrk1PtErr<=trk1PtErr_max){
       weight_Btrk1PtErr = h_trk1PtErr->GetBinContent(h_trk1PtErr->FindBin(Btrk1PtErr));
     }
     if(Bchi2cl>=chi2cl_min && Bchi2cl<=chi2cl_max){
       weight_Bchi2cl = h_chi2cl->GetBinContent(h_chi2cl->FindBin(Bchi2cl));
     }
     if(BsvpvDistance>=svpvDistance_min && BsvpvDistance<=svpvDistance_max){
       weight_BsvpvDistance = h_svpvDistance->GetBinContent(h_svpvDistance->FindBin(BsvpvDistance));
     }
     if(BsvpvDisErr>=svpvDisErr_min && BsvpvDisErr<=svpvDisErr_max){
       weight_BsvpvDisErr = h_svpvDisErr->GetBinContent(h_svpvDisErr->FindBin(BsvpvDisErr));
     }
     if(Bmumumass>=mumumass_min && Bmumumass<=mumumass_max){
       weight_Bmumumass = h_mumumass->GetBinContent(h_mumumass->FindBin(Bmumumass));
     }
     if(Bmu1eta>=mu1eta_min && Bmu1eta<=mu1eta_max){
       weight_Bmu1eta = h_mu1eta->GetBinContent(h_mu1eta->FindBin(Bmu1eta));
     }
     if(Bmu2eta>=mu2eta_min && Bmu2eta<=mu2eta_max){
       weight_Bmu2eta = h_mu2eta->GetBinContent(h_mu2eta->FindBin(Bmu2eta));
     }
     if(Bmu1pt>=mu1pt_min && Bmu1pt<=mu1pt_max){
       weight_Bmu1pt = h_mu1pt->GetBinContent(h_mu1pt->FindBin(Bmu1pt));
     }
     if(Bmu2pt>=mu2pt_min && Bmu2pt<=mu2pt_max){
       weight_Bmu2pt = h_mu2pt->GetBinContent(h_mu2pt->FindBin(Bmu2pt));
     }
     if(Bmu1dxyPV>=mu1dxyPV_min && Bmu1dxyPV<=mu1dxyPV_max){
       weight_Bmu1dxyPV = h_mu1dxyPV->GetBinContent(h_mu1dxyPV->FindBin(Bmu1dxyPV));
     }
     if(Bmu2dxyPV>=mu2dxyPV_min && Bmu2dxyPV<=mu2dxyPV_max){
       weight_Bmu2dxyPV = h_mu2dxyPV->GetBinContent(h_mu2dxyPV->FindBin(Bmu2dxyPV));
     }
     if(Bmu1dzPV>=mu1dzPV_min && Bmu1dzPV<=mu1dzPV_max){
       weight_Bmu1dzPV = h_mu1dzPV->GetBinContent(h_mu1dzPV->FindBin(Bmu1dzPV));
     }
     if(Bmu2dzPV>=mu2dzPV_min && Bmu2dzPV<=mu2dzPV_max){
       weight_Bmu2dzPV = h_mu2dzPV->GetBinContent(h_mu2dzPV->FindBin(Bmu2dzPV));
     }
     if(Bd0>=d0_min && Bd0<=d0_max){
       weight_Bd0 = h_d0->GetBinContent(h_d0->FindBin(Bd0));
     }
     if(Bd0Err>=d0Err_min && Bd0Err<=d0Err_max){
       weight_Bd0Err = h_d0Err->GetBinContent(h_d0Err->FindBin(Bd0Err));
     }
     if(Bdtheta>=dtheta_min && Bdtheta<=dtheta_max){
       weight_Bdtheta = h_dtheta->GetBinContent(h_dtheta->FindBin(Bdtheta));
     }
     if(Balpha>=alpha_min && Balpha<=alpha_max){
       weight_Balpha = h_alpha->GetBinContent(h_alpha->FindBin(Balpha));
     }
     if(Btrk1Dz1>=trk1Dz1_min && Btrk1Dz1<=trk1Dz1_max){
       weight_Btrk1Dz1 = h_trk1Dz1->GetBinContent(h_trk1Dz1->FindBin(Btrk1Dz1));
     }
     if(Btrk1DzError1>=trk1DzError1_min && Btrk1DzError1<=trk1DzError1_max){
       weight_Btrk1DzError1 = h_trk1DzError1->GetBinContent(h_trk1DzError1->FindBin(Btrk1DzError1));
     }
     if(Btrk1Dxy1>=trk1Dxy1_min && Btrk1Dxy1<=trk1Dxy1_max){
       weight_Btrk1Dxy1 = h_trk1Dxy1->GetBinContent(h_trk1Dxy1->FindBin(Btrk1Dxy1));
     }
     if(Btrk1DxyError1>=trk1DxyError1_min && Btrk1DxyError1<=trk1DxyError1_max){
       weight_Btrk1DxyError1 = h_trk1DxyError1->GetBinContent(h_trk1DxyError1->FindBin(Btrk1DxyError1));
     }
     if(Bmumueta>=mumueta_min && Bmumueta<=mumueta_max){
       weight_Bmumueta = h_mumueta->GetBinContent(h_mumueta->FindBin(Bmumueta));
     }
     if(Bmumuphi>=mumuphi_min && Bmumuphi<=mumuphi_max){
       weight_Bmumuphi = h_mumuphi->GetBinContent(h_mumuphi->FindBin(Bmumuphi));
     }
     if(Bmumupt>=mumupt_min && Bmumupt<=mumupt_max){
       weight_Bmumupt = h_mumupt->GetBinContent(h_mumupt->FindBin(Bmumupt));
     }
     else{
       weight_Bpt = 1.;
       weight_By = 1.;
       weight_Btrk1Pt = 1.;
       weight_Btrk1Eta = 1.;
       weight_Btrk1PtErr = 1.;
       weight_Bchi2cl = 1.;
       weight_BsvpvDistance = 1.;
       weight_BsvpvDisErr = 1.;
       weight_Bmumumass = 1.;
       weight_Bmu1eta = 1.;
       weight_Bmu2eta = 1.;
       weight_Bmu1pt = 1.;
       weight_Bmu2pt = 1.;
       weight_Bmu1dxyPV = 1.;
       weight_Bmu2dxyPV = 1.;
       weight_Bmu1dzPV = 1.;
       weight_Bmu2dzPV = 1.;
       weight_Bd0 = 1.;
       weight_Bd0Err = 1.;
       weight_Bdtheta = 1.;
       weight_Balpha = 1.;
       weight_Btrk1Dz1 = 1.;
       weight_Btrk1DzError1 = 1.;
       weight_Btrk1Dxy1 = 1.;
       weight_Btrk1DxyError1 = 1.;
       weight_Bmumueta = 1.;
       weight_Bmumuphi = 1.;
       weight_Bmumupt = 1.;

     }
     tw->Fill();
  }

  cout << "Ending Cycle" << endl;

  f_tree->Write();
  tw->Show();
  f_tree->Close();
  f_wei->Close();
}
//Weights were added to TTree

void set_up_workspace_variables(RooWorkspace& w)
{

  if(particle == 0){

    double mass_min, mass_max;
    double size_min, size_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double trk1pt_min, trk1pt_max;
    double trk1eta_min, trk1eta_max;
    double trk1pterr_min, trk1pterr_max;
    double chi2cl_min, chi2cl_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDistanceErr_min, svpvDistanceErr_max;
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

    mass_min = 5.;
    mass_max = 6.;

    size_min = 0.;
    size_max = 15.;    

    y_min = -2.4;
    y_max = 2.4;

    pt_min = 5.;
    pt_max = 40.;

    trk1pt_min = 5.;
    trk1pt_max = 15.;

    trk1eta_min = -2.4;
    trk1eta_max = 2.4;

    trk1pterr_min = 0.;  
    trk1pterr_max = 0.5;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 22.;

    svpvDistanceErr_min = 0.;
    svpvDistanceErr_max = 0.2;

    mumumass_min = 2.95;
    mumumass_max = 3.25; 
 
    mu1eta_min = -2.4;
    mu1eta_max = 2.4;

    mu2eta_min = -2.4;
    mu2eta_max = 2.4;

    mu1pt_min = 5.;
    mu1pt_max = 20.;

    mu2pt_min = 5.;
    mu2pt_max = 20.;

    mu1dxyPV_min = -0.3; 
    mu1dxyPV_max = 0.3;

    mu2dxyPV_min = -0.3;
    mu2dxyPV_max = 0.3;

    mu1dzPV_min = -20.; 
    mu1dzPV_max = 20.;

    mu2dzPV_min = -20.;
    mu2dzPV_max = 20.;

    d0_min = 0.; 
    d0_max = 0.5;

    d0err_min = 0.;
    d0err_max = 0.0005;

    dtheta_min = 0.;
    dtheta_max = 3.5;

    alpha_min = 0.;
    alpha_max = 3.5;

    trk1Dz1_min = -50.;
    trk1Dz1_max = 50.;

    trk1DzError1_min = -0.1;
    trk1DzError1_max = 0.1;

    trk1Dxy1_min = -10.;
    trk1Dxy1_max = 10.;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.1;

    mumueta_min = -3.5;
    mumueta_max = 3.5;

    mumuphi_min = -3.5;
    mumuphi_max = 3.5;

    mumupt_min = 0.;
    mumupt_max = 50.;


    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);  
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
    RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDistanceErr_min, svpvDistanceErr_max);
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
    
  }
    
  
  if(particle == 1){
    double mass_min, mass_max;
    double size_min, size_max;
    double y_min, y_max;
    double pt_min, pt_max;
    double trk1pt_min, trk1pt_max;
    double trk2pt_min, trk2pt_max;
    double trk1eta_min, trk1eta_max;
    double trk1pterr_min, trk1pterr_max;
    double chi2cl_min, chi2cl_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDistanceErr_min, svpvDistanceErr_max;
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
    
    double BDT_5_10_min, BDT_5_10_max;
    double BDT_10_15_min, BDT_10_15_max;
    double BDT_15_20_min, BDT_15_20_max;
    double BDT_20_50_min, BDT_20_50_max;
    double BDT_total_min, BDT_total_max;
 
    mass_min = 5.0;
    mass_max = 6.0;

    size_min = 0.5;
    size_max = 3.5;

    y_min = -2.4;
    y_max = 2.4;

    pt_min = 5.;
    pt_max = 100.;

    trk1pt_min = 1.5;
    trk1pt_max = 17.5;

    trk2pt_min = 1.5;
    trk2pt_max = 17.5;

    trk1eta_min = -2.3;
    trk1eta_max = 2.3;

    trk1pterr_min = 0.;
    trk1pterr_max = 0.2;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    svpvDistance_min = 0.;
    svpvDistance_max = 22.;

    svpvDistanceErr_min = 0.;
    svpvDistanceErr_max = 0.1;
   
    mumumass_min = 2.95;
    mumumass_max = 3.2;

    mu1eta_min = -2.3;
    mu1eta_max = 2.3;

    mu2eta_min = -2.3;
    mu2eta_max = 2.3;

    mu1pt_min = 2.;
    mu1pt_max = 31.;

    mu2pt_min = 2.;
    mu2pt_max = 40.;

    mu1dxyPV_min = -0.2;
    mu1dxyPV_max = 0.21;

    mu2dxyPV_min = -0.2;
    mu2dxyPV_max = 0.2;

    mu1dzPV_min = -25.;
    mu1dzPV_max = 20.;

    mu2dzPV_min = -23.;
    mu2dzPV_max = 20.;

    d0_min = 0.;
    d0_max = 1.;

    d0err_min = 0.;
    d0err_max = 0.0004;

    dtheta_min = 0.;
    dtheta_max = 3.5;

    alpha_min = 0.;
    alpha_max = 3.3;

    trk1Dz1_min = -25.;
    trk1Dz1_max = 20.;

    trk1DzError1_min = -0.1;
    trk1DzError1_max = 0.1;

    trk1Dxy1_min = -0.2;
    trk1Dxy1_max = 0.25;

    trk1DxyError1_min = 0.;
    trk1DxyError1_max = 0.04;

    mumueta_min = -3.5;
    mumueta_max = 4.0;

    mumuphi_min = -3.5;
    mumuphi_max = 3.5;

    mumupt_min = 0.;
    mumupt_max = 80.;

    BDT_5_10_min = -0.4;
    BDT_5_10_max = 0.8;

    BDT_10_15_min = -0.3;
    BDT_10_15_max = 0.8;

    BDT_15_20_min = -0.2;
    BDT_15_20_max = 0.8;

    BDT_20_50_min = 0.05;
    BDT_20_50_max = 0.6;

    //BDT_total_min = 0.29;
    // BDT_total_max = DATA_CUT? 0.70 : 0.73;


    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar Btrk1Pt("Btrk1Pt","Btrk1Pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk1Eta("Btrk1Eta","Btrk1Eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1PtErr("Btrk1PtErr","Btrk1PtErr",trk1pterr_min,trk1pterr_max);
    RooRealVar Bchi2cl("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDisErr("BsvpvDisErr", "BsvpvDisErr", svpvDistanceErr_min, svpvDistanceErr_max);
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
    RooRealVar Bmumueta("Bmumueta", "Bmumueta", mumueta_min, mumueta_max);
    RooRealVar Bmumuphi("Bmumuphi", "Bmumuphi", mumuphi_min, mumuphi_max);
    RooRealVar Bmumupt("Bmumupt", "Bmumupt", mumupt_min, mumupt_max);
    RooRealVar BDT_pt_5_10("BDT_pt_5_10", "BDT_pt_5_10", BDT_5_10_min, BDT_5_10_max);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", BDT_10_15_min, BDT_10_15_max);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", BDT_15_20_min, BDT_15_20_max);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", BDT_20_50_min, BDT_20_50_max);
    //RooRealVar BDT_total("BDT_total", "BDT_total", BDT_total_min, BDT_total_max);
    

    w.import(Bmass);
    w.import(By);
    w.import(Bpt);
    w.import(Btrk1Pt);
    w.import(Btrk1Eta);
    w.import(Btrk1PtErr);
    w.import(Bchi2cl);
    w.import(BsvpvDistance);
    w.import(BsvpvDisErr);
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
    w.import(BDT_pt_5_10);
    w.import(BDT_pt_10_15);
    w.import(BDT_pt_15_20);
    w.import(BDT_pt_20_50);
    //w.import(BDT_total);
          
  }
}
