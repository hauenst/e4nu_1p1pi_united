#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

using namespace std;

void hist_draw()
{
  TFile *file_in;
  TFile *file_in1;

  file_in = new TFile("data_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root");
  file_in1 = new TFile("C12_2261.root");


  TH1F *old_1p2pi = (TH1F*)file_in1->Get("2pi_1p_rot_cal_pimi");
  TH1F *new_1p2pi = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pimi");
  TH1F *old_1p1pi1phot = (TH1F*)file_in1->Get("1pi_1p_1phot_rot_cal_pimi");
  //old_1p2pi->Add(old_1p1pi1phot, 1);

  double N_E_bins = old_1p2pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_1p2pi->SetBinContent(i,old_1p2pi->GetBinContent(i)/old_1p2pi->GetBinWidth(i));
  }

  double N_E_bins1 = new_1p2pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins;j++){
    new_1p2pi->SetBinContent(j,new_1p2pi->GetBinContent(j)/new_1p2pi->GetBinWidth(j));
  }

  TCanvas *c1 = new TCanvas("c1","",567,370);

  old_1p2pi->GetXaxis()->SetTitle("E_{cal}");
  old_1p2pi->GetXaxis()->SetTitleSize(0.05);
  old_1p2pi->GetXaxis()->SetLabelSize(0.05);
  old_1p2pi->GetYaxis()->SetLabelSize(0.05);
  old_1p2pi->GetXaxis()->SetTitleOffset(0.7);
  old_1p2pi->Draw();

  c1->Print("old_1p2pi.svg");

  TCanvas *c2 = new TCanvas("c2","",567,370);

  new_1p2pi->GetXaxis()->SetTitle("E_{cal}");
  new_1p2pi->GetXaxis()->SetTitleSize(0.05);
  new_1p2pi->GetXaxis()->SetLabelSize(0.05);
  new_1p2pi->GetYaxis()->SetLabelSize(0.05);
  new_1p2pi->GetXaxis()->SetTitleOffset(0.7);
  new_1p2pi->Draw();

  c2->Print("new_1p2pi.svg");

  TCanvas *c3 = new TCanvas("c3","",567,370);


  old_1p2pi->GetXaxis()->SetTitle("E_{cal}");
  old_1p2pi->GetXaxis()->SetTitleSize(0.05);
  old_1p2pi->GetXaxis()->SetLabelSize(0.05);
  old_1p2pi->GetYaxis()->SetLabelSize(0.05);
  old_1p2pi->GetXaxis()->SetTitleOffset(0.7);
  new_1p2pi->SetLineColor(2);
  new_1p2pi->Draw();
  old_1p2pi->Draw("SAME");

  c3->Print("comparison.svg");
}
