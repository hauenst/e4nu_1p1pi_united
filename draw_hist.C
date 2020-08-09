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

double N_E_bins;

void draw_hist(){
  TFile *file_in;
  TFile *file_in1;

  file_in = new TFile("56Fe_4461.root");
  file_in1 = new TFile("C12_2261.root");

TH1F *sub_cal_all_pimi = (TH1F*)file_in->Get("sub_cal_all_pimi");
TH1F *sub_kin_e_all_pimi = (TH1F*)file_in->Get("sub_kin_e_all_pimi");
TH1F *sub_cal_all_pipl = (TH1F*)file_in->Get("sub_cal_all_pipl");
TH1F *sub_kin_e_all_pipl = (TH1F*)file_in->Get("sub_kin_e_all_pipl");
TH2F *mult_1p = (TH2F*)file_in->Get("h2_phot_pi_1p");
TH2F *mult_2p = (TH2F*)file_in->Get("h2_phot_pi_2p");
TH2F *mult_3p = (TH2F*)file_in->Get("h2_phot_pi_3p");
TH1F *process1 = (TH1F*)file_in1->Get("process1");
TH1F *process2 = (TH1F*)file_in1->Get("process2");
TH1F *process3 = (TH1F*)file_in1->Get("process3");
TH1F *process4 = (TH1F*)file_in1->Get("process4");
TH1F *process5 = (TH1F*)file_in1->Get("process5");
TH1F *h1_cal_p_slice1 = (TH1F*)file_in1->Get("h1_cal_p_slice1");
TH1F *h1_cal_p_slice2 = (TH1F*)file_in1->Get("h1_cal_p_slice2");
TH1F *h1_cal_p_slice3 = (TH1F*)file_in1->Get("h1_cal_p_slice3");
TH1F *h1_cal_p_slice1_sub = (TH1F*)file_in1->Get("h1_cal_p_slice1_sub");
TH1F *h1_cal_p_slice2_sub = (TH1F*)file_in1->Get("h1_cal_p_slice2_sub");
TH1F *h1_cal_p_slice3_sub = (TH1F*)file_in1->Get("h1_cal_p_slice3_sub");

gStyle->SetOptStat(0);

N_E_bins = sub_cal_all_pimi->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  sub_cal_all_pimi->SetBinContent(i,sub_cal_all_pimi->GetBinContent(i)/sub_cal_all_pimi->GetBinWidth(i));
}

N_E_bins = sub_kin_e_all_pimi->GetNbinsX();

for(int j = 1; j<=N_E_bins;j++){
  sub_kin_e_all_pimi->SetBinContent(j,sub_kin_e_all_pimi->GetBinContent(j)/sub_kin_e_all_pimi->GetBinWidth(j));
}

N_E_bins = sub_cal_all_pipl->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  sub_cal_all_pipl->SetBinContent(i,sub_cal_all_pipl->GetBinContent(i)/sub_cal_all_pipl->GetBinWidth(i));
}

N_E_bins = sub_kin_e_all_pipl->GetNbinsX();

for(int j = 1; j<=N_E_bins;j++){
  sub_kin_e_all_pipl->SetBinContent(j,sub_kin_e_all_pipl->GetBinContent(j)/sub_kin_e_all_pipl->GetBinWidth(j));
}

N_E_bins = h1_cal_p_slice1->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice1->SetBinContent(i,h1_cal_p_slice1->GetBinContent(i)/h1_cal_p_slice1->GetBinWidth(i));
}

N_E_bins = h1_cal_p_slice2->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice2->SetBinContent(i,h1_cal_p_slice2->GetBinContent(i)/h1_cal_p_slice2->GetBinWidth(i));
}

N_E_bins = h1_cal_p_slice3->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice3->SetBinContent(i,h1_cal_p_slice3->GetBinContent(i)/h1_cal_p_slice3->GetBinWidth(i));
}

N_E_bins = h1_cal_p_slice1_sub->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice1_sub->SetBinContent(i,h1_cal_p_slice1_sub->GetBinContent(i)/h1_cal_p_slice1_sub->GetBinWidth(i));
}

N_E_bins = h1_cal_p_slice2_sub->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice2_sub->SetBinContent(i,h1_cal_p_slice2_sub->GetBinContent(i)/h1_cal_p_slice2_sub->GetBinWidth(i));
}

N_E_bins = h1_cal_p_slice3_sub->GetNbinsX();

for(int i = 1; i<=N_E_bins;i++){
  h1_cal_p_slice3_sub->SetBinContent(i,h1_cal_p_slice3_sub->GetBinContent(i)/h1_cal_p_slice3_sub->GetBinWidth(i));
}

TCanvas *c1 = new TCanvas("c1","",567,370);

sub_cal_all_pimi->SetAxisRange(0, 6, "X");
sub_cal_all_pimi->SetMinimum(0);
sub_kin_e_all_pimi->SetMinimum(0);
sub_cal_all_pimi->SetLineWidth(3);
sub_kin_e_all_pimi->SetLineWidth(3);
sub_kin_e_all_pimi->SetLineColor(1);
sub_cal_all_pimi->GetXaxis()->SetTitle("E_{rec}");
sub_cal_all_pimi->GetXaxis()->SetTitleSize(0.05);
sub_cal_all_pimi->GetXaxis()->SetLabelSize(0.05);
sub_cal_all_pimi->GetYaxis()->SetLabelSize(0.05);
sub_cal_all_pimi->GetXaxis()->SetTitleOffset(0.7);
sub_cal_all_pimi->Draw("HIST");
sub_kin_e_all_pimi->Draw("HIST SAME");


c1->Print("56Fe_4461_pimi.svg");

TCanvas *c2 = new TCanvas("c2","",567,370);

sub_cal_all_pipl->SetAxisRange(0, 6, "X");
sub_cal_all_pipl->SetMinimum(0);
sub_kin_e_all_pipl->SetMinimum(0);
sub_cal_all_pipl->SetLineWidth(3);
sub_kin_e_all_pipl->SetLineWidth(3);
sub_kin_e_all_pipl->SetLineColor(1);
sub_cal_all_pipl->GetXaxis()->SetTitle("E_{rec}");
sub_cal_all_pipl->GetXaxis()->SetTitleSize(0.05);
sub_cal_all_pipl->GetXaxis()->SetLabelSize(0.05);
sub_cal_all_pipl->GetYaxis()->SetLabelSize(0.05);
sub_cal_all_pipl->GetXaxis()->SetTitleOffset(0.7);
sub_cal_all_pipl->Draw("HIST");
sub_kin_e_all_pipl->Draw("HIST SAME");



c2->Print("56Fe_4461_pipl.svg");

TCanvas *c3 = new TCanvas("c3","",500,500);

mult_1p->GetXaxis()->SetTitle("#pi^{#pm} Multiplicity");
mult_1p->GetYaxis()->SetTitle("#gamma Multiplicity");
mult_1p->SetTitleOffset(1.05,"X");
mult_1p->SetTitleOffset(1.05,"Y");
mult_1p->GetXaxis()->SetTitleSize(0.05);
mult_1p->GetXaxis()->SetLabelSize(0.05);
mult_1p->GetYaxis()->SetTitleSize(0.05);
mult_1p->GetYaxis()->SetLabelSize(0.05);
mult_1p->SetMarkerColor(0);
mult_1p->Draw("COL TEXT");

c3->Print("56Fe_4461_mult_1p.svg");

TCanvas *c4 = new TCanvas("c4","",500,500);

mult_2p->GetXaxis()->SetTitle("#pi^{#pm} Multiplicity");
mult_2p->GetYaxis()->SetTitle("#gamma Multiplicity");
mult_2p->SetTitleOffset(1.05,"X");
mult_2p->SetTitleOffset(1.05,"Y");
mult_2p->GetXaxis()->SetTitleSize(0.05);
mult_2p->GetXaxis()->SetLabelSize(0.05);
mult_2p->GetYaxis()->SetTitleSize(0.05);
mult_2p->GetYaxis()->SetLabelSize(0.05);
mult_2p->SetMarkerColor(0);
mult_2p->Draw("COL TEXT");

c4->Print("56Fe_4461_mult_2p.svg");

TCanvas *c5 = new TCanvas("c5","",500,500);

mult_3p->GetXaxis()->SetTitle("#pi^{#pm} Multiplicity");
mult_3p->GetYaxis()->SetTitle("#gamma Multiplicity");
mult_3p->SetTitleOffset(1.05,"X");
mult_3p->SetTitleOffset(1.05,"Y");
mult_3p->GetXaxis()->SetTitleSize(0.05);
mult_3p->GetXaxis()->SetLabelSize(0.05);
mult_3p->GetYaxis()->SetTitleSize(0.05);
mult_3p->GetYaxis()->SetLabelSize(0.05);
mult_3p->SetMarkerColor(0);
mult_3p->Draw("COL TEXT");

c5->Print("56Fe_4461_mult_3p.svg");

TCanvas *c6 = new TCanvas("c6","",500,500);

process1->SetLineWidth(3);
process1->SetLineColor(1);
process1->Draw("HIST");
process2->SetLineWidth(3);
process2->SetLineColor(2);
process2->Draw("HIST SAME");
process3->SetLineWidth(3);
process3->SetLineColor(3);
process3->Draw("HIST SAME");
process4->SetLineWidth(3);
process4->SetLineColor(4);
process4->Draw("HIST SAME");
process5->SetLineWidth(3);
process5->SetLineColor(5);
process5->Draw("HIST SAME");

c6->Print("C12_2261_process.svg");

TCanvas *c7 = new TCanvas("c7","",567,370);

h1_cal_p_slice1->SetAxisRange(0, 3.5, "X");
h1_cal_p_slice1->SetMinimum(0);
h1_cal_p_slice1->SetLineWidth(3);
h1_cal_p_slice1_sub->SetLineWidth(3);
h1_cal_p_slice1_sub->SetLineColor(2);
h1_cal_p_slice1->GetXaxis()->SetTitle("E_{cal}");
h1_cal_p_slice1->GetXaxis()->SetTitleSize(0.05);
h1_cal_p_slice1->GetXaxis()->SetLabelSize(0.05);
h1_cal_p_slice1->GetYaxis()->SetLabelSize(0.05);
h1_cal_p_slice1->GetXaxis()->SetTitleOffset(0.7);
h1_cal_p_slice1->Draw("HIST");
h1_cal_p_slice1_sub->Draw("HIST SAME");

c7->Print("C12_2261_slice1.svg");

TCanvas *c8 = new TCanvas("c8","",567,370);

h1_cal_p_slice2->SetAxisRange(0, 3.5, "X");
h1_cal_p_slice2->SetMinimum(0);
h1_cal_p_slice2->SetLineWidth(3);
h1_cal_p_slice2_sub->SetLineWidth(3);
h1_cal_p_slice2_sub->SetLineColor(2);
h1_cal_p_slice2->GetXaxis()->SetTitle("E_{cal}");
h1_cal_p_slice2->GetXaxis()->SetTitleSize(0.05);
h1_cal_p_slice2->GetXaxis()->SetLabelSize(0.05);
h1_cal_p_slice2->GetYaxis()->SetLabelSize(0.05);
h1_cal_p_slice2->GetXaxis()->SetTitleOffset(0.7);
h1_cal_p_slice2->Draw("HIST");
h1_cal_p_slice2_sub->Draw("HIST SAME");

c8->Print("C12_2261_slice2.svg");

TCanvas *c9 = new TCanvas("c9","",567,370);

h1_cal_p_slice3->SetAxisRange(0, 3.5, "X");
h1_cal_p_slice3->SetMinimum(0);
h1_cal_p_slice3->SetLineWidth(3);
h1_cal_p_slice3_sub->SetLineWidth(3);
h1_cal_p_slice3_sub->SetLineColor(2);
h1_cal_p_slice3->GetXaxis()->SetTitle("E_{cal}");
h1_cal_p_slice3->GetXaxis()->SetTitleSize(0.05);
h1_cal_p_slice3->GetXaxis()->SetLabelSize(0.05);
h1_cal_p_slice3->GetYaxis()->SetLabelSize(0.05);
h1_cal_p_slice3->GetXaxis()->SetTitleOffset(0.7);
h1_cal_p_slice3->Draw("HIST");
h1_cal_p_slice3_sub->Draw("HIST SAME");

c9->Print("C12_2261_slice3.svg");
}
