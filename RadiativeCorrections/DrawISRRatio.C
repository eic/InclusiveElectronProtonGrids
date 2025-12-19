#include "ePIC_style.C"

void DrawISRRatio() {
  gROOT->ProcessLine("set_ePIC_style()");
  
  // Create a canvas with two pads (top for histograms, bottom for ratio)
  TCanvas *canvas = new TCanvas("c", "ISR Histograms with Ratio", 900, 600);
  canvas->Divide(1,2);
  
  // Adjust pad sizes
  TPad *pad1 = (TPad*)canvas->cd(1);
  pad1->SetPad(0,0.4,1,1);   // top pad
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.005);
  pad1->SetLogy();
  
  TPad *pad2 = (TPad*)canvas->cd(2);
  pad2->SetPad(0,0,1,0.4);   // bottom pad
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  
  // Load data
  TChain *rad_chain = new TChain("events");
  rad_chain->Add("new_track_ele_rad_output_tree.root");
  
  // Define histograms
  TH1F *hE_gamma_noCut   = new TH1F("hE_gamma_noCut",   ";E_{#gamma} [GeV];Events", 45, 0, 18);
  TH1F *hE_gamma_EmpzCut = new TH1F("hE_gamma_EmpzCut", ";E_{#gamma} [GeV];Events", 45, 0, 18);
  
  // Fill histograms
  rad_chain->Draw("E_gamma >> hE_gamma_noCut",
                  "weight*(y_ele>0.01 && y_ele<0.95 && Q2_ele>10)*(channel==1)");
  rad_chain->Draw("E_gamma >> hE_gamma_EmpzCut",
                  "weight*((E*(1-cos(theta))+sigma_h > 32)&&(E*(1-cos(theta))+sigma_h < 40))"
                  "*(y_ele>0.01 && y_ele<0.95 && Q2_ele>10)*(channel==1)");
  
  // Style
  hE_gamma_noCut->SetLineColor(kBlack);
  hE_gamma_EmpzCut->SetLineColor(kRed);
  hE_gamma_noCut->SetLineWidth(2);
  hE_gamma_EmpzCut->SetLineWidth(2);
  
  // === TOP PAD ===
  pad1->cd();
  hE_gamma_noCut->Draw("hist");
  hE_gamma_EmpzCut->Draw("histsame");

  hE_gamma_noCut->GetYaxis()->SetNdivisions(505);
  hE_gamma_noCut->GetYaxis()->SetTitleSize(0.09);
  hE_gamma_noCut->GetYaxis()->SetTitleOffset(0.5);
  hE_gamma_noCut->GetYaxis()->SetLabelSize(0.08);
  
  auto legend = new TLegend(0.65, 0.72, 0.9, 0.9);
  legend->AddEntry(hE_gamma_noCut, "No E-p_{z} cut", "l");
  legend->AddEntry(hE_gamma_EmpzCut, "32 < E-p_{z} < 40 GeV", "l");
  legend->Draw();

  TLatex Text_com;
  Text_com.SetTextAlign(13);
  // Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 140 GeV");
  // Text_com.DrawLatexNDC(.15,.8,"L_{proj} = 1 fb^{-1}");
  // with kinematic cuts
  Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 140 GeV");
  Text_com.DrawLatexNDC(.15,.8,"Q^{2}_{e} > 10 GeV^{2}, 0.01 < y_{e} < 0.95");
  Text_com.DrawLatexNDC(.15,.75,"L_{proj} = 1 fb^{-1}");
  
  TLatex Text_ePIC;
  Text_ePIC.SetTextSize(0.05);
  Text_ePIC.SetTextFont(62);
  Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");

  TLatex Text_date;
  Text_date.SetTextSize(0.035);
  Text_date.SetTextFont(52);
  // Text_date.DrawLatexNDC(.65,.96,"Simu campaign: MM/YYYY");

  // === BOTTOM PAD (RATIO) ===
  pad2->cd();

  TH1F *h_ratio = (TH1F*)hE_gamma_EmpzCut->Clone("h_ratio");
  h_ratio->Divide(hE_gamma_noCut);
  
  h_ratio->SetTitle("");
  // h_ratio->GetYaxis()->SetTitle("Ratio (Cut / No cut)");
  h_ratio->GetYaxis()->SetTitle("Ratio");
  h_ratio->GetYaxis()->SetRangeUser(0,1);
  // h_ratio->SetLineColor(kBlack);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetTitleSize(0.14);
  h_ratio->GetYaxis()->SetTitleOffset(0.34);
  h_ratio->GetYaxis()->SetLabelSize(0.08);

  h_ratio->GetXaxis()->SetTitleSize(0.1);
  h_ratio->GetXaxis()->SetLabelSize(0.09);
  // h_ratio->SetMarkerStyle(20);
  // h_ratio->SetMarkerSize(0.8);
  
  h_ratio->Draw("E1");

  canvas->cd();
  canvas->Update();
}
