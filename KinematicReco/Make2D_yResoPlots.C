int Make2D_yResoPlots(bool useLogScale = true){

  // Settings 
  bool usePercent = true;
  double beam_s = 4*18*275;

  int nPointsX {20};
  int nPointsQ2 {20};



  char *histname = new char[30];
  TH1D *y_reso_e_track[nPointsX][nPointsQ2];
  TH1D *y_reso_jb[nPointsX][nPointsQ2];
  TH1D *y_reso_da[nPointsX][nPointsQ2];
  TH1D *y_reso_sig[nPointsX][nPointsQ2];
  TH1D *y_reso_esig[nPointsX][nPointsQ2];

  for (int i{0}; i<nPointsX; i++){
    for (int j{0}; j<nPointsQ2; j++){

      sprintf(histname, "y_reso_e_track%d_%d", i, j);
      y_reso_e_track[i][j] = new TH1D(histname,
				      ";(y_{rec}-y_{true})/y_{true};",
				      2000, -1, 1);

      sprintf(histname, "y_reso_jb%d_%d", i, j);
      y_reso_jb[i][j] = new TH1D(histname,
				      ";(y_{rec}-y_{true})/y_{true};",
				      2000, -1, 1);

      sprintf(histname, "y_reso_da%d_%d", i, j);
      y_reso_da[i][j] = new TH1D(histname,
				      ";(y_{rec}-y_{true})/y_{true};",
				      2000, -1, 1);

      sprintf(histname, "y_reso_sig%d_%d", i, j);
      y_reso_sig[i][j] = new TH1D(histname,
				      ";(y_{rec}-y_{true})/y_{true};",
				      2000, -1, 1);

      sprintf(histname, "y_reso_esig%d_%d", i, j);
      y_reso_esig[i][j] = new TH1D(histname,
				      ";(y_{rec}-y_{true})/y_{true};",
				      2000, -1, 1);

    }
  }

  TChain *tree = new TChain("event_var");
  tree->Add("./TreeOutput.root");
  int nentries = tree->GetEntries();

  Double_t x_true, x_e_ecal, x_e_track, x_jb, x_da, x_sig, x_esig;
  Double_t y_true, y_e_ecal, y_e_track, y_jb, y_da, y_sig, y_esig;
  Double_t Q2_true, Q2_e_ecal, Q2_e_track, Q2_jb, Q2_da, Q2_sig, Q2_esig;

  tree->SetBranchAddress("x_true", &x_true);
  tree->SetBranchAddress("y_true", &y_true);
  tree->SetBranchAddress("Q2_true", &Q2_true);

  tree->SetBranchAddress("x_e_track", &x_e_track);
  tree->SetBranchAddress("y_e_track", &y_e_track);
  tree->SetBranchAddress("Q2_e_track", &Q2_e_track);

  tree->SetBranchAddress("x_jb", &x_jb);
  tree->SetBranchAddress("y_jb", &y_jb);
  tree->SetBranchAddress("Q2_jb", &Q2_jb);

  tree->SetBranchAddress("x_da", &x_da);
  tree->SetBranchAddress("y_da", &y_da);
  tree->SetBranchAddress("Q2_da", &Q2_da);

  tree->SetBranchAddress("x_sig", &x_sig);
  tree->SetBranchAddress("y_sig", &y_sig);
  tree->SetBranchAddress("Q2_sig", &Q2_sig);

  tree->SetBranchAddress("x_esig", &x_esig);
  tree->SetBranchAddress("y_esig", &y_esig);
  tree->SetBranchAddress("Q2_esig", &Q2_esig);


  Double_t x_bins[nPointsX+1];
  Double_t y_bins[nPointsQ2+1];

  Double_t xmin_exponent = -4;// bin from 1e-4 to 1e0 in base 10
  Double_t xmax_exponent = 0;

  Double_t ymin_exponent = 0;// bin from 1 to 1e3 in base 10
  Double_t ymax_exponent = 4;

  Double_t xmin = 0;
  Double_t xmax = 1;

  Double_t ymin = 1;
  Double_t ymax = 1000;

  if (useLogScale){
    for (int i {0}; i <= nPointsX; i++){
      x_bins[i] = TMath::Power(10, xmin_exponent + i*(xmax_exponent - xmin_exponent)/nPointsX);
    }
  }
  else {
    for (int i {0}; i <= nPointsX; i++){
      x_bins[i] = xmin + i*(xmax - xmin)/nPointsX;
    }
  }
  if (useLogScale){
    for (int i {0}; i <= nPointsQ2; i++){
      y_bins[i] = TMath::Power(10, ymin_exponent + i*(ymax_exponent - ymin_exponent)/nPointsQ2);
    }
  }
  else {
    for (int i {0}; i <= nPointsQ2; i++){
      y_bins[i] = ymin + i*(ymax - ymin)/nPointsQ2;
    }
  }

  for (auto item : x_bins){cout << item << endl;}
  for (auto item : y_bins){cout << item << endl;}

  int histogramIndexX, histogramIndexQ2;
  for (int i {0}; i<nentries; i++){// loop over events
    tree->GetEntry(i);

    if (Q2_true > 1 && Q2_true < 10000 && y_true < 1.00 && y_true > 0.001){
      for (int j{0}; j<nPointsX;j++){
	if (x_true > x_bins[j] && x_true < x_bins[j+1]){
	  histogramIndexX = j;
	  for (int k{0}; k<nPointsQ2;k++){
	    if (Q2_true > y_bins[k] && Q2_true < y_bins[k+1]){
	      histogramIndexQ2 = k;

	      y_reso_e_track[histogramIndexX][histogramIndexQ2]->Fill((y_e_track-y_true)/y_true);
	      y_reso_jb[histogramIndexX][histogramIndexQ2]->Fill((y_jb-y_true)/y_true);
	      y_reso_da[histogramIndexX][histogramIndexQ2]->Fill((y_da-y_true)/y_true);
	      y_reso_sig[histogramIndexX][histogramIndexQ2]->Fill((y_sig-y_true)/y_true);
	      y_reso_esig[histogramIndexX][histogramIndexQ2]->Fill((y_esig-y_true)/y_true);
	    }
	  }
	}
      }
    }
  }


  TH2D *y_reso_2D_e_track = new TH2D("y_reso_2D_e_track", ";x;Q^{2} [GeV^{2}]", nPointsX, x_bins, nPointsQ2, y_bins);
  TH2D *y_reso_2D_jb = new TH2D("y_reso_2D_jb", ";x;Q^{2} [GeV^{2}]", nPointsX, x_bins, nPointsQ2, y_bins);
  TH2D *y_reso_2D_da = new TH2D("y_reso_2D_da", ";x;Q^{2} [GeV^{2}]", nPointsX, x_bins, nPointsQ2, y_bins);
  TH2D *y_reso_2D_sig = new TH2D("y_reso_2D_sig", ";x;Q^{2} [GeV^{2}]", nPointsX, x_bins, nPointsQ2, y_bins);
  TH2D *y_reso_2D_esig = new TH2D("y_reso_2D_esig", ";x;Q^{2} [GeV^{2}]", nPointsX, x_bins, nPointsQ2, y_bins);

  if (usePercent){
    for (int i{0}; i<nPointsX; i++){
      for (int j{0}; j<nPointsQ2; j++){
	if (y_reso_e_track[i][j]->GetEntries() > 10 && abs(y_reso_e_track[i][j]->GetMean()) < 0.5){
	  y_reso_2D_e_track->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, 100*y_reso_e_track[i][j]->GetRMS()); 
	}
	if (y_reso_jb[i][j]->GetEntries() > 10 && abs(y_reso_jb[i][j]->GetMean()) < 0.5){
	  y_reso_2D_jb->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, 100*y_reso_jb[i][j]->GetRMS()); 
	}
	if (y_reso_da[i][j]->GetEntries() > 10 && abs(y_reso_da[i][j]->GetMean()) < 0.5){
	  y_reso_2D_da->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, 100*y_reso_da[i][j]->GetRMS()); 
	}
	if (y_reso_sig[i][j]->GetEntries() > 10 && abs(y_reso_sig[i][j]->GetMean()) < 0.5){
	  y_reso_2D_sig->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, 100*y_reso_sig[i][j]->GetRMS()); 
	}
	if (y_reso_esig[i][j]->GetEntries() > 10 && abs(y_reso_esig[i][j]->GetMean()) < 0.5){
	  y_reso_2D_esig->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, 100*y_reso_esig[i][j]->GetRMS()); 
	}
      }
    }
  }
  else{
    for (int i{0}; i<nPointsX; i++){
      for (int j{0}; j<nPointsQ2; j++){
	if (y_reso_e_track[i][j]->GetEntries() > 10){
	  y_reso_2D_e_track->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, y_reso_e_track[i][j]->GetRMS()); 
	}
	if (y_reso_jb[i][j]->GetEntries() > 10){
	  y_reso_2D_jb->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, y_reso_jb[i][j]->GetRMS()); 
	}
	if (y_reso_da[i][j]->GetEntries() > 10){
	  y_reso_2D_da->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, y_reso_da[i][j]->GetRMS()); 
	}
	if (y_reso_sig[i][j]->GetEntries() > 10){
	  y_reso_2D_sig->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, y_reso_sig[i][j]->GetRMS()); 
	}
	if (y_reso_esig[i][j]->GetEntries() > 10){
	  y_reso_2D_esig->Fill((x_bins[i]+x_bins[i+1])/2, (y_bins[j]+y_bins[j+1])/2, y_reso_esig[i][j]->GetRMS()); 
	}
      }
    }
  }


  TFile *myHistFile = new TFile("./HistogramsLog.root", "RECREATE");

  myHistFile->mkdir("xQ2_hists");
  myHistFile->cd("xQ2_hists");
  y_reso_2D_e_track->Write();
  y_reso_2D_jb->Write();
  y_reso_2D_da->Write();
  y_reso_2D_sig->Write();
  y_reso_2D_esig->Write();
  myHistFile->cd("../");
  myHistFile->mkdir("eMethod");
  myHistFile->mkdir("JBMethod");
  myHistFile->mkdir("DAMethod");
  myHistFile->mkdir("SigmaMethod");
  myHistFile->mkdir("eSigmaMethod");

  for (int i{0}; i<nPointsX; i++){
    for (int j{0}; j<nPointsQ2; j++){
      myHistFile->cd("eMethod");
      y_reso_e_track[i][j]->Write();
      myHistFile->cd("../");
      myHistFile->cd("JBMethod");
      y_reso_jb[i][j]->Write();
      myHistFile->cd("../");
      myHistFile->cd("DAMethod");
      y_reso_da[i][j]->Write();
      myHistFile->cd("../");
      myHistFile->cd("SigmaMethod");
      y_reso_sig[i][j]->Write();
      myHistFile->cd("../");
      myHistFile->cd("eSigmaMethod");
      y_reso_esig[i][j]->Write();
      myHistFile->cd("../");
    }
  }
  myHistFile->Close();

  std::ofstream resolutions_18_275_outfile;
  resolutions_18_275_outfile.open("18_275_resolutions.csv");
  resolutions_18_275_outfile << "x," << "Q2," << "y," << "SE," << "JB," << "DA," << "SIG," << "ESIG\n";

  if (useLogScale){
    for (int binX{1}; binX<=y_reso_2D_jb->GetNbinsX(); binX++){
      for (int binQ2{1}; binQ2<=y_reso_2D_jb->GetNbinsY(); binQ2++){
	resolutions_18_275_outfile << y_reso_2D_e_track->GetXaxis()->GetBinCenterLog(binX) << ","
				   << y_reso_2D_e_track->GetYaxis()->GetBinCenterLog(binQ2) << ","
				   << y_reso_2D_e_track->GetYaxis()->GetBinCenterLog(binQ2)/(beam_s*y_reso_2D_e_track->GetXaxis()->GetBinCenterLog(binX)) << ","
				   << y_reso_2D_e_track->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_jb->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_da->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_sig->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_esig->GetBinContent(binX, binQ2) << "\n";
      }
    }

  }
  else{
    for (int binX{1}; binX<=y_reso_2D_jb->GetNbinsX(); binX++){
      for (int binQ2{1}; binQ2<=y_reso_2D_jb->GetNbinsY(); binQ2++){
	resolutions_18_275_outfile << y_reso_2D_e_track->GetXaxis()->GetBinCenter(binX) << ","
				   << y_reso_2D_e_track->GetYaxis()->GetBinCenter(binQ2) << ","
				   << y_reso_2D_e_track->GetYaxis()->GetBinCenter(binQ2)/(beam_s*y_reso_2D_e_track->GetXaxis()->GetBinCenter(binX)) << ","
				   << y_reso_2D_e_track->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_jb->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_da->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_sig->GetBinContent(binX, binQ2) << ","
				   << y_reso_2D_esig->GetBinContent(binX, binQ2) << "\n";
      }
    }
  }

  return 0;
}
