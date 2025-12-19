
void RetrieveKinematics(TString outfile="TreeOutput.root",int fileNum=-1){

  // string extension = "";
  // if (fileNum >= 0) extension.append(std::to_string(fileNum));
  // extension.append(".eicrecon.tree.edm4eic.root");

  TChain *chain = new TChain("events");
  chain->Add("Campaign_25.10.0/18x275/merged_18x275_Q2_1*");


  int nentries = chain->GetEntries();
  cout << nentries << " Entries found !" << endl;

  TTreeReader tr(chain);
  // Inclusive Kinematics Truth
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_x(tr, "InclusiveKinematicsTruth.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_y(tr, "InclusiveKinematicsTruth.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsTruth_Q2(tr, "InclusiveKinematicsTruth.Q2");
  
  // Inclusive Kinematics Electron 
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_x(tr, "InclusiveKinematicsElectron.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_y(tr, "InclusiveKinematicsElectron.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsElectron_Q2(tr, "InclusiveKinematicsElectron.Q2");

  // Inclusive Kinematics JB
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_x(tr, "InclusiveKinematicsJB.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_y(tr, "InclusiveKinematicsJB.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsJB_Q2(tr, "InclusiveKinematicsJB.Q2");

  // Inclusive Kinematics DA
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_x(tr, "InclusiveKinematicsDA.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_y(tr, "InclusiveKinematicsDA.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsDA_Q2(tr, "InclusiveKinematicsDA.Q2");

  // Inclusive Kinematics Sigma
  TTreeReaderArray<Float_t> InclusiveKinematicsSigma_x(tr, "InclusiveKinematicsSigma.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsSigma_y(tr, "InclusiveKinematicsSigma.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsSigma_Q2(tr, "InclusiveKinematicsSigma.Q2");

  // Inclusive Kinematics e-Sigma
  TTreeReaderArray<Float_t> InclusiveKinematicsESigma_x(tr, "InclusiveKinematicsESigma.x");
  TTreeReaderArray<Float_t> InclusiveKinematicsESigma_y(tr, "InclusiveKinematicsESigma.y");
  TTreeReaderArray<Float_t> InclusiveKinematicsESigma_Q2(tr, "InclusiveKinematicsESigma.Q2");

  // Setup tree to store kinematic variables
  Double_t x_true, x_e_ecal, x_e_track, x_jb, x_da, x_sig, x_esig;
  Double_t y_true, y_e_ecal, y_e_track, y_jb, y_da, y_sig, y_esig;
  Double_t Q2_true, Q2_e_ecal, Q2_e_track, Q2_jb, Q2_da, Q2_sig, Q2_esig;
  Double_t W_true;
  Double_t E_scat_true, theta_scat_true, E_scat_reco, theta_scat_reco;
  Double_t E_ele_reco, theta_ele_reco, sigmah_reco, pt_had_reco;
  Double_t E_ele_true, theta_ele_true, sigmah_true, pt_had_true;
  
  TFile *varTreeFile = new TFile(outfile.Data(), "RECREATE");
  TTree *varTree = new TTree("event_var", "event_var");
  
  // Create branches
  varTree->Branch("x_true", &x_true, "x_true/D");
  varTree->Branch("x_e_track", &x_e_track, "x_e_track/D");
  varTree->Branch("x_jb", &x_jb, "x_jb/D");
  varTree->Branch("x_da", &x_da, "x_da/D");
  varTree->Branch("x_sig", &x_sig, "x_sig/D");
  varTree->Branch("x_esig", &x_esig, "x_esig/D");
  
  varTree->Branch("y_true", &y_true, "y_true/D");
  varTree->Branch("y_e_track", &y_e_track, "y_e_track/D");
  varTree->Branch("y_jb", &y_jb, "y_jb/D");
  varTree->Branch("y_da", &y_da, "y_da/D");
  varTree->Branch("y_sig", &y_sig, "y_sig/D");
  varTree->Branch("y_esig", &y_esig, "y_esig/D");
  
  varTree->Branch("Q2_true", &Q2_true, "Q2_true/D");
  varTree->Branch("Q2_e_track", &Q2_e_track, "Q2_e_track/D");
  varTree->Branch("Q2_jb", &Q2_jb, "Q2_jb/D");
  varTree->Branch("Q2_da", &Q2_da, "Q2_da/D");
  varTree->Branch("Q2_sig", &Q2_sig, "Q2_sig/D");
  varTree->Branch("Q2_esig", &Q2_esig, "Q2_esig/D");
  
  varTree->Branch("W_true", &W_true, "W_true/D");
  varTree->Branch("E_scat_true", &E_scat_true, "E_scat_true/D");
  varTree->Branch("theta_scat_true", &theta_scat_true, "theta_scat_true/D");
  varTree->Branch("E_ele_true", &E_ele_true, "E_ele_true/D");
  varTree->Branch("theta_ele_true", &theta_ele_true, "theta_ele_true/D");
  varTree->Branch("sigmah_true", &sigmah_true, "sigmah_true/D");
  varTree->Branch("pt_had_true", &pt_had_true, "pt_had_true/D");
  
  varTree->Branch("E_scat_reco", &E_scat_reco, "E_scat_reco/D");
  varTree->Branch("theta_scat_reco", &theta_scat_reco, "theta_scat_reco/D");
  varTree->Branch("E_ele_reco", &E_ele_reco, "E_ele_reco/D");
  varTree->Branch("theta_ele_reco", &theta_ele_reco, "theta_ele_reco/D");
  varTree->Branch("sigmah_reco", &sigmah_reco, "sigmah_reco/D");
  varTree->Branch("pt_had_reco", &pt_had_reco, "pt_had_reco/D");

  // Event loop
  int counter = 0;
  while (tr.Next()){

    if (!InclusiveKinematicsTruth_y.GetSize()) continue;
    if (!InclusiveKinematicsElectron_y.GetSize()) continue;
    if (!InclusiveKinematicsJB_y.GetSize()) continue;

    // if (counter > 500) break;
    counter++;
    if (counter % 10000 == 0){
      cout << "Event: " << counter << endl;
    }

    x_true = InclusiveKinematicsTruth_x[0];
    y_true = InclusiveKinematicsTruth_y[0];
    Q2_true = InclusiveKinematicsTruth_Q2[0];

    x_e_track = InclusiveKinematicsElectron_x[0];
    y_e_track = InclusiveKinematicsElectron_y[0];
    Q2_e_track = InclusiveKinematicsElectron_Q2[0];

    x_jb = InclusiveKinematicsJB_x[0];
    y_jb = InclusiveKinematicsJB_y[0];
    Q2_jb = InclusiveKinematicsJB_Q2[0];

    x_da = InclusiveKinematicsDA_x[0];
    y_da = InclusiveKinematicsDA_y[0];
    Q2_da = InclusiveKinematicsDA_Q2[0];

    x_sig = InclusiveKinematicsSigma_x[0];
    y_sig = InclusiveKinematicsSigma_y[0];
    Q2_sig = InclusiveKinematicsSigma_Q2[0];

    x_esig = InclusiveKinematicsESigma_x[0];
    y_esig = InclusiveKinematicsESigma_y[0];
    Q2_esig = InclusiveKinematicsESigma_Q2[0];

    E_scat_reco = 0;
    E_ele_reco = 0;
    theta_scat_reco = 0;
    theta_ele_reco = 0;
    sigmah_reco = 2*18.*y_jb;
    pt_had_reco = TMath::Sqrt(Q2_jb*(1-y_jb));

    E_scat_true = 0;
    E_ele_true = 0;
    theta_scat_true = 0;
    theta_ele_true = 0;
    // TODO: update to use true beam energy for event (or better, the actual particles)
    sigmah_true = y_true*2*18.;
    pt_had_true = TMath::Sqrt(Q2_true*(1-y_true));

    // Fill the tree for this event
    varTree->Fill();
  }// End event loop

  // Write the tree to the file and close it
  varTree->Write();
  varTreeFile->Close();
}
