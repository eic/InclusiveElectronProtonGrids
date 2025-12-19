#include "Boost.h"
#include "Beam.h"

// PODIO
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

// DATA MODEL
#include "edm4eic/InclusiveKinematicsCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoClusterParticleAssociationCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/HadronicFinalStateCollection.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/utils/vector_utils.h"

// GENERAL
#include <filesystem>
#include <vector>

using ROOT::Math::PxPyPzEVector;

std::vector<float> calc_elec_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_jb_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_da_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_sig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);
std::vector<float> calc_esig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam);

TLorentzVector TransformLabToHeadOnFrame(TLorentzVector eBeam, TLorentzVector pBeam, TLorentzVector Lvec);
std::vector<std::string> get_filevector(std::string filename);
double getPolarAngle(const edm4hep::MCParticle& mcParticle);

void Write_RadCorr_Tree() {
  // Settings
  Float_t E_ebeam = 18;
  Float_t E_pbeam = 275;
  Float_t target_lumi = 1.;
  Float_t m_e = 0.000511;
  Float_t m_p = 0.938;
  double xAngle = 25e-3;
  PxPyPzEVector pni, ei;
  ei.SetPxPyPzE(0, 0, -E_ebeam, sqrt(E_ebeam*E_ebeam+m_e*m_e));
  pni.SetPxPyPzE(-1*E_pbeam*TMath::Sin(xAngle), 0, E_pbeam*TMath::Cos(xAngle), sqrt(E_pbeam*E_pbeam+m_p*m_p));
  
  // paths to file lists
  std::vector<TString> filenames = {
    "../RECO/djangoh_Rad_ep_18x275_Q2_1_10/filelist.txt",
    "../RECO/djangoh_Rad_ep_18x275_Q2_10_100/filelist.txt",
    "../RECO/djangoh_Rad_ep_18x275_Q2_100_1000/filelist.txt",
    "../RECO/djangoh_Rad_ep_18x275_Q2_1000_10000/filelist.txt"
  };
  
  std::vector<double> gen_xsecs = {0.6733778912,0.07405000239,0.3760005602e-2,0.8615325533e-4}; // ub
  
  const int nBinsX = 21;
  const int nBinsQ2 = 21;
  double xBins[nBinsX + 1];
  double Q2Bins[nBinsQ2 + 1];
  for (int i = 0; i <= nBinsX; ++i) {
    xBins[i] = TMath::Power(10, -4 + i * 0.2); // Logarithmic bins for x
  }
  for (int i = 0; i <= nBinsQ2; ++i) {
    Q2Bins[i] = TMath::Power(10, i * 0.2); // Logarithmic bins for Q2
  }
  
  Float_t x_born_truth, x_rad_truth, x_ele, x_jb, x_da, x_sig, x_esig;
  Float_t y_born_truth, y_rad_truth, y_ele, y_jb, y_da, y_sig, y_esig;
  Float_t Q2_born_truth, Q2_rad_truth, Q2_ele, Q2_jb, Q2_da, Q2_sig, Q2_esig;
  Float_t W_truth, nu_truth;
  Float_t E, theta, sigma_h, pt_had;
  Float_t sigma_h_true, pt_had_true;
  Float_t weight;
  Int_t channel; // 0=nonrad, 1=ISR, 2=FSR
  Float_t E_gamma, theta_gamma;
  
  
  int file_index = 0;

  TFile *file = new TFile("new_track_ele_rad_output_tree.root", "RECREATE");
  TTree *tree = new TTree("events", "QED Rad Variables Tree");

  // Set branches
  tree->Branch("x_born_truth", &x_born_truth);
  tree->Branch("x_rad_truth", &x_rad_truth);
  tree->Branch("x_ele", &x_ele);
  tree->Branch("x_jb", &x_jb);
  tree->Branch("x_da", &x_da);
  tree->Branch("x_sig", &x_sig);
  tree->Branch("x_esig", &x_esig);
  
  tree->Branch("y_born_truth", &y_born_truth);
  tree->Branch("y_rad_truth", &y_rad_truth);
  tree->Branch("y_ele", &y_ele);
  tree->Branch("y_jb", &y_jb);
  tree->Branch("y_da", &y_da);
  tree->Branch("y_sig", &y_sig);
  tree->Branch("y_esig", &y_esig);
  
  tree->Branch("Q2_born_truth", &Q2_born_truth);
  tree->Branch("Q2_rad_truth", &Q2_rad_truth);
  tree->Branch("Q2_ele", &Q2_ele);
  tree->Branch("Q2_jb", &Q2_jb);
  tree->Branch("Q2_da", &Q2_da);
  tree->Branch("Q2_sig", &Q2_sig);
  tree->Branch("Q2_esig", &Q2_esig);

  tree->Branch("W_truth", &W_truth);
  tree->Branch("nu_truth", &nu_truth);
  
  tree->Branch("E", &E);
  tree->Branch("theta", &theta);
  tree->Branch("sigma_h", &sigma_h);
  tree->Branch("pt_had", &pt_had);
  tree->Branch("sigma_h_true", &sigma_h_true);
  tree->Branch("pt_had_true", &pt_had_true);
  
  tree->Branch("weight", &weight);
  tree->Branch("channel", &channel); // Int_t
  tree->Branch("E_gamma", &E_gamma);
  tree->Branch("theta_gamma", &theta_gamma);
  
  for (auto filename : filenames){
    auto inFiles = get_filevector(filenames[file_index].Data());
    
    auto reader = podio::ROOTReader();
    reader.openFiles(inFiles);
    double nEntries = reader.getEntries("events");
    cout << nEntries << " events found" << endl;
    
    weight = target_lumi/((nEntries)/(gen_xsecs[file_index]*1e9));
    cout << "Weight for file set " << file_index << ": " << weight << endl;
    
    for (size_t i = 0; i < reader.getEntries("events"); i++) {// begin event loop
      const auto event = podio::Frame(reader.readNextEntry("events"));
      if (i%100==0) cout << i << " events processed" << endl;
      
      // Retrieve Inclusive Kinematics Collections
      auto& kin_truth = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
      auto& kin_electron = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
      auto& kin_jb = event.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsJB");
      
      // Retrieve Scattered electron and HFS
      auto& eleCollection = event.get<edm4eic::ReconstructedParticleCollection>("ScatteredElectronsTruth");
      // auto& hfsCollection = event.get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");
      auto& ecalClusters = event.get<edm4eic::ClusterCollection>("EcalClusters");
      
      auto& rcparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
      auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticlesHeadOnFrameNoBeamFX");
      
      // Store kinematics from InclusiveKinematics branches
      if (kin_truth.empty() || kin_electron.empty() || kin_jb.empty()) continue;
      
      x_rad_truth = kin_truth.x()[0];
      y_rad_truth = kin_truth.y()[0];
      Q2_rad_truth = kin_truth.Q2()[0];
      
      auto boost = eicrecon::determine_boost(ei, pni);
      
      PxPyPzEVector scat_ele;
      E = eleCollection[0].getEnergy();
      // if (!eleCollection[0].getClusters().empty()){
      // E = eleCollection[0].getClusters()[0].getEnergy();
      // }
      auto& ele_momentum = eleCollection[0].getMomentum();
      scat_ele.SetPxPyPzE(ele_momentum.x, ele_momentum.y, ele_momentum.z, E);
      theta = scat_ele.Theta();
      
      // Comment these out since we're getting the HFS manually
      // sigma_h = hfsCollection[0].getSigma();
      // pt_had = hfsCollection[0].getPT();

      // Reset variables
      E_gamma=0;theta_gamma=0;channel=0;
    
      PxPyPzEVector photon_fv;
      for (const auto mcp : mcparts) {
	if (mcp.getPDG() == 22 && mcp.getGeneratorStatus() == 1){
	  // cout << mcp.getParents()[0].getObjectID().index << endl;
	  int index = mcp.getParents()[0].getObjectID().index;
	  if (index == 0){
	    channel = 1;
	    E_gamma = mcp.getEnergy();
	    theta_gamma = getPolarAngle(mcp);
	    auto mcp_p           = mcp.getMomentum();
	    auto mcp_p_mag       = edm4hep::utils::magnitude(mcp_p);
	    photon_fv.SetPxPyPzE(mcp_p.x,mcp_p.y,mcp_p.z,mcp_p_mag);
	    break;
	  }
	  if (index == 3){
	    channel = 2;
	    E_gamma = mcp.getEnergy();
	    theta_gamma = getPolarAngle(mcp);
	    auto mcp_p           = mcp.getMomentum();
	    auto mcp_p_mag       = edm4hep::utils::magnitude(mcp_p);
	    photon_fv.SetPxPyPzE(mcp_p.x,mcp_p.y,mcp_p.z,mcp_p_mag);
	    break;
	  }
	}
      }
      // cout << channel << " " << E_gamma << " " << theta_gamma << endl;

      // Get incoming electron beam
      const auto ei_coll = eicrecon::find_first_beam_electron(&mcparts);
      if (ei_coll.empty()) {
	cout << "No beam electron found" << endl;
	return;
      }
      const auto ei_p           = ei_coll[0].getMomentum();
      const auto ei_p_mag       = edm4hep::utils::magnitude(ei_p);
      static const auto ei_mass = ei_coll[0].getMass();
      PxPyPzEVector ei(ei_p.x, ei_p.y, ei_p.z, std::hypot(ei_p_mag, ei_mass));
      
      // Get incoming hadron beam
      const auto pi_coll = eicrecon::find_first_beam_hadron(&mcparts);
      if (pi_coll.empty()) {
	cout << "No beam hadron found" << endl;
	return;
      }
      const auto pi_p     = pi_coll[0].getMomentum();
      const auto pi_p_mag = edm4hep::utils::magnitude(pi_p);
      const auto pi_mass  = pi_coll[0].getMass();
      PxPyPzEVector pi(pi_p.x, pi_p.y, pi_p.z, std::hypot(pi_p_mag, pi_mass));
      
      // Get first scattered electron
      // Scattered electron. Currently taken as first status==1 electron in HEPMC record,
      // which seems to be correct based on a cursory glance at the Pythia8 output. In the future,
      // it may be better to trace back each final-state electron and see which one originates from
      // the beam.
      const auto ef_coll = eicrecon::find_first_scattered_electron(&mcparts);
      if (ef_coll.empty()) {
	cout << "No truth scattered electron found" << endl;
	return;
      }
      const auto ef_p           = ef_coll[0].getMomentum();
      const auto ef_p_mag       = edm4hep::utils::magnitude(ef_p);
      static const auto ef_mass = ef_coll[0].getMass();
      PxPyPzEVector ef(ef_p.x, ef_p.y, ef_p.z, std::hypot(ef_p_mag, ef_mass));

      // correct the electron four vectors
      // TODO

      // DIS kinematics calculations
      if (channel == 1) ei -= photon_fv;
      if (channel == 2) ef += photon_fv;
      const auto q        = ei - ef;
      const auto q_dot_pi = q.Dot(pi);
      Q2_born_truth       = -q.Dot(q);
      y_born_truth        = q_dot_pi / ei.Dot(pi);
      nu_truth      = q_dot_pi / pi_mass;
      x_born_truth        = Q2_born_truth / (2. * q_dot_pi);
      W_truth        = sqrt(pi_mass * pi_mass + 2. * q_dot_pi - Q2_born_truth);
      
      double pxsum = 0;
      double pysum = 0;
      double pzsum = 0;
      double Esum  = 0;
      const auto ef_rc_id{eleCollection[0].getObjectID().index};
      for (const auto p : rcparts) {
	bool isHadron = true;
	// Check if it's the scattered electron
      if (p.getObjectID().index != ef_rc_id) {
	// Lorentz vector in lab frame
	PxPyPzEVector hf_lab(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
	// Boost to colinear frame
	PxPyPzEVector hf_boosted = boost(hf_lab);
	// PxPyPzEVector hf_boosted = hf_lab;

	pxsum += hf_boosted.Px();
	pysum += hf_boosted.Py();
	pzsum += hf_boosted.Pz();
	Esum += hf_boosted.E();
      }
    }
    sigma_h = Esum - pzsum;
    pt_had = sqrt(pxsum*pxsum + pysum*pysum);
    sigma_h_true = y_born_truth*2*E_ebeam;
    pt_had_true = sqrt(Q2_born_truth*(1 - y_born_truth));
   
    // Calculate kinematics manually
    std::vector<float> elec_reco = calc_elec_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> jb_reco = calc_jb_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> da_reco = calc_da_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> sigma_reco = calc_sig_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);
    std::vector<float> esigma_reco = calc_esig_method(E, theta, pt_had, sigma_h, E_ebeam, E_pbeam);

    x_ele = elec_reco[0];
    x_jb = jb_reco[0];
    x_da = da_reco[0];
    x_sig = sigma_reco[0];
    x_esig = esigma_reco[0];
   
    y_ele = elec_reco[1];
    y_jb = jb_reco[1];
    y_da = da_reco[1];
    y_sig = sigma_reco[1];
    y_esig = esigma_reco[1];
   
    Q2_ele = elec_reco[2];
    Q2_jb = jb_reco[2];
    Q2_da = da_reco[2];
    Q2_sig = sigma_reco[2];
    Q2_esig = esigma_reco[2];
   
    tree->Fill();
   
  }// end event loop
  file_index++;
  }
  // Writing the tree
  tree->Write();
  file->Close();
 
 
  cout << "Done!" << endl;
}

// electron method
std::vector<float> calc_elec_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
  float y = 1. - (E/E_ebeam)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// jb method
std::vector<float> calc_jb_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float y = sigma_h/(2*E_ebeam);
  float Q2 = pt_had*pt_had / (1-y);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// float angle method
std::vector<float> calc_da_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float alpha_h = sigma_h/pt_had;
  float alpha_e = TMath::Tan(theta/2);
  float y = alpha_h / (alpha_e + alpha_h);
  float Q2 = 4*E_ebeam*E_ebeam / (alpha_e * (alpha_h + alpha_e));
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// sigma method
std::vector<float> calc_sig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float y = sigma_h/(sigma_h + E*(1 - TMath::Cos(theta))); 
  float Q2 = E*E*TMath::Sin(theta)*TMath::Sin(theta) / (1-y);
  float x = Q2/(4*E_ebeam*E_pbeam*y);
  return {x, y, Q2};
}

// e-sigma method
std::vector<float> calc_esig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
  float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
  float x = calc_sig_method(E,theta,pt_had,sigma_h,E_ebeam,E_pbeam)[0];
  float y = Q2/(4*E_ebeam*E_pbeam*x);
  return {x, y, Q2};
}

TLorentzVector TransformLabToHeadOnFrame(TLorentzVector eBeam, TLorentzVector pBeam, TLorentzVector Lvec) {
  TLorentzVector cmBoost = eBeam + pBeam;
  // Define boost for back to back beams
  TLorentzVector boost(-cmBoost[0],-cmBoost[1],-cmBoost[2],cmBoost[3]);
  TVector3 b = boost.BoostVector();
 
  // Define boost to restore original energies
  // TLorentzVector boostBack(0.0,0.0,cmBoost[2],cmBoost[3]);
  TLorentzVector boostBack(0.0,0.0,sqrt(cmBoost[0]*cmBoost[0]+cmBoost[1]*cmBoost[1]+cmBoost[2]*cmBoost[2]),cmBoost[3]);
  TVector3 bb = boostBack.BoostVector();
 
  // Define angles for rotations
  TLorentzVector p = pBeam;
  p.Boost(b);
  double rotAboutY = -1.0*TMath::ATan2(p.Px(),p.Pz());
  double rotAboutX = +1.0*TMath::ATan2(p.Py(),p.Pz());
 
  Lvec.Boost(b);
  Lvec.RotateY(rotAboutY);
  Lvec.RotateX(rotAboutX);
  Lvec.Boost(bb);
  return Lvec;
}

std::vector<std::string> get_filevector(std::string filename){
  std::vector<std::string> fileVector;
  std::ifstream in(filename);
  std::string file("");
  while (in >> file) fileVector.push_back(file.data());
  return fileVector;
}

double getPolarAngle(const edm4hep::MCParticle& mcParticle) {
    double px = mcParticle.getMomentum().x;
    double py = mcParticle.getMomentum().y;
    double pz = mcParticle.getMomentum().z;

    double pMag = std::sqrt(px*px + py*py + pz*pz);

    // Guard against division by zero
    if (pMag == 0.0) return 0.0;

    double cosTheta = pz / pMag;
    
    // Clamp value to [-1,1] to avoid domain error in acos
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;

    return std::acos(cosTheta);  // returns theta in radians
}
