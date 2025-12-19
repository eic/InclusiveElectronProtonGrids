#ifndef TREE_WRITER_H
#define TREE_WRITER_H

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// ROOT globals
static TFile *tw_outfile = nullptr;
static TTree *tw_tree    = nullptr;

// Branch variables
static Double_t tw_x_true, tw_x_ele, tw_x_jb, tw_x_da, tw_x_sig, tw_x_esig;
static Double_t tw_y_true, tw_y_ele, tw_y_jb, tw_y_da, tw_y_sig, tw_y_esig;
static Double_t tw_Q2_true, tw_Q2_ele, tw_Q2_jb, tw_Q2_da, tw_Q2_sig, tw_Q2_esig;

static Double_t tw_E, tw_theta;
static Double_t tw_sigma_h, tw_pt_had;
static Double_t tw_sigma_h_true, tw_pt_had_true;

static Double_t tw_weight;

// Initialise the output file & create the tree
inline void treewriter_init(const char *filename = "event_variables.root")
{
    tw_outfile = new TFile(filename, "RECREATE");
    tw_tree    = new TTree("event_var", "event_var");

    // Kinematics
    tw_tree->Branch("x_true", &tw_x_true, "x_true/D");
    tw_tree->Branch("x_ele",  &tw_x_ele,  "x_ele/D");
    tw_tree->Branch("x_jb",   &tw_x_jb,   "x_jb/D");
    tw_tree->Branch("x_da",   &tw_x_da,   "x_da/D");
    tw_tree->Branch("x_sig",  &tw_x_sig,  "x_sig/D");
    tw_tree->Branch("x_esig", &tw_x_esig, "x_esig/D");

    tw_tree->Branch("y_true", &tw_y_true, "y_true/D");
    tw_tree->Branch("y_ele",  &tw_y_ele,  "y_ele/D");
    tw_tree->Branch("y_jb",   &tw_y_jb,   "y_jb/D");
    tw_tree->Branch("y_da",   &tw_y_da,   "y_da/D");
    tw_tree->Branch("y_sig",  &tw_y_sig,  "y_sig/D");
    tw_tree->Branch("y_esig", &tw_y_esig, "y_esig/D");

    tw_tree->Branch("Q2_true", &tw_Q2_true, "Q2_true/D");
    tw_tree->Branch("Q2_ele",  &tw_Q2_ele,  "Q2_ele/D");
    tw_tree->Branch("Q2_jb",   &tw_Q2_jb,   "Q2_jb/D");
    tw_tree->Branch("Q2_da",   &tw_Q2_da,   "Q2_da/D");
    tw_tree->Branch("Q2_sig",  &tw_Q2_sig,  "Q2_sig/D");
    tw_tree->Branch("Q2_esig", &tw_Q2_esig, "Q2_esig/D");

    // Inclusive FS info
    tw_tree->Branch("E",            &tw_E,            "E/D");
    tw_tree->Branch("theta",        &tw_theta,        "theta/D");
    tw_tree->Branch("sigma_h",      &tw_sigma_h,      "sigma_h/D");
    tw_tree->Branch("pt_had",       &tw_pt_had,       "pt_had/D");
    tw_tree->Branch("sigma_h_true", &tw_sigma_h_true, "sigma_h_true/D");
    tw_tree->Branch("pt_had_true",  &tw_pt_had_true,  "pt_had_true/D");

    // Event weight
    tw_tree->Branch("weight", &tw_weight, "weight/D");
}

// Fill tree
inline void treewriter_fill(
    Float_t x_truth, Float_t x_ele, Float_t x_jb, Float_t x_da, Float_t x_sig, Float_t x_esig,
    Float_t y_truth, Float_t y_ele, Float_t y_jb, Float_t y_da, Float_t y_sig, Float_t y_esig,
    Float_t Q2_truth, Float_t Q2_ele, Float_t Q2_jb, Float_t Q2_da, Float_t Q2_sig, Float_t Q2_esig,
    Float_t E, Float_t theta, Float_t sigma_h, Float_t pt_had,
    Float_t sigma_h_true, Float_t pt_had_true,
    Float_t weight )
{
    // Assign values
    tw_x_true = x_truth;    tw_x_ele = x_ele;    tw_x_jb = x_jb;
    tw_x_da   = x_da;       tw_x_sig = x_sig;    tw_x_esig = x_esig;

    tw_y_true = y_truth;    tw_y_ele = y_ele;    tw_y_jb = y_jb;
    tw_y_da   = y_da;       tw_y_sig = y_sig;    tw_y_esig = y_esig;

    tw_Q2_true = Q2_truth;  tw_Q2_ele = Q2_ele;  tw_Q2_jb = Q2_jb;
    tw_Q2_da   = Q2_da;     tw_Q2_sig = Q2_sig;  tw_Q2_esig = Q2_esig;

    tw_E = E;
    tw_theta = theta;
    tw_sigma_h = sigma_h;
    tw_pt_had  = pt_had;
    tw_sigma_h_true = sigma_h_true;
    tw_pt_had_true  = pt_had_true;

    tw_weight = weight;

    tw_tree->Fill();
}

// Save file
inline void treewriter_close()
{
    if (tw_outfile) {
        tw_outfile->cd();
        tw_tree->Write();
        tw_outfile->Close();
        delete tw_outfile;
        tw_outfile = nullptr;
    }
}

#endif
