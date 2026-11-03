#pragma once

#include "common.hxx"

// Namespaces

using namespace std;
using namespace RooFit;

/*
/////////////////////////////////////////////////////////////
INPUT FILES
/////////////////////////////////////////////////////////////
*/

// Data file
TFile * f_IN;

// File locations
TString f_BHT1 = "../../data/BHT1.masssp.root";
TString f_BHT2 = "../../data/BHT2_251128.unfolded_trees.root";

/*
/////////////////////////////////////////////////////////////
NAME HANDLERS
/////////////////////////////////////////////////////////////
*/

TString s_tag = "default";	// tag for output files

// Curent date (YYMMDD) for output purposes
time_t timestamp = time(NULL);
struct tm datetime = *localtime(&timestamp);

char time_output[7];

strftime(time_output, 7, "%y%m%d", &datetime);

TString s_date = time_output;

/*
/////////////////////////////////////////////////////////////
RUN PARAMETERS
/////////////////////////////////////////////////////////////
*/

Double_t kMASS_MIN = 7.6;	// x-axis lower limit for fitting
Double_t kMASS_MAX = 13.0;	// x-axis upper limit for fitting

Double_t kBIN_WIDTH = 0.2;	// Bin width for drawing

Int_t kSIGMA = 1;	// sigma for error drawing

Bool_t kLINEAR_ERR = true;	// Use linear error visualisation (faster, less accurate for strong correlations)

/*
/////////////////////////////////////////////////////////////
GLOBAL VARIABLES
/////////////////////////////////////////////////////////////
*/

// Global maps to store parameters for different bins/triggers
map<TString, vector<vector<Double_t>>> map_par_sig;
map<TString, vector<vector<Double_t>>> map_par_bg;
map<TString, vector<vector<Double_t>>> map_frac;

// Variables for fitting
RooRealVar mass("m", "mass", kMASS_MIN, kMASS_MAX, "GeV/c^{2}");
mass.setRange("full", kMASS_MIN, kMASS_MAX);

RooRealVar charge("charge", "charge", -3, 3, "e");

Bool_t bool_bin = false;
vector<TString> vec_bins = {"int", "0mult1", "1mult2", "2mult3", "3mult4", "4mult5", "4mult8" ,"5mult8"};


/*
/////////////////////////////////////////////////////////////
FUNCTIONS
/////////////////////////////////////////////////////////////
*/

// This function fits the signal via RooFit
// ds_US and ds_LS are the datasets for unlike-sign and like-sign pairs respectively
// trigger (BHT1 or BHT2) and bin (e.g. "0mult3", "int", etc.) are used to load the correct parameters from the header files
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS, TString kTRIG, TString kBIN, TString kBG, TString kSIG);
