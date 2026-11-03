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
TString f_BHT2_full = "../../data/BHT2_251128.unfolded_trees.root";
TString f_BHT2_trim = "../../data/BHT2_251128.unfolded_trees.trim.root";	// made by TrimTree.cxx
TString f_BHT2 = f_BHT2_trim;	// file used by UpsFitter.cxx (f_BHT2_full for the untrimmed tree)

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
if (getenv("UPS_DCA")) s_date += Form("_dca%s", getenv("UPS_DCA"));

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

Bool_t kOFFSET = true;	// likelihood offset in the fits (better numerical precision, same minimum)

Int_t kMAX_POLISH = 5;	// maximum number of extra migrad passes while the likelihood still improves
Double_t kPOLISH_TOL = 0.001;	// improvement of the likelihood below which the polishing stops
Int_t kN_STARTS = 6;	// number of start points of the fit, the lowest likelihood of the converged fits is kept
UInt_t kSEED = 12345;	// seed of the random start points
Int_t kMAX_FIT_ATTEMPTS = 3;	// number of fit repetitions (from the last fitted values) after a failed fit

Double_t kMIN_US_ENTRIES = 100;	// fewer unlike-sign entries than this flags the fit as FAILED
Double_t kCHI2NDF_MAX = 3.0;	// chi2/ndf (US) above this flags the fit as WARN
Double_t kMEAN_MAX = 0;	// diagnostic variation meanmax_<x>: upper bound of the 1S mean when above zero
Bool_t kSEP_SCALE = false;	// variation sepscale: the 2S and 3S mass distances to the 1S mean are scaled by a floating parameter instead of fixed
Bool_t kSEP_CON = false;	// variation sepcon: the distance scale gets a Gaussian constraint from the sepscale fit of the integrated bin
Double_t kSEP_2S = 0.560180;	// 2S - 1S mass distance from embedding
Double_t kSEP_3S = 0.890502;	// 3S - 1S mass distance from embedding
Double_t kSEP_MIN = 0.5;	// bounds of the distance scale
Double_t kSEP_MAX = 1.5;
Double_t kCOUNT_LO = 8.4;	// bin counting window (clipped to the fit range)
Double_t kCOUNT_HI = 11.2;
Int_t kSHAPE_FIX = 0;	// lineshape variation: 0 floating, +1 sigma/alpha/n fixed at the BHT2_int upper bound, -1 at the lower bound
Int_t kRATIO_FIX = 0;	// 2S/3S ratio variation: same, for var_sig1_frac and var_sig2_frac

Double_t kTRIM_DCA = 3;	// the trimmed tree holds tracks up to this DCA (cm), TrimTree.cxx must be run with UPS_DCA=3
Double_t kLIMIT_FRAC = 0.001;	// parameter closer to a bound than this fraction of its range is AT_LIMIT

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
vector<TString> vec_bins = {"int", "0mult1", "1mult2", "2mult3", "3mult4", "3mult8", "4mult5", "4mult8", "5mult8"};


/*
/////////////////////////////////////////////////////////////
FUNCTIONS
/////////////////////////////////////////////////////////////
*/

// This function fits the signal via RooFit
// ds_US and ds_LS are the datasets for unlike-sign and like-sign pairs respectively
// trigger (BHT1 or BHT2) and bin (e.g. "0mult3", "int", etc.) are used to load the correct parameters from the header files
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS, TString kTRIG, TString kBIN, TString kBG, TString kSIG, TString kMULT, TString kVAR);
