#include <iostream>
#include "map"

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooHist.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"

// Namespaces
using namespace std;
using namespace RooFit;

/*
/////////////////////////////////////////////////////////////
INPUT FILES
/////////////////////////////////////////////////////////////
*/

// File locations
TString infile = "data.root";

/*
/////////////////////////////////////////////////////////////
NAME HANDLERS
/////////////////////////////////////////////////////////////
*/

TString tag = "default";	// tag for output files

// Curent date (YYMMDD) for output purposes
time_t timestamp = time(NULL);
struct tm datetime = *localtime(&timestamp);

char output[6];

strftime(output, 8, "%y%m%d", &datetime);

TString date = output;

/*
/////////////////////////////////////////////////////////////
RUN PARAMETERS
/////////////////////////////////////////////////////////////
*/

Double_t min_mass = 7.6;	// x-axis lower limit for fitting
Double_t max_mass = 13.0;	// x-axis upper limit for fitting

Double_t bin_width = 0.2;	// Bin width for drawing

Int_t sigma = 1;	// sigma for error drawing

Bool_t useLinear = true;	// Use linear error visualisation (faster, less accurate for strong correlations)

/*
/////////////////////////////////////////////////////////////
GLOBAL VARIABLES
/////////////////////////////////////////////////////////////
*/

// Data file
TFile * f_data;

// Variables for fitting
RooRealVar mass("mass", "mass", min_mass, max_mass, "GeV/c^{2}");
RooRealVar charge("charge", "charge", -3, 3, "e");

/*
/////////////////////////////////////////////////////////////
FUNCTIONS
/////////////////////////////////////////////////////////////
*/

// This function fits the signal via RooFit
// ds_US and ds_LS are the datasets for unlike-sign and like-sign pairs respectively
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS);
