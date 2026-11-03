#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "TString.h"
#include "TDatime.h"

using namespace std;

/*
/////////////////////////////////////////////////////////////
RESULT FILES (written by UpsFitter.cxx, read by Summary.cxx and Syst.cxx)
/////////////////////////////////////////////////////////////
*/

struct Par
{
	TString name;
	Double_t val, err, lo, hi;
	Bool_t at_limit;
};

struct Result
{
	Bool_t found = false;
	TString flag = "MISSING";	// OK | WARN | FAILED | MISSING
	TString reason = "";
	Double_t US = 0, LS = 0, chi2 = 0, edm = 0, SoverRoot = 0;
	Int_t attempts = 0, status = -1, covq = -1;
	Double_t S = 0, dS = 0, S1 = 0, dS1 = 0, S2 = 0, dS2 = 0, S3 = 0, dS3 = 0;
	Double_t C = 0, dC = 0;	// bin counting yield (corrected for the signal fraction in the window)
	vector<Par> pars;
};

// Bins: int plus the two multiplicity schemes
vector<TString> vec_scheme5 = {"0mult1", "1mult2", "2mult3", "3mult4", "4mult8"};
vector<TString> vec_scheme4 = {"0mult1", "1mult2", "2mult3", "3mult8"};
vector<TString> vec_sys_bins = {"int", "0mult1", "1mult2", "2mult3", "3mult4", "4mult8", "3mult8"};
vector<TString> vec_mult_bins = {"0mult1", "1mult2", "2mult3", "3mult4", "4mult8", "3mult8"};
vector<TString> vec_mult_variants = {"subcorr", "sub", "corr", "none"};

// A systematic source: one or more variations (the deviation is the maximum), cut = selection cut (Barlow criterion applies)
struct Source
{
	TString name;
	vector<TString> vars;
	Bool_t cut;
	Bool_t diag = false;	// shown in the table and plots, not part of the total
};

// Variation names as in headers/variations.hxx, "cheby" is the nominal fit with the Chebyshev background
static vector<Source> MakeSources()
{
	vector<Source> sources;
	sources.push_back({"bg", {"cheby"}, false});
	sources.push_back({"range", {"range_lo_up", "range_lo_dn", "range_hi_up", "range_hi_dn"}, false});
	sources.push_back({"lineshape", {"lineshape_up", "lineshape_dn"}, false});
	sources.push_back({"ratio", {"ratio_up", "ratio_dn"}, false});
	sources.push_back({"mean", {"meanfree"}, false});
	sources.push_back({"sep", {"sepscale"}, false, true});
	sources.push_back({"sepcon", {"sepcon"}, false});
	sources.push_back({"count", {"count"}, false});

	vector<TString> cuts = {"nHitsFit", "nHitsRatio", "DCA", "nSigmaE_lo", "nSigmaE_hi", "ETOWoECLU", "EoP_lo", "EoP_hi", "Dphi", "DZ", "ADC"};
	vector<TString> noloose = {"DCA", "nSigmaE_hi", "ADC"};
	for (auto c : cuts)
	{
		Source s;
		s.name = c;
		s.cut = true;
		s.vars.push_back("cut_" + c + "_tight");
		if (find(noloose.begin(), noloose.end(), c) == noloose.end()) s.vars.push_back("cut_" + c + "_loose");
		sources.push_back(s);
	}
	return sources;
}

// All fitted variations, nominal first ("count" is the nominal fit read as bin counting, not a fit of its own)
static vector<TString> MakeVars()
{
	vector<TString> vars = {"nom"};
	for (auto s : MakeSources())
		for (auto v : s.vars)
			if (v != "count") vars.push_back(v);
	return vars;
}

// Today as YYMMDD, the date used in the result file names
static TString TodayYYMMDD()
{
	TString s = Form("%06d", TDatime().GetDate() % 1000000);
	if (getenv("UPS_DCA")) s += Form("_dca%s", getenv("UPS_DCA"));
	return s;
}

// Same name as built in UpsFitter.cxx: <trig>_<bin>_<date>_<tag>.txt
static TString ResultFile(TString bin, TString var, TString bg, TString mult, TString date)
{
	TString tag = "default";
	if (bin != "int") tag += "_" + mult;
	if (bg != "exp") tag += "_" + bg;
	if (var != "nom") tag += "_" + var;
	return "output/results/BHT2_" + bin + "_" + date + "_" + tag + ".txt";
}

static Result ReadResult(TString file)
{
	Result r;
	ifstream in(file.Data());
	if (!in.is_open()) return r;
	r.found = true;
	r.flag = "FAILED";
	r.reason = " no_FIT_STATUS";

	string line;
	while (getline(in, line))
	{
		istringstream ss(line);
		string key, a, b;
		ss >> key >> a >> b;
		if (key == "US_entries:") r.US = atof(a.c_str());
		else if (key == "LS_entries:") r.LS = atof(a.c_str());
		else if (key == "chi2ndf_US:") r.chi2 = atof(a.c_str());
		else if (key == "fit_attempts:") r.attempts = atoi(a.c_str());
		else if (key == "minuit_status:") r.status = atoi(a.c_str());
		else if (key == "cov_quality:") r.covq = atoi(a.c_str());
		else if (key == "edm:") r.edm = atof(a.c_str());
		else if (key == "S:") { r.S = atof(a.c_str()); r.dS = atof(b.c_str()); }
		else if (key == "S_1S:") { r.S1 = atof(a.c_str()); r.dS1 = atof(b.c_str()); }
		else if (key == "S_2S:") { r.S2 = atof(a.c_str()); r.dS2 = atof(b.c_str()); }
		else if (key == "S_3S:") { r.S3 = atof(a.c_str()); r.dS3 = atof(b.c_str()); }
		else if (key == "COUNT:") { r.C = atof(a.c_str()); r.dC = atof(b.c_str()); }
		else if (key == "S/sqrt(S+B):") r.SoverRoot = atof(a.c_str());
		else if (key == "PAR:")
		{
			Par p;
			string c, d, e;
			istringstream sp(line);
			sp >> key >> a >> b >> c >> d >> e;
			p.name = a;
			p.val = atof(b.c_str());
			p.err = atof(c.c_str());
			p.lo = atof(d.c_str());
			p.hi = atof(e.c_str());
			string f;
			sp >> f;
			p.at_limit = (f == "AT_LIMIT");
			r.pars.push_back(p);
		}
		else if (key == "FIT_STATUS:")
		{
			r.flag = TString(a).ReplaceAll(":", "");
			size_t pos = line.find(a);
			r.reason = (pos == string::npos) ? "" : TString(line.substr(pos + a.size()));
		}
	}
	return r;
}

// Results of all systematics bins for one variation
static map<TString, Result> LoadVariation(TString var, TString mult, TString date)
{
	TString bg = "exp";
	TString variation = var;
	if (var == "cheby") { bg = "cheby"; variation = "nom"; }

	if (var == "count") variation = "nom";

	map<TString, Result> res;
	for (auto bin : vec_sys_bins)
	{
		res[bin] = ReadResult(ResultFile(bin, variation, bg, mult, date));
		// Bin counting instead of the fitted yield, same fit status
		if (var == "count") { res[bin].S = res[bin].C; res[bin].dS = res[bin].dC; }
	}
	return res;
}

// A result that can enter the systematics
static Bool_t Usable(const Result &r)
{
	return r.found && r.flag != "FAILED" && r.flag != "MISSING" && std::isfinite(r.S) && std::isfinite(r.dS) && r.S > 0 && r.dS > 0;
}

/*
/////////////////////////////////////////////////////////////
DEVIATION OF THE SELF-NORMALISED YIELD
/////////////////////////////////////////////////////////////
*/

// The yield fraction f = N_bin/N_int is formed within each variation (int is run for every variation)
// UpsFitter has no event counts, so the event normalisation of r_i is not applied: it is the same for nominal and variation,
// and the relative deviation rel = f_var/f_nom - 1 does not depend on it
// sig is the Barlow sigma of the relative deviation, sqrt(|sigma_nom^2 - sigma_var^2|)/f_nom, with N_bin and N_int treated as uncorrelated
struct Delta
{
	Bool_t ok = false;
	Double_t rel = 0;
	Double_t sig = 0;
};

static Delta MakeDelta(const Result &int_nom, const Result &bin_nom, const Result &int_var, const Result &bin_var)
{
	Delta d;
	if (!Usable(int_nom) || !Usable(bin_nom) || !Usable(int_var) || !Usable(bin_var)) return d;

	Double_t f_nom = bin_nom.S/int_nom.S;
	Double_t f_var = bin_var.S/int_var.S;
	Double_t s_nom = f_nom*sqrt(pow(bin_nom.dS/bin_nom.S, 2) + pow(int_nom.dS/int_nom.S, 2));
	Double_t s_var = f_var*sqrt(pow(bin_var.dS/bin_var.S, 2) + pow(int_var.dS/int_var.S, 2));

	d.ok = true;
	d.rel = f_var/f_nom - 1;
	d.sig = sqrt(fabs(s_nom*s_nom - s_var*s_var))/f_nom;
	return d;
}

/////////////////////////////////////////////////////////////
