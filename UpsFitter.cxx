#include "headers/common.hxx"
#include "headers/setup.hxx"
#include "headers/cuts.hxx"
#include "headers/variations.hxx"
#include "headers/results.hxx"
//#include "pars/parsBHT1.hxx"
#include "pars/parsBHT2.hxx"

/////////////////////////////////////////////////////////////
void UpsFitter(TString kTRIG = "BHT2", TString kBIN = "int", TString kBG = "exp", TString kSIG = "CBShape", TString kMULT = "subcorr", TString kVAR = "nom")
{
	// Load the data file
	if (kTRIG == "BHT1")
	{
		f_IN = new TFile(f_BHT1, "READ");
	}
	else if (kTRIG == "BHT2")
	{
		f_IN = new TFile(f_BHT2, "READ");
	}
	else
	{
		cout << "Invalid trigger specified. Exiting..." << endl;
		return;
	}

	// Check if the bin makes sense
	for (auto b : vec_bins)
		if (b == kBIN) bool_bin = true;
	if (!bool_bin)
	{
		cout << "Invalid bin specified. Exiting..." << endl;
		return;
	}

	if (kMULT == "none") s_mult = s_mult_none;
	else if (kMULT == "sub") s_mult = s_mult_sub;
	else if (kMULT == "corr") s_mult = s_mult_corr;
	else if (kMULT == "subcorr") s_mult = s_mult_subcorr;
	else
	{
		cout << "Invalid TofMult variant (none|sub|corr|subcorr). Exiting..." << endl;
		return;
	}

	if (kBIN == "int")
		s_cut_mult = "";
	else if (kBIN.Contains("mult"))
	{
		s_tag += "_" + kMULT;
		s_cut_mult = Form(" && (%s > 3*%s) && (%s <= 3*%s)",
						s_mult.Data(), TString(kBIN(0, kBIN.Index("mult"))).Data(),
						s_mult.Data(), TString(kBIN(kBIN.Index("mult") + 4, kBIN.Length())).Data());
	}
	else
	{
		cout << "Invalid bin format. Exiting..." << endl;
		return;
	}

	if (kBG != "exp") s_tag += "_" + kBG;

	// Apply the systematic variation (nom = nominal)
	if (!SetVariation(kVAR))
	{
		cout << "Invalid or not implemented variation " << kVAR << ". Exiting..." << endl;
		return;
	}
	if (kVAR != "nom") s_tag += "_" + kVAR;

	if (f_BHT2 == f_BHT2_trim && kDCA_HI > kTRIM_DCA)
	{
		cout << "DCA cut " << kDCA_HI << " is looser than the trimmed tree (" << kTRIM_DCA << "). Exiting..." << endl;
		return;
	}

	// Initialize cuts

	TString s_cuts = MakeCuts(kTRIG);

	s_cuts += s_cut_mult;

	cout << "Applied cuts: " << s_cuts << endl;

	// Load datasets for the specified bin
	TTree * tree = (TTree *)f_IN->Get("tree");

	TTree * tree_cut = (TTree *)tree->CopyTree(s_cuts);

	tree_cut->Print();
		
	RooDataSet * ds_US = new RooDataSet("ds_US_" + kBIN, "ds_US_" + kBIN, RooArgSet(mass, charge), Import(*tree_cut), Cut("charge==0"));
	ds_US = (RooDataSet *)ds_US->reduce(RooArgSet(mass));
	ds_US->Print();
	RooDataSet * ds_LS = new RooDataSet("ds_LS_" + kBIN, "ds_LS_" + kBIN, RooArgSet(mass, charge), Import(*tree_cut), Cut("charge!=0"));
	ds_LS = (RooDataSet *)ds_LS->reduce(RooArgSet(mass));
	ds_LS->Print();

	// This calls the fitting function
	fit_bin(ds_US, ds_LS, kTRIG, kBIN, kBG, kSIG, kMULT, kVAR);

	// Clean up
	delete ds_US;
	delete ds_LS;
	f_IN->Close();
	cout << "Fit for " << kTRIG << "_" << kBIN << " completed!" << endl;

	return;
}

/////////////////////////////////////////////////////////////
// TODO: rename the parameters sensibly and clean this up
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS, TString kTRIG, TString kBIN, TString kBG, TString kSIG, TString kMULT, TString kVAR)
{
	// Save number of entries
	double US_entries = ds_US->sumEntries();
	double LS_entries = ds_LS->sumEntries();

	// Define a sample for simultaneous fitting
	RooCategory sample("sample", "sample");
	sample.defineType("US");
	sample.defineType("LS");

	// Simultaneous model
	RooSimultaneous * simPdf = new RooSimultaneous("simPdf", "simPdf", sample);

	// Define a combined dataset 
	RooDataSet * ds_comb = new RooDataSet("ds_comb", "Combined dataset", RooArgList(mass), Index(sample), Import("US", *ds_US), Import("LS", *ds_LS));
	
	// Get parameters for this bin and trigger
	TString key_bin = kTRIG + "_" + kBIN;
	if (!map_par_sig.count(key_bin) || !map_frac.count(key_bin))
	{
		cout << "No parameters for " << key_bin << ", falling back to " << kTRIG << "_int" << endl;
		key_bin = kTRIG + "_int";
	}
	vector<vector<Double_t>> par_sig = map_par_sig[key_bin];
	TString key_bg = kTRIG + "_" + kBG + "_" + kBIN;
	if (!map_par_bg.count(key_bg)) key_bg = kTRIG + "_" + kBG;
	vector<vector<Double_t>> par_bg = map_par_bg[key_bg];
	vector<vector<Double_t>> frac = map_frac[key_bin];


	// Load parameters and create PDFs
	// Background
	RooAbsPdf * pdf_bg_comb;
	RooAbsPdf * pdf_bg_res;
	if (kBG == "exp")
	{
		// par_bg: comb exponent, res exponent
		RooRealVar * var_bg_comb = new RooRealVar("var_bg_comb","Combiantorial background exponent", par_bg[0][0], par_bg[0][1], par_bg[0][2]);
		RooRealVar * var_bg_res = new RooRealVar("var_bg_res","Residual background exponent", par_bg[1][0], par_bg[1][1], par_bg[1][2]);

		pdf_bg_comb = new RooExponential("pdf_bg_comb","Background, combinatorial, exponential", mass, *var_bg_comb);
		pdf_bg_res = new RooExponential("pdf_bg_res","Background, residual, exponential", mass, *var_bg_res);
	}
	else if (kBG == "cheby")
	{
		// par_bg: comb c1, comb c2, res c1, res c2 (second order Chebyshev)
		RooRealVar * var_bg_comb_c1 = new RooRealVar("var_bg_comb_c1","Combiantorial background Chebyshev c1", par_bg[0][0], par_bg[0][1], par_bg[0][2]);
		RooRealVar * var_bg_comb_c2 = new RooRealVar("var_bg_comb_c2","Combiantorial background Chebyshev c2", par_bg[1][0], par_bg[1][1], par_bg[1][2]);
		RooRealVar * var_bg_res_c1 = new RooRealVar("var_bg_res_c1","Residual background Chebyshev c1", par_bg[2][0], par_bg[2][1], par_bg[2][2]);
		RooRealVar * var_bg_res_c2 = new RooRealVar("var_bg_res_c2","Residual background Chebyshev c2", par_bg[3][0], par_bg[3][1], par_bg[3][2]);

		pdf_bg_comb = new RooChebychev("pdf_bg_comb","Background, combinatorial, Chebyshev", mass, RooArgList(*var_bg_comb_c1, *var_bg_comb_c2));
		pdf_bg_res = new RooChebychev("pdf_bg_res","Background, residual, Chebyshev", mass, RooArgList(*var_bg_res_c1, *var_bg_res_c2));
	}
	else
	{
		cout << "Invalid background type (exp|cheby). Exiting..." << endl;
		return;
	}

	Double_t mean_int = 0;	// common 1S mean from the nominal integrated bin, 0 in the integrated bin

	// Signal
	//if (kSIG == "CBShape")
	//{
		RooRealVar var_sig1_mu("var_sig1_mu","Signal 1 #mu", par_sig[0][0], par_sig[0][1], (kMEAN_MAX > 0) ? kMEAN_MAX : par_sig[0][2]);
		// The 1S mean is the same in all bins: held at the value fitted in the nominal integrated bin (the integrated bin itself floats).
		// Variation meanfree: the mean floats in every bin, started from the integrated-bin value
		if (kBIN != "int")
		{
			Result r_int = ReadResult(ResultFile("int", "nom", kBG, kMULT, s_date));
			for (auto p : r_int.pars)
				if (p.name == "var_sig1_mu") mean_int = p.val;
			if (!r_int.found || r_int.flag == "FAILED" || mean_int <= 0)
			{
				cout << "No usable nominal fit of the integrated bin for the common 1S mean. Exiting..." << endl;
				return;
			}
			var_sig1_mu.setVal(mean_int);
			if (kVAR != "meanfree") var_sig1_mu.setConstant(true);
		}
		RooRealVar var_sig1_sigma("var_sig1_sigma","Signal 1 #sigma", par_sig[1][0], par_sig[1][1], par_sig[1][2]);
		RooRealVar var_sig1_alpha("var_sig1_alpha","Signal 1 #alpha", par_sig[2][0], par_sig[2][1], par_sig[2][2]);
		RooRealVar var_sig1_n("var_sig1_n","Signal 1 #n", par_sig[3][0], par_sig[3][1], par_sig[3][2]);

		RooRealVar var_sig2_mu("var_sig2_mu","Signal 2 #mu", par_sig[4][0], par_sig[4][1], par_sig[4][2]);
		RooRealVar var_sig2_sigma("var_sig2_sigma","Signal 2 #sigma", par_sig[5][0], par_sig[5][1], par_sig[5][2]);
		RooRealVar var_sig2_alpha("var_sig2_alpha","Signal 2 #alpha", par_sig[6][0], par_sig[6][1], par_sig[6][2]);
		RooRealVar var_sig2_n("var_sig2_n","Signal 2 #n", par_sig[7][0], par_sig[7][1], par_sig[7][2]);

		RooRealVar var_sig3_mu("var_sig3_mu","Signal 3 #mu", par_sig[8][0], par_sig[8][1], par_sig[8][2]);
		RooRealVar var_sig3_sigma("var_sig3_sigma","Signal 3 #sigma", par_sig[9][0], par_sig[9][1], par_sig[9][2]);
		RooRealVar var_sig3_alpha("var_sig3_alpha","Signal 3 #alpha", par_sig[10][0], par_sig[10][1], par_sig[10][2]);
		RooRealVar var_sig3_n("var_sig3_n","Signal 3 #n", par_sig[11][0], par_sig[11][1], par_sig[11][2]);

		/*
		/////////////////////
		!!!!   WARNING   !!!!
		/////////////////////
		*/

		// Bound means for 2S and 3S based on embedding
		// THIS SUPERCEDES SIG2_MU AND SIG3_MU VALUES WITH BOUND ONES

		// Option sepscale: the distances are scaled by one floating parameter (linear), otherwise they are fixed
		RooRealVar var_sep_scale("var_sep_scale", "Scale of the 2S and 3S mass distances", 1.0, kSEP_MIN, kSEP_MAX);
		var_sep_scale.setConstant(!kSEP_SCALE);
		// Variation sepcon: Gaussian constraint on the scale from the sepscale fit of the integrated bin (the integrated bin itself is left unconstrained)
		RooArgSet set_constraint;
		if (kSEP_CON && kBIN != "int")
		{
			Result r_sep = ReadResult(ResultFile("int", "sepscale", kBG, kMULT, s_date));
			Double_t sep_mu = 0, sep_sigma = 0;
			for (auto p : r_sep.pars)
				if (p.name == "var_sep_scale") { sep_mu = p.val; sep_sigma = p.err; }
			if (!r_sep.found || r_sep.flag == "FAILED" || sep_sigma <= 0)
			{
				cout << "No usable sepscale fit of the integrated bin for sepcon. Exiting..." << endl;
				return;
			}
			RooConstVar * con_mu = new RooConstVar("con_sep_mu", "Constraint mean", sep_mu);
			RooConstVar * con_sigma = new RooConstVar("con_sep_sigma", "Constraint width", sep_sigma);
			set_constraint.add(*(new RooGaussian("con_sep", "Gaussian constraint on the distance scale", var_sep_scale, *con_mu, *con_sigma)));
			cout << "sepcon: scale constrained to " << sep_mu << " +- " << sep_sigma << endl;
		}
		RooFormulaVar form_sig2_mu_bound("form_sig2_mu", "Signal 2 #mu - bound", Form("@0 + @1*%.6f", kSEP_2S), RooArgList(var_sig1_mu, var_sep_scale));
		RooFormulaVar form_sig3_mu_bound("form_sig3_mu", "Signal 3 #mu - bound", Form("@0 + @1*%.6f", kSEP_3S), RooArgList(var_sig1_mu, var_sep_scale));

		/*
		/////////////////////
		YOU WERE WARNED ¯\_(ツ)_/¯
		/////////////////////
		*/

		RooCBShape * pdf_sig1 = new RooCBShape("pdf_sig1", "Signal 1 - Crystall Ball, onesided", mass, var_sig1_mu, var_sig1_sigma, var_sig1_alpha, var_sig1_n);
		RooCBShape * pdf_sig2 = new RooCBShape("pdf_sig2", "Signal 2 - Crystall Ball, onesided", mass, form_sig2_mu_bound, var_sig2_sigma, var_sig2_alpha, var_sig2_n);
		RooCBShape * pdf_sig3 = new RooCBShape("pdf_sig3", "Signal 3 - Crystall Ball, onesided", mass, form_sig3_mu_bound, var_sig3_sigma, var_sig3_alpha, var_sig3_n);
	//}
	//else
	//{
	//	cout << "Invalid signal shape. Exiting..." << endl;
	//	return;
	//}

	/*
	RooCrystalBall *gauss1 = new RooCrystalBall("gauss1", "Gaussian signal 1", mass, sig1_mu, sig1_sigma, sig1_alphaL, sig1_alphaR, sig1_nL, sig1_nR);
	RooCrystalBall *gauss2 = new RooCrystalBall("gauss2", "Gaussian signal 2", mass, sig2_mu, sig2_sigma, sig2_alphaL, sig2_alphaR, sig2_nL, sig2_nR);
	RooCrystalBall *gauss3 = new RooCrystalBall("gauss3", "Gaussian signal 3", mass, sig3_mu, sig3_sigma, sig3_alphaL, sig3_alphaR, sig3_nL, sig3_nR);
	*/
	
	// Recursive fractions
	RooRealVar var_bg_comb_frac("var_bg_comb_frac","Fraction of combinatorial background", frac[0][0], frac[0][1], frac[0][2]);
	RooRealVar var_bg_res_frac("var_bg_res_frac","Fraction of residual background", frac[1][0], frac[1][1], frac[1][2]);
	RooRealVar var_sig1_frac("var_sig1_frac","Fraction of signal 1", frac[2][0], frac[2][1], frac[2][2]);
	RooRealVar var_sig2_frac("var_sig2_frac","Fraction of signal 2", frac[3][0], frac[3][1], frac[3][2]);

	/*
	/////////////////////
	!!!!   WARNING   !!!!
	/////////////////////

	THESE SETTING CHANGE PARAMETER VALUES FOUND IN THE PAR HEADERS
	DO NOT BE SURPRISED THE CODE DOES NOT PICK THEM UP!!!
	*/

	// Ensure the combinatorial background describes the LS data
	// THIS SUPERCEDES FRACTION SETTING FOR COMBINATORIAL BACKGROUND

	var_bg_comb_frac.setVal(ds_LS->sumEntries()/ds_US->sumEntries());
	var_bg_comb_frac.setConstant(true);

	/*
	/////////////////////
	YOU WERE WARNED ¯\_(ツ)_/¯
	/////////////////////
	*/

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	//if (kBG == "exp" && kSIG == "CBShape")
		RooAddPdf * pdf_total_US = new RooAddPdf("pdf_total_US", "Total PDF for US data", RooArgList(*pdf_bg_comb, *pdf_bg_res, *pdf_sig1, *pdf_sig2, *pdf_sig3), RooArgList(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac, var_sig2_frac), kTRUE);
	//else
	//{
	//	cout << "Invalid signal/background commbination. Exitting..." << endl;
	//	return;
	//}

	// Add PDFs to the simultaneous model
	simPdf->addPdf(*pdf_bg_comb, "LS");
	simPdf->addPdf(*pdf_total_US, "US");

	// Variations start from the fitted values of the nominal fit of the same bin (same date, background and TofMult)
	// Fixed parameters keep their values, the lineshape and ratio variations set theirs below
	if (kVAR != "nom")
	{
		Result r_nom = ReadResult(ResultFile(kBIN, "nom", kBG, kMULT, s_date));
		if (r_nom.found && r_nom.flag != "FAILED")
		{
			RooArgSet * pars_all = simPdf->getParameters(*ds_comb);
			int n_seed = 0;
			for (auto p : r_nom.pars)
			{
				RooRealVar * v = (RooRealVar *)pars_all->find(p.name);
				if (!v || v->isConstant()) continue;
				// The result files store six digits, keep the value inside the bounds
				v->setVal(std::min(std::max(p.val, v->getMin()), v->getMax()));
				n_seed++;
			}
			cout << "Start values of " << n_seed << " parameters taken from the nominal fit" << endl;
			delete pars_all;
		}
		else cout << "No usable nominal fit for " << kBIN << ", starting from the parameter file values" << endl;
	}

	// Fit the model to the combined dataset
	// Stage 1: lineshape and signal fractions held at their start values, only backgrounds and the 1S mean float
	// Stage 2 (below) releases them within their bounds, starting from the stage 1 result
	// Lineshape and ratio variations: fixed at the lower/upper bound of the BHT2_int (embedding) entry and never released
	vector<vector<Double_t>> par_int = map_par_sig[kTRIG + "_int"];
	vector<vector<Double_t>> frac_int = map_frac[kTRIG + "_int"];
	if (kSHAPE_FIX != 0)
	{
		int bound = (kSHAPE_FIX > 0) ? 2 : 1;
		vector<RooRealVar *> vars_shape = {&var_sig1_sigma, &var_sig1_alpha, &var_sig1_n, &var_sig2_sigma, &var_sig2_alpha, &var_sig2_n, &var_sig3_sigma, &var_sig3_alpha, &var_sig3_n};
		vector<int> idx_shape = {1, 2, 3, 5, 6, 7, 9, 10, 11};
		for (int i = 0; i < 9; i++)
		{
			vars_shape[i]->setVal(par_int[idx_shape[i]][bound]);
			vars_shape[i]->setConstant(true);
		}
	}
	if (kRATIO_FIX != 0)
	{
		int bound = (kRATIO_FIX > 0) ? 2 : 1;
		var_sig1_frac.setVal(frac_int[2][bound]);
		var_sig1_frac.setConstant(true);
		var_sig2_frac.setVal(frac_int[3][bound]);
		var_sig2_frac.setConstant(true);
	}

	RooArgSet set_shape;
	for (auto var : {&var_sig1_sigma, &var_sig1_alpha, &var_sig1_n, &var_sig2_sigma, &var_sig2_alpha, &var_sig2_n, &var_sig3_sigma, &var_sig3_alpha, &var_sig3_n, &var_sig1_frac, &var_sig2_frac})
		if (!var->isConstant()) set_shape.add(*var);

	RooCmdArg cmd_con = (set_constraint.getSize() > 0) ? ExternalConstraints(set_constraint) : RooCmdArg::none();

	// Likelihood without offset, the same function for every pass so that the passes and starts can be compared
	RooAbsReal * nll = simPdf->createNLL(*ds_comb, Extended(true), cmd_con);

	// Minuit status 1 (covariance made positive-definite) with covariance quality 3 is a usable fit
	auto usable = [](RooFitResult * r) { return (r->status() == 0 || r->status() == 1) && r->covQual() == 3; };

	// Full fit, all parameters floating
	auto run_migrad = [&]()
	{
		return simPdf->fitTo(*ds_comb,
								Extended(true),
								Minimizer("Minuit2", "migrad"),
								Offset(kOFFSET),
								MaxCalls(100000),
								Strategy(2),
								EvalErrorWall(1),
								RecoverFromUndefinedRegions(1000),
								cmd_con,
								Save()
							);
	};

	// One complete fit from the current start values: stage 1, full fit with retries, polishing
	Int_t fit_attempts = 0;
	Int_t polish_passes = 0;
	auto fit_chain = [&](Double_t &nll_out)
	{
		for (auto par : set_shape) ((RooRealVar *)par)->setConstant(true);
		simPdf->fitTo(*ds_comb,
						Extended(true),
						Minimizer("Minuit2", "migrad"),
						Offset(kOFFSET),
						MaxCalls(100000),
						Strategy(2),
						EvalErrorWall(1),
						RecoverFromUndefinedRegions(1000),
						cmd_con,
						PrintLevel(-1)
					);
		for (auto par : set_shape) ((RooRealVar *)par)->setConstant(false);

		// A failed fit is repeated from the last fitted values (the parameters keep them), up to kMAX_FIT_ATTEMPTS times
		RooFitResult * res = 0;
		fit_attempts = 0;
		while (fit_attempts < kMAX_FIT_ATTEMPTS)
		{
			delete res;

			// Retry: a gradient-free simplex pass first moves the parameters off a point where the numerical gradient is unusable
			if (fit_attempts > 0)
				simPdf->fitTo(*ds_comb, Extended(true), Minimizer("Minuit2", "simplex"), Offset(kOFFSET), MaxCalls(100000), cmd_con, PrintLevel(-1));

			res = run_migrad();
			fit_attempts++;
			if (usable(res)) break;
			cout << "Fit attempt " << fit_attempts << " failed with status " << res->status() << endl;
		}

		// Polish: the likelihood is flat along some lineshape directions and migrad can stop early, so repeat it from the last point
		// until the likelihood no longer improves; a pass that does not converge is dropped
		nll_out = nll->getVal();
		polish_passes = 0;
		while (usable(res) && polish_passes < kMAX_POLISH)
		{
			RooArgSet * pars_now = simPdf->getParameters(*ds_comb);
			RooArgSet pars_saved;
			pars_now->snapshot(pars_saved);
			RooFitResult * res_next = run_migrad();
			polish_passes++;
			Double_t nll_next = nll->getVal();
			if (!usable(res_next) || nll_next > nll_out)
			{
				// Back to the last good point, its covariance is the one in res
				pars_now->assignValueOnly(pars_saved);
				delete res_next;
				delete pars_now;
				break;
			}
			Double_t gain = nll_out - nll_next;
			delete res;
			res = res_next;
			nll_out = nll_next;
			delete pars_now;
			if (gain < kPOLISH_TOL) break;
		}
		return res;
	};

	// The likelihood has nearly degenerate minima along the lineshape directions, so the fit is repeated from kN_STARTS start points
	// (the first is the parameter file or nominal fit, the others are random points inside the bounds of the floating lineshape
	// parameters, fixed seed) and the lowest likelihood among the converged fits is kept
	RooArgSet * pars_all = simPdf->getParameters(*ds_comb);
	vector<RooRealVar *> vars_start;
	for (auto par : set_shape) vars_start.push_back((RooRealVar *)par);
	if (!var_sig1_mu.isConstant()) vars_start.push_back(&var_sig1_mu);
	if (!var_sep_scale.isConstant()) vars_start.push_back(&var_sep_scale);

	TRandom3 rnd(kSEED);
	RooFitResult * result = 0;
	Double_t nll_final = 0;
	Int_t best_start = -1;
	Int_t attempts_best = 0;
	Int_t polish_best = 0;
	for (int k = 0; k < kN_STARTS; k++)
	{
		if (k > 0)
			for (auto v : vars_start) v->setVal(v->getMin() + (0.1 + 0.8*rnd.Rndm())*(v->getMax() - v->getMin()));

		Double_t nll_k;
		RooFitResult * res_k = fit_chain(nll_k);
		bool good_k = usable(res_k);
		cout << "Start " << k << ": status " << res_k->status() << ", covQual " << res_k->covQual() << ", nll " << Form("%.4f", nll_k) << endl;

		bool better = !result || (good_k && (best_start < 0 || nll_k < nll_final - kPOLISH_TOL));
		if (better)
		{
			delete result;
			result = res_k;
			nll_final = nll_k;
			attempts_best = fit_attempts;
			polish_best = polish_passes;
			if (good_k) best_start = k;
		}
		else delete res_k;
	}

	// The parameters hold the last start, restore the kept fit
	RooArgList pars_best = result->floatParsFinal();
	for (int i = 0; i < pars_best.getSize(); i++)
	{
		RooRealVar * p = (RooRealVar *)pars_best.at(i);
		RooRealVar * v = (RooRealVar *)pars_all->find(p->GetName());
		if (v) { v->setVal(p->getVal()); v->setError(p->getError()); }
	}
	delete pars_all;
	fit_attempts = attempts_best;
	polish_passes = polish_best;
	cout << "Kept start " << best_start << ", nll = " << Form("%.4f", nll_final) << endl;
	delete nll;

	// Compute yields
	RooGenericPdf yield1("yield1", "yield from signal fit", "(1.0-var_bg_comb_frac)*(1.0-var_bg_res_frac)", RooArgSet(var_bg_comb_frac, var_bg_res_frac));
	RooGenericPdf yield2("yield2", "yield from signal 2 fit", "(1.0-var_bg_comb_frac)*(1.0-var_bg_res_frac)*(1.0-var_sig1_frac)", RooArgSet(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac));
	RooGenericPdf yield3("yield3", "yield from signal 3 fit", "(1.0-var_bg_comb_frac)*(1.0-var_bg_res_frac)*(1.0-var_sig1_frac)*(1.0-var_sig2_frac)", RooArgSet(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac, var_sig2_frac));
	RooGenericPdf yield_1S("yield_1S", "yield of Upsilon(1S)", "(1.0-var_bg_comb_frac)*(1.0-var_bg_res_frac)*var_sig1_frac", RooArgSet(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac));
	RooGenericPdf yield_2S("yield_2S", "yield of Upsilon(2S)", "(1.0-var_bg_comb_frac)*(1.0-var_bg_res_frac)*(1.0-var_sig1_frac)*var_sig2_frac", RooArgSet(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac, var_sig2_frac));

	/*
	cout << "US entries: " << US_entries << ", LS entries: " << LS_entries << endl << endl;
	cout << "yield frac: " << yield1.getVal() << ", yield: " << (yield1.getVal()*US_entries) << endl;
	cout << "error frac: " << yield1.getPropagatedError(*result) << ", error: " << (yield1.getPropagatedError(*result)*US_entries) << endl << endl;
	cout << "yield frac: " << yield1.getVal() << ", yield: " << ((yield1.getVal()-yield2.getVal()-yield3.getVal())*US_entries) << endl;
	cout << "yield2 frac: " << yield2.getVal() << ", yield2: " << ((yield2.getVal()-yield3.getVal()) * US_entries) << endl;
	cout << "yield3 frac: " << yield3.getVal() << ", yield3: " << (yield3.getVal() * US_entries) << endl;
	*/

	// Create a canvas for plotting
	TCanvas *cnv = new TCanvas("cnv", "cnv", 1600, 1200);

	// Plot simulated data
	cnv->cd();
	RooPlot* frame = mass.frame(((kMASS_MAX - kMASS_MIN)/kBIN_WIDTH));
	frame->SetTitle("");
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitleOffset(0.9);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.05);
	frame->GetXaxis()->SetTitleOffset(0.85);
	frame->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.02);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.05);
	ds_LS->plotOn(frame, MarkerColor(kRed), LineColor(kRed));
	ds_US->plotOn(frame);
	pdf_total_US->plotOn(frame, VisualizeError(*result, kSIGMA, kLINEAR_ERR), DrawOption("F"), FillStyle(3001));
	pdf_total_US->plotOn(frame);
	/*
	total_sim->plotOn(frame, Components(*expLS), LineColor(kRed), LineStyle(kDashed), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kRed), FillStyle(3001));
	total_sim->plotOn(frame, Components(*expLS), LineColor(kRed), LineStyle(kDashed));
	total_sim->plotOn(frame, Components(*expUS), LineColor(kGray), LineStyle(kDashed), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kGray), FillStyle(3001));
	total_sim->plotOn(frame, Components(*expUS), LineColor(kGray), LineStyle(kDashed));
	total_sim->plotOn(frame, Components(*gauss1), LineColor(kMagenta), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kMagenta), FillStyle(3001));
	total_sim->plotOn(frame, Components(*gauss1), LineColor(kMagenta), LineStyle(kSolid));
	total_sim->plotOn(frame, Components(*gauss2), LineColor(kMagenta+2), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kMagenta+2), FillStyle(3001));
	total_sim->plotOn(frame, Components(*gauss2), LineColor(kMagenta+2), LineStyle(kSolid));
	total_sim->plotOn(frame, Components(*gauss3), LineColor(kMagenta+4), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kMagenta+4), FillStyle(3001));
	total_sim->plotOn(frame, Components(*gauss3), LineColor(kMagenta+4), LineStyle(kSolid));
	total_sim->plotOn(frame, Components(RooArgSet(*gauss1, *gauss2, *gauss3)), LineColor(kCyan), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kTeal), FillStyle(3001));
	total_sim->plotOn(frame, Components(RooArgSet(*gauss1, *gauss2, *gauss3)), LineColor(kCyan), LineStyle(kSolid));
	total_sim->plotOn(frame, Components(RooArgSet(*expLS, *expUS)), LineColor(kOrange), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kOrange), FillStyle(3001));
	total_sim->plotOn(frame, Components(RooArgSet(*expLS, *expUS)), LineColor(kOrange), LineStyle(kSolid));
	total_sim->plotOn(frame);
	*/
	frame->Draw("F");
	//--------------------------------------------

	// Calculate chi2 after extracting the number of parameters - used to calculate Ndof (replaced below by the bin-by-bin Pearson chi2)
	Int_t nparam_US = simPdf->getParameters(*ds_US)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_US = " << nparam_US << endl;
	Double_t chi2ndf_US = frame->chiSquare(nparam_US);

	// LS: combinatorial background pdf against the LS data, own frame (the chi2 of a frame uses its last curve and histogram)
	RooPlot * frame_LS = mass.frame((kMASS_MAX - kMASS_MIN)/kBIN_WIDTH);
	ds_LS->plotOn(frame_LS);
	pdf_bg_comb->plotOn(frame_LS, Normalization(LS_entries, RooAbsReal::NumEvent));
	Int_t nparam_LS = pdf_bg_comb->getParameters(*ds_LS)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_LS = " << nparam_LS << endl;
	Double_t chi2ndf_LS = frame_LS->chiSquare(nparam_LS);
	delete frame_LS;

	// Diagnostic: US chi2 bin by bin, expected counts from the integral of the fitted pdf over each bin
	Int_t nb_diag = (Int_t)TMath::Ceil((kMASS_MAX - kMASS_MIN)/kBIN_WIDTH - 1e-6);	// the last bin is clipped to the fit range
	vector<Double_t> n_obs(nb_diag, 0.0);
	for (int i = 0; i < ds_US->numEntries(); i++)
	{
		Double_t m_i = ds_US->get(i)->getRealValue("m");
		Int_t ib = (Int_t)((m_i - kMASS_MIN)/kBIN_WIDTH);
		if (ib >= 0 && ib < nb_diag) n_obs[ib] += ds_US->weight();
	}
	Double_t chi2_pearson = 0, chi2_dev = 0, mu_tot = 0;
	Int_t n_low = 0;
	cout << "Bin-by-bin US: m_lo, observed, expected, pull" << endl;
	for (int ib = 0; ib < nb_diag; ib++)
	{
		mass.setRange("diag", kMASS_MIN + ib*kBIN_WIDTH, std::min(kMASS_MIN + (ib + 1)*kBIN_WIDTH, kMASS_MAX));
		RooAbsReal * integ = pdf_total_US->createIntegral(mass, NormSet(mass), Range("diag"));
		Double_t mu = integ->getVal()*US_entries;
		delete integ;
		mu_tot += mu;
		if (mu < 5) n_low++;
		if (mu > 0) chi2_pearson += pow(n_obs[ib] - mu, 2)/mu;
		if (mu > 0) chi2_dev += 2*(mu - n_obs[ib] + (n_obs[ib] > 0 ? n_obs[ib]*log(n_obs[ib]/mu) : 0));
		cout << Form("  %5.2f %7.0f %9.2f %6.2f", kMASS_MIN + ib*kBIN_WIDTH, n_obs[ib], mu, mu > 0 ? (n_obs[ib] - mu)/sqrt(mu) : 0) << endl;
	}
	// Parameters pinned at a bound are effectively fixed and do not remove degrees of freedom
	Int_t n_pinned = 0;
	for (int i = 0; i < result->floatParsFinal().getSize(); i++)
	{
		RooRealVar * p = (RooRealVar *)result->floatParsFinal().at(i);
		Double_t range = p->getMax() - p->getMin();
		if ((p->getVal() - p->getMin() < kLIMIT_FRAC*range) || (p->getMax() - p->getVal() < kLIMIT_FRAC*range)) n_pinned++;
	}
	Int_t ndf_diag = nb_diag - (nparam_US - n_pinned);
	Double_t chi2ndf_frame = chi2ndf_US;	// frame chi2 with all floating parameters subtracted, kept for reference
	chi2ndf_US = ndf_diag > 0 ? chi2_pearson/ndf_diag : 0;
	cout << "Diagnostic US chi2: Pearson " << chi2_pearson << ", deviance " << chi2_dev << ", bins " << nb_diag << ", parameters " << nparam_US << ", pinned at a bound " << n_pinned << ", ndf " << ndf_diag << ", chi2/ndf " << chi2ndf_US << " (frame: " << chi2ndf_frame << "), bins with expected < 5: " << n_low << ", expected total " << mu_tot << " of " << US_entries << endl;

	/*
	// pulls
	TCanvas *cnv_pull = new TCanvas("cnv_pull", "cnv_pull");
	cnv_pull->cd();

	RooHist* hpull;
	hpull = frame->pullHist("ds_US", "total_Norm[mass]");

	hpull->Draw();
	*/

	gStyle->SetOptStat(0) ;
    TH2* hcorr = result->correlationHist();

	result->correlationMatrix().Print();
	result->covarianceMatrix().Print();

	/*
	// Optional: Run Minos only on critical parameters
	for (auto pname : {"sig1_yield", "sig2_yield", "sig3_yield", "bgUS_yield"}) {
		RooAbsArg* p = result->floatParsFinal().find(pname);
		if (p) result->minos(*p);
	}
	*/

	// Print yields
	cout << "US entries: " << US_entries << ", LS entries: " << LS_entries << endl;
	cout << "frac: " << yield1.getVal() << ", yield: " << (yield1.getVal()*US_entries) << endl;
	cout << "error frac: " << yield1.getPropagatedError(*result) << ", error: " << (yield1.getPropagatedError(*result)*US_entries) << endl << endl;

	// Print fit results
	cout << "Results!" << endl;
	result->Print("v");

	
	// Signal yields (1S+2S+3S, 1S, 2S, 3S) and fit quality flags
	Double_t N_S = yield1.getVal()*US_entries;
	Double_t dN_S = yield1.getPropagatedError(*result)*US_entries;
	Double_t N_1S = yield_1S.getVal()*US_entries;
	Double_t dN_1S = yield_1S.getPropagatedError(*result)*US_entries;
	Double_t N_2S = yield_2S.getVal()*US_entries;
	Double_t dN_2S = yield_2S.getPropagatedError(*result)*US_entries;
	Double_t N_3S = yield3.getVal()*US_entries;
	Double_t dN_3S = yield3.getPropagatedError(*result)*US_entries;

	// Bin counting: US entries in the window, minus the like-sign entries (combinatorial) and the fitted residual background,
	// divided by the fraction of the fitted signal inside the window (the tails outside are model dependent)
	Double_t c_lo = std::max(kCOUNT_LO, kMASS_MIN);
	Double_t c_hi = std::min(kCOUNT_HI, kMASS_MAX);
	mass.setRange("count", c_lo, c_hi);
	TString s_win = Form("m>=%.6f&&m<%.6f", c_lo, c_hi);
	Double_t N_win = ds_US->sumEntries(s_win);
	Double_t L_win = ds_LS->sumEntries(s_win);
	RooAbsReal * int_res = pdf_bg_res->createIntegral(mass, NormSet(mass), Range("count"));
	RooAbsReal * int_sig1 = pdf_sig1->createIntegral(mass, NormSet(mass), Range("count"));
	RooAbsReal * int_sig2 = pdf_sig2->createIntegral(mass, NormSet(mass), Range("count"));
	RooAbsReal * int_sig3 = pdf_sig3->createIntegral(mass, NormSet(mass), Range("count"));
	RooFormulaVar form_B_res("form_B_res", "residual background fraction in the window", "(1.0-@0)*@1*@2", RooArgList(var_bg_comb_frac, var_bg_res_frac, *int_res));
	RooFormulaVar form_eff("form_eff", "signal fraction in the window", "@0*@2+(1.0-@0)*@1*@3+(1.0-@0)*(1.0-@1)*@4", RooArgList(var_sig1_frac, var_sig2_frac, *int_sig1, *int_sig2, *int_sig3));
	Double_t B_res = form_B_res.getVal()*US_entries;
	Double_t dB_res = form_B_res.getPropagatedError(*result)*US_entries;
	Double_t eff_win = form_eff.getVal();
	Double_t deff_win = form_eff.getPropagatedError(*result);
	Double_t N_raw = N_win - L_win - B_res;
	Double_t dN_raw = sqrt(N_win + L_win + dB_res*dB_res);
	Double_t N_count = N_raw/eff_win;
	Double_t dN_count = sqrt(pow(dN_raw/eff_win, 2) + pow(N_raw*deff_win/(eff_win*eff_win), 2));
	cout << "Counting in [" << c_lo << ", " << c_hi << "]: US " << N_win << ", LS " << L_win << ", residual " << B_res << ", signal fraction " << eff_win << ", count " << N_count << " +- " << dN_count << " (fit " << N_S << ")" << endl;

	// FAILED: fit did not converge or yield is unusable. WARN: suspicious but usable.
	// Lineshape and fraction parameters (var_sig*, except the mean) are bounded by embedding, sitting at those bounds is expected
	TString s_fail = "";
	TString s_warn = "";
	TString s_info = "";	// noted in the result file, does not change the status
	// Minuit status 1 (covariance made positive-definite) with covariance quality 3 is a usable fit
	if (result->status() != 0 && !(result->status() == 1 && result->covQual() == 3)) s_fail += Form(" status=%d", result->status());
	if (result->covQual() != 3) s_fail += Form(" covQual=%d", result->covQual());
	if (!std::isfinite(N_S) || !std::isfinite(dN_S) || dN_S <= 0) s_fail += " yield_not_finite";
	if (US_entries < kMIN_US_ENTRIES) s_fail += Form(" US_entries=%.0f", US_entries);
	if (chi2ndf_US > kCHI2NDF_MAX) s_warn += Form(" chi2ndf_US=%.2f", chi2ndf_US);

	// Save results in a text file
	ofstream out;
	out.open("output/results/" + kTRIG + "_" + kBIN + "_" + s_date + "_" + s_tag + ".txt");
	out << "RUN: " << kTRIG << " " << kBIN << " " << kBG << " " << kSIG << " " << s_tag << endl;
	out << "US_entries: " << US_entries << endl;
	out << "LS_entries: " << LS_entries << endl;
	out << "chi2ndf_US: " << chi2ndf_US << endl;
	out << "chi2ndf_US_frame: " << chi2ndf_frame << endl;
	out << "chi2ndf_LS: " << chi2ndf_LS << endl;
	out << "chi2_diag: " << chi2_pearson << " " << chi2_dev << " " << nb_diag << " " << nparam_US << " " << n_low << " " << n_pinned << endl;
	out << "fit_attempts: " << fit_attempts << endl;
	out << "polish_passes: " << polish_passes << endl;
	out << "kept_start: " << best_start << " of " << kN_STARTS << endl;
	out << "minuit_status: " << result->status() << endl;
	out << "cov_quality: " << result->covQual() << endl;
	out << "edm: " << result->edm() << endl;
	out << "nll: " << Form("%.4f", nll_final) << endl;
	out << "mean_fixed: " << (var_sig1_mu.isConstant() ? mean_int : 0) << endl;
	out << "S: " << N_S << " " << dN_S << endl;
	out << "S_1S: " << N_1S << " " << dN_1S << endl;
	out << "S_2S: " << N_2S << " " << dN_2S << endl;
	out << "S_3S: " << N_3S << " " << dN_3S << endl;
	out << "COUNT_WINDOW: " << c_lo << " " << c_hi << endl;
	out << "COUNT_RAW: " << N_win << " " << L_win << " " << B_res << " " << dB_res << endl;
	out << "COUNT_EFF: " << eff_win << " " << deff_win << endl;
	out << "COUNT: " << N_count << " " << dN_count << endl;
	out << "B: " << (US_entries - N_S) << endl;
	out << "S/sqrt(S+B): " << N_S/sqrt(US_entries) << endl;

	// One line per floating parameter: name, value, error, lower bound, upper bound, flag
	RooArgList pars_final = result->floatParsFinal();
	for (int i = 0; i < pars_final.getSize(); i++)
	{
		RooRealVar * p = (RooRealVar *)pars_final.at(i);
		Double_t range = p->getMax() - p->getMin();
		bool at_limit = (p->getVal() - p->getMin() < kLIMIT_FRAC*range) || (p->getMax() - p->getVal() < kLIMIT_FRAC*range);
		TString s_name = p->GetName();
		bool embedding = s_name.BeginsWith("var_sig") && !s_name.EndsWith("_mu");
		if (at_limit && s_name == "var_sig1_mu") s_info += " " + s_name + "_at_limit";
		else if (at_limit && !embedding) s_warn += " " + s_name + "_at_limit";
		out << "PAR: " << s_name << " " << p->getVal() << " " << p->getError() << " " << p->getMin() << " " << p->getMax() << " " << (at_limit ? "AT_LIMIT" : "OK") << endl;
	}

	if (s_fail.Length() > 0) out << "FIT_STATUS: FAILED:" << s_fail << endl;
	else if (s_warn.Length() > 0) out << "FIT_STATUS: WARN:" << s_warn << endl;
	else out << "FIT_STATUS: OK" << endl;
	if (s_info.Length() > 0) out << "INFO:" << s_info << endl;
	out.close();

	if (s_fail.Length() > 0) cout << "!!! FIT FAILED:" << s_fail << " !!!" << endl;
	else if (s_warn.Length() > 0) cout << "!!! FIT WARNING:" << s_warn << " !!!" << endl;

	TLatex *latex1 = new TLatex();
	latex1->SetTextSize(0.045);
	latex1->SetNDC();
	latex1->DrawLatex(0.55, 0.85, Form("#frac{#chi^{2}}{N^{US}_{dof}} = %.2f, #frac{#chi^{2}}{N^{LS}_{dof}} = %.2f", chi2ndf_US, chi2ndf_LS));
	latex1->DrawLatex(0.55, 0.75, Form("#Upsilon(1S+2S+3S) = %.0f #pm %.0f", (yield1.getVal()*US_entries), (yield1.getPropagatedError(*result)*US_entries)));
	latex1->DrawLatex(0.55, 0.65, Form("S/#sqrt{S+B} = %.2f", N_S/sqrt(US_entries)));
	latex1->DrawLatex(0.55, 0.55, Form("bin counting = %.0f #pm %.0f", N_count, dN_count));

	cout << "S = " << (yield1.getVal())*US_entries << ", B = " << (1- yield1.getVal())*US_entries << endl;
	cout << "S/sqrt(S+B) = " << (yield1.getVal()*US_entries)/ sqrt(US_entries) << endl;

	/*
	TLatex *latex3 = new TLatex();
	latex3->SetTextSize(0.05);
	latex3->SetNDC();
	latex3->DrawLatex(0.55, 0.6, "1 < #LT TofMult #GT < 2");
	*/

	// TODO rewrite as with pT bins!!!
	cnv->cd()->SaveAs("output/fig/" + kTRIG + "_" + kBIN + "_" + s_date + "_" + s_tag + ".png");
	cnv->cd()->SaveAs("output/rootfiles/" + kTRIG + "_" + kBIN + "_" + s_date + "_" + s_tag + ".root");
	
	gStyle->SetPaintTextFormat("4.2f");
	hcorr->Draw("colz text");
	gPad->SetRightMargin(0.15);
	cnv->cd()->SaveAs("output/corr/" + kTRIG + "_" + kBIN + "_" + s_date + "_" + s_tag + "_corr.png");

	// Clean up
	delete cnv;
	delete frame;
	delete result;
}
/////////////////////////////////////////////////////////////