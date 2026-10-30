#include "setup.hxx"
//#include "pars/parsBHT1.hxx"
#include "pars/parsBHT2.hxx"

/////////////////////////////////////////////////////////////
void UpsFitter(TString kTRIG = "BHT2", TString kBIN = "int", TString kBG = "exp", TString kSIG = "CBShape")
{
	// Load the data file
	if (kTRIG == "BHT1")
		f_IN = new TFile(f_BHT1, "READ");
	else if (kTRIG == "BHT2")
		f_IN = new TFile(f_BHT2, "READ");
	else
	{
		cout << "Invalid trigger specified. Exiting..." << endl;
		return 1;
	}

	// Check if the bin makes sense
	for (auto b : vec_bins)
		if (b == kBIN) bool_bin = true;
	if (!bool_bin)
	{
		cout << "Invalid bin specified. Exiting..." << endl;
		return 2;
	}

	// Load datasets for the specified bin
	TTree * tree = (TTree *)f_IN->Get("t_" + kBIN);
		
	RooDataSet * ds_US = new RooDataSet("ds_US_" + kBIN, "ds_US_" + kBIN, RooArgSet(mass, charge), Import(*tree), Cut("charge==0"));
	ds_US = (RooDataSet *)ds_US->reduce(RooArgSet(mass));
	ds_US->Print();
	RooDataSet * ds_LS = new RooDataSet("ds_LS_" + kBIN, "ds_LS_" + kBIN, RooArgSet(mass, charge), Import(*tree), Cut("charge!=0"));
	ds_LS = (RooDataSet *)ds_LS->reduce(RooArgSet(mass));
	ds_LS->Print();

	// This calls the fitting function
	fit_bin(ds_US, ds_LS, kTRIG, kBIN, kBG, kSIG);

	// Clean up
	delete ds_US;
	delete ds_LS;
	f_IN->Close();
	cout << "Fit for " << kTRIG << "_" << kBIN << " completed!" << endl;

	return 0;
}

/////////////////////////////////////////////////////////////
// TODO: rename the parameters sensibly and clean this up
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS, TString kTRIG, TString kBIN, TString kBG, TString kSIG)
{
	// Save number of entries
	double US_entries = ds_US->sumEntries();
	double LS_entries = ds_LS->sumEntries();

	// Define a sample for simultaneous fitting
	RooCategory sample("sample", "sample");
	sample.defineType("US");
	sample.defineType("LS");

	// Simultaneous model
	RooSimultaneous simPdf("simPdf", "simPdf", sample);

	// Define a combined dataset 
	RooDataSet ds_comb("ds_comb", "Combined dataset", RooArgList(mass), Index(sample), Import("US", *ds_US), Import("LS", *ds_LS));
	
	// Get parameters for this bin and trigger
	vector<vector<Double_t>> par_sig = map_par[kTRIG + "_" + kBIN];
	vector<vector<Double_t>> par_bg = map_par_bg[kTRIG + "_" + kBG];
	vector<vector<Double_t>> frac = map_frac[kTRIG + "_" + kBIN];

	// Load parameters and create PDFs
	// Background
	if (kBG == "exp")
	{
		RooRealVar var_bg_comb1("var_bg_comb","Combiantorial background exponent", par_bg[0][0], par_bg[0][1], par_bg[0][2]);
		RooRealVar var_bg_res1("var_bg_res","Residual background exponent", par_bg[1][0], par_bg[1][1], par_bg[1][2]);

		RooExponential pdf_bg_comb("pdf_bg_comb","Background, combinatorial, exponential", mass, var_bg_comb1);
		RooExponential pdf_bg_res("pdf_bg_res","Background, residual, exponential", mass, var_bg_res1);
	}
	else if (kBG == "cheby")
	{
		cout << "Not implemented yet :(. Exitting..." << endl;
		return;
	}
	else
	{
		cout << "Invalid background type. Exitting..." << endl;
		return;
	}

	// Signal
	if (kSIG == "CBShape")
	{
		RooRealVar var_sig1_mu("var_sig1_mu","Signal 1 #mu", par_sig[0][0], par_sig[0][1], par_sig[0][2]);
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

		RooFormulaVar form_sig2_mu_bound("form_sig2_mu", "Signal 2 #mu - bound", "@0 + 0.560180", RooArgList(sig1_mu));
		RooFormulaVar form_sig3_mu_bound("form_sig3_mu", "Signal 3 #mu - bound", "@0 + 0.890502", RooArgList(sig1_mu));

		/*
		/////////////////////
		YOU WERE WARNED ¯\_(ツ)_/¯
		/////////////////////
		*/

		RooCBShape pdf_sig1("pdf_sig1", "Signal 1 - Crystall Ball, onesided", mass, var_sig1_mu, var_sig1_sigma, var_sig1_alpha, var_sig1_n);
		RooCBShape pdf_sig2("pdf_sig2", "Signal 2 - Crystall Ball, onesided", mass, form_sig2_mu_bound, var_sig2_sigma, var_sig2_alpha, var_sig2_n);
		RooCBShape pdf_sig3("pdf_sig3", "Signal 3 - Crystall Ball, onesided", mass, form_sig3_mu_bound, var_sig3_sigma, var_sig3_alpha, var_sig3_n);
	}
	else
	{
		cout << "Invalid signal shape. Exiting..." << endl;
		return;
	}

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
	if (kBG == "exp" && kSIG == "CBShape")
		RooAddPdf pdf_total_US("pdf_total_US", "Total PDF for US data", RooArgList(pdf_bg_comb, pdf_bg_res, pdf_sig1, pdf_sig2, pdf_sig3), RooArgList(var_bg_comb_frac, var_bg_res_frac, var_sig1_frac, var_sig2_frac), kTRUE);
	else
	{
		cout << "Invalid signal/background commbination. Exitting..." << endl;
		return;
	}

	// Add PDFs to the simultaneous model
	simPdf->addPdf(*expLS, "LS");
	simPdf->addPdf(pdf_total_US, "US");

	// Fit the model to the combined dataset
	RooFitResult * result = simPdf->fitTo(*ds_comb,
										Extended(true),
										//Minimizer("Minuit2", "scan"),
										Minimizer("Minuit2", "migrad"),
										MaxCalls(100000),
										Strategy(2),
										EvalErrorWall(1),
										RecoverFromUndefinedRegions(1000),
										Save()
									);

	// Compute yields
	/*
	RooGenericPdf yield1("yield1", "yield from signal fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)", RooArgSet(bgLS_frac, bgUS_frac));
	RooGenericPdf yield2("yield2", "yield from signal 2 fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)*(1.0-sig1_frac)", RooArgSet(bgLS_frac, bgUS_frac, sig1_frac));
	RooGenericPdf yield3("yield3", "yield from signal 3 fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)*(1.0-sig1_frac)*(1.0-sig2_frac)", RooArgSet(bgLS_frac, bgUS_frac, sig1_frac, sig2_frac));
	*/

	cout << "US entries: " << US_entries << ", LS entries: " << LS_entries << endl << endl;
	cout << "yield frac: " << yield1.getVal() << ", yield: " << (yield1.getVal()*US_entries) << endl;
	cout << "error frac: " << yield1.getPropagatedError(*result) << ", error: " << (yield1.getPropagatedError(*result)*US_entries) << endl << endl;
	cout << "yield frac: " << yield1.getVal() << ", yield: " << ((yield1.getVal()-yield2.getVal()-yield3.getVal())*US_entries) << endl;
	cout << "yield2 frac: " << yield2.getVal() << ", yield2: " << ((yield2.getVal()-yield3.getVal()) * US_entries) << endl;
	cout << "yield3 frac: " << yield3.getVal() << ", yield3: " << (yield3.getVal() * US_entries) << endl;

	// Create a canvas for plotting
	TCanvas *cnv = new TCanvas("cnv", "cnv", 1600, 1200);

	// Plot simulated data
	cnv->cd();
	RooPlot* frame = mass.frame(((max_mass - min_mass)/bin_width));
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
	total_sim->plotOn(frame, VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillStyle(3001));
	total_sim->plotOn(frame);
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
	frame->Draw("F");
	//--------------------------------------------

	// Calculate chi2 after extracting the number of parameters - used to calculate Ndof
	Int_t nparam_US = total_sim->getParameters(*ds_US)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_US = " << nparam_US << endl;
	Double_t chi2ndf_US = frame->chiSquare(nparam_US);

	Int_t nparam_LS = expLS->getParameters(*ds_LS)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_LS = " << nparam_LS << endl;
	Double_t chi2ndf_LS = frame->chiSquare(nparam_LS);

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

	//RooGenericPdf yield1("yield1", "yield from signal fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)", RooArgSet(bgLS_frac, bgUS_frac));

	// Compute signal yields
	RooAbsReal * sig1_int = gauss1->createIntegral(mass, NormSet(mass), Range("full"));
	RooProduct sig1_yield("sig1_yield", "sig1_yield", RooArgSet(*sig1_int, US_entries));
	RooAbsReal * sig2_int = gauss2->createIntegral(mass, NormSet(mass), Range("full"));
	RooProduct sig2_yield("sig2_yield", "sig2_yield", RooArgSet(*sig2_int, US_entries));
	RooAbsReal * sig3_int = gauss3->createIntegral(mass, NormSet(mass), Range("full"));
	RooProduct sig3_yield("sig3_yield", "sig3_yield", RooArgSet(*sig3_int, US_entries));


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

	cout << "sig1 yield: " << sig1_yield.getVal() << endl;
	cout << "sig2 yield: " << sig2_yield.getVal() << endl;
	cout << "sig3 yield: " << sig3_yield.getVal() << endl;
	cout << "combined: " << (sig1_yield.getVal() + sig2_yield.getVal() + sig3_yield.getVal()) << endl;

	// Print fit results
	cout << "Results!" << endl;
	result->Print("v");

	// Save results in a text file
	ofstream out;
	out.open("output/results/" + trigger + "_" + bin + "_" + date + "_" + tag + ".txt");
	out << "US entries: " << US_entries << ", LS entries: " << LS_entries << endl;
	out << "#chi^{2}/N^{US}_{dof} = " << chi2ndf_US << ", #chi^{2}/N^{LS}_{dof} = " << chi2ndf_LS << endl;
	out << "frac: " << yield1.getVal() << ", yield: " << (yield1.getVal()*US_entries) << endl;
	out << "error frac: " << yield1.getPropagatedError(*result) << ", error: " << (yield1.getPropagatedError(*result)*US_entries) << endl;
	out << "yield frac: " << yield1.getVal() << ", yield: " << ((yield1.getVal()-yield2.getVal()-yield3.getVal())*US_entries) << endl;
	out << "yield2 frac: " << yield2.getVal() << ", yield2: " << ((yield2.getVal()-yield3.getVal()) * US_entries) << endl;
	out << "yield3 frac: " << yield3.getVal() << ", yield3: " << (yield3.getVal() * US_entries) << endl;
	out << "bgLS: " << bgLS.getVal() << " #pm " << bgLS.getError() << endl;
	out << "bgUS: " << bgUS.getVal() << " #pm " << bgUS.getError() << endl;
	out << "sig1 mu: " << sig1_mu.getVal() << " #pm " << sig1_mu.getError() << endl;
	out << "sig2 mu: " << sig2_mu.getVal() << " #pm " << sig2_mu.getError() << endl;
	out << "sig3 mu: " << sig3_mu.getVal() << " #pm " << sig3_mu.getError() << endl;
	out << "sig1 sigma: " << sig1_sigma.getVal() << " #pm " << sig1_sigma.getError() << endl;
	out << "sig2 sigma: " << sig2_sigma.getVal() << " #pm " << sig2_sigma.getError() << endl;
	out << "sig3 sigma: " << sig3_sigma.getVal() << " #pm " << sig3_sigma.getError() << endl;
	out << "sig1 alpha: " << sig1_alpha.getVal() << " #pm " << sig1_alpha.getError() << endl;
	out << "sig2 alpha: " << sig2_alpha.getVal() << " #pm " << sig2_alpha.getError() << endl;
	out << "sig3 alpha: " << sig3_alpha.getVal() << " #pm " << sig3_alpha.getError() << endl;
	out << "sig1 n: " << sig1_n.getVal() << " #pm " << sig1_n.getError() << endl;
	out << "sig2 n: " << sig2_n.getVal() << " #pm " << sig2_n.getError() << endl;
	out << "sig3 n: " << sig3_n.getVal() << " #pm " << sig3_n.getError() << endl;
	out << "bgLS_frac: " << bgLS_frac.getVal() << " #pm " << bgLS_frac.getError() << endl;
	out << "bgUS_frac: " << bgUS_frac.getVal() << " #pm " << bgUS_frac.getError() << endl;
	out << "sig1_frac: " << sig1_frac.getVal() << " #pm " << sig1_frac.getError() << endl;
	out << "sig2_frac: " << sig2_frac.getVal() << " #pm " << sig2_frac.getError() << endl;
	out << "S = " << (yield1.getVal()*US_entries) << ", B = " << US_entries << endl;
	out << "S/sqrt(S+B) = " << (yield1.getVal()*US_entries)/ sqrt((yield1.getVal()*US_entries + US_entries)) << endl;
	out.close();

	
	TLatex *latex1 = new TLatex();
	latex1->SetTextSize(0.045);
	latex1->SetNDC();
	latex1->DrawLatex(0.55, 0.85, Form("#frac{#chi^{2}}{N^{US}_{dof}} = %.2f, #frac{#chi^{2}}{N^{LS}_{dof}} = %.2f", chi2ndf_US, chi2ndf_LS));
	latex1->DrawLatex(0.55, 0.75, Form("#Upsilon(1S+2S+3S) = %.0f #pm %.0f", (yield1.getVal()*US_entries), (yield1.getPropagatedError(*result)*US_entries)));
	latex1->DrawLatex(0.55, 0.65, Form("S/#sqrt{S+B} = %.2f", (yield1.getVal()*US_entries)/ sqrt((yield1.getVal()*US_entries + US_entries))));

	cout << "S = " << (yield1.getVal()*US_entries) << ", B = " << US_entries << endl;
	cout << "S/sqrt(S+B) = " << (yield1.getVal()*US_entries)/ sqrt((yield1.getVal()*US_entries + US_entries)) << endl;

	/*
	TLatex *latex3 = new TLatex();
	latex3->SetTextSize(0.05);
	latex3->SetNDC();
	latex3->DrawLatex(0.55, 0.6, "1 < #LT TofMult #GT < 2");
	*/

	// TODO rewrite as with pT bins!!!
	cnv->cd()->SaveAs("output/fig/" + trigger + "_" + bin + "_" + date + "_" + tag + ".png");
	cnv->cd()->SaveAs("output/rootfiles/" + trigger + "_" + bin + "_" + date + "_" + tag + ".root");
	
	gStyle->SetPaintTextFormat("4.2f");
	hcorr->Draw("colz text");
	gPad->SetRightMargin(0.15);
	cnv->cd()->SaveAs("output/corr/" + trigger + "_" + bin + "_" + date + "_" + tag + "_corr.png");

	// Clean up
	delete cnv;
	delete frame;
	delete result;
}
/////////////////////////////////////////////////////////////