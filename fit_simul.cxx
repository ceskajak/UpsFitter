#include "setup.hxx"
#include "pars/parsBHT1.hxx"
#include "pars/parsBHT2.hxx"

/////////////////////////////////////////////////////////////
void fit_simul(TString trigger = "BHT2", TString bin = "9mult20")
{
	// Load the data file
	if (trigger == "BHT1")
		f_data = new TFile(infile_BHT1, "READ");
	else if (trigger == "BHT2")
		f_data = new TFile(infile_BHT2, "READ");
	else
	{
		cout << "Invalid trigger specified. Exiting." << endl;
		return;
	}

	// Load datasets for the specified bin
	TTree * tree = (TTree *)f_data->Get("t_" + bin);
		
	RooDataSet * ds_US = new RooDataSet("ds_US_" + bin, "ds_US_" + bin, RooArgSet(mass, charge), Import(*tree), Cut("charge==0"));
	ds_US = (RooDataSet *)ds_US->reduce(RooArgSet(mass));
	ds_US->Print();
	RooDataSet * ds_LS = new RooDataSet("ds_LS_" + bin, "ds_LS_" + bin, RooArgSet(mass, charge), Import(*tree), Cut("charge!=0"));
	ds_LS = (RooDataSet *)ds_LS->reduce(RooArgSet(mass));
	ds_LS->Print();

	// This calls the fitting function
	fit_bin(ds_US, ds_LS, trigger, bin);

	// Clean up
	delete ds_US;
	delete ds_LS;
	f_data->Close();
	cout << "Fit for " << trigger << "_" << bin << " completed!" << endl;
}

/////////////////////////////////////////////////////////////
// TODO: rename the parameters sensibly and clean this up
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS, TString trigger, TString bin)
{
	// Define variable + ranges
	RooRealVar mass("mass", "mass", min_mass, max_mass, "GeV/c^{2}");
	mass.setRange("full", min_mass, max_mass);

	// Save number of entries
	double US_entries = ds_US->sumEntries();
	double LS_entries = ds_LS->sumEntries();

	// Define a sample for simultaneous fitting
	RooCategory sample("sample", "sample");
	sample.defineType("US");
	sample.defineType("LS");

	// Define a combined dataset 
	RooDataSet * ds_comb = new RooDataSet("ds_comb", "Combined dataset", RooArgList(mass), Index(sample), Import("US", *ds_US), Import("LS", *ds_LS));

	TString trigger_bin = trigger + "_" + bin;

	vector<vector<Double_t>> par = map_par[trigger_bin];
	vector<vector<Double_t>> frac = map_frac[trigger_bin];

	// Load the parameters from the array (saved in the header)
	RooRealVar bgLS("bgLS","LS", par[0][0], par[0][1], par[0][2]);
	RooRealVar bgUS("bgUS","US", par[1][0], par[1][1], par[1][2]);
	RooRealVar sig1_mu("sig1_mu","#mu", par[2][0], par[2][1], par[2][2]);
	RooRealVar sig2_mu("sig2_mu","#mu", par[3][0], par[3][1], par[3][2]);
	RooRealVar sig3_mu("sig3_mu","#mu", par[4][0], par[4][1], par[4][2]);
	RooRealVar sig1_sigma("sig1_sigma","#sigma", par[5][0], par[5][1], par[5][2]);
	RooRealVar sig2_sigma("sig2_sigma","#sigma", par[6][0], par[6][1], par[6][2]);
	RooRealVar sig3_sigma("sig3_sigma","#sigma", par[7][0], par[7][1], par[7][2]);
	RooRealVar sig1_alpha("sig1_alpha","#alpha", par[8][0], par[8][1], par[8][2]);
	RooRealVar sig2_alpha("sig2_alpha","#alpha", par[9][0], par[9][1], par[9][2]);
	RooRealVar sig3_alpha("sig3_alpha","#alpha", par[10][0], par[10][1], par[10][2]);
	RooRealVar sig1_n("sig1_n","#n", par[11][0], par[11][1], par[11][2]);
	RooRealVar sig2_n("sig2_n","#n", par[12][0], par[12][1], par[12][2]);
	RooRealVar sig3_n("sig3_n","#n", par[13][0], par[13][1], par[13][2]);

	// used for 2D CB
	/*
	RooRealVar sig1_alphaL("sig1_alphaL","#alpha_{L}", par[10][0], par[10][1], par[10][2]);
	RooRealVar sig2_alphaL("sig2_alphaL","#alpha_{L}", par[11][0], par[11][1], par[11][2]);
	RooRealVar sig3_alphaL("sig3_alphaL","#alpha_{L}", par[12][0], par[12][1], par[12][2]);
	RooRealVar sig1_alphaR("sig1_alphaR","#alpha_{R}", par[13][0], par[13][1], par[13][2]);
	RooRealVar sig2_alphaR("sig2_alphaR","#alpha_{R}", par[14][0], par[14][1], par[14][2]);
	RooRealVar sig3_alphaR("sig3_alphaR","#alpha_{R}", par[15][0], par[15][1], par[15][2]);
	RooRealVar sig1_nL("sig1_nL","#n_{L}", par[16][0], par[16][1], par[16][2]);
	RooRealVar sig2_nL("sig2_nL","#n_{L}", par[17][0], par[17][1], par[17][2]);
	RooRealVar sig3_nL("sig3_nL","#n_{L}", par[18][0], par[18][1], par[18][2]);
	RooRealVar sig1_nR("sig1_nR","#n_{R}", par[19][0], par[19][1], par[19][2]);
	RooRealVar sig2_nR("sig2_nR","#n_{R}", par[20][0], par[20][1], par[20][2]);
	RooRealVar sig3_nR("sig3_nR","#n_{R}", par[21][0], par[21][1], par[21][2]);
	*/

	// Load recursive fractions for the components
	RooRealVar bgLS_frac("bgLS_frac","bgLS_frac", frac[0][0], frac[0][1], frac[0][2]);
	RooRealVar bgUS_frac("bgUS_frac","bgUS_frac", frac[1][0], frac[1][1], frac[1][2]);
	RooRealVar sig1_frac("sig1_frac","sig1_frac", frac[2][0], frac[2][1], frac[2][2]);
	RooRealVar sig2_frac("sig2_frac","sig2_frac", frac[3][0], frac[3][1], frac[3][2]);

	// Ensure the combinatorial background describes the LS data
	bgLS_frac.setVal(ds_LS->sumEntries()/ds_US->sumEntries());
	bgLS_frac.setConstant(true);

	// Mass delta from embedding:
	RooFormulaVar sig2_mu_bound("sig2_mu", "sig2_mu", "@0 + 0.560180", RooArgList(sig1_mu));
	RooFormulaVar sig3_mu_bound("sig3_mu", "sig3_mu", "@0 + 0.890502", RooArgList(sig1_mu));

	// Create PDFs
	RooExponential *expLS = new RooExponential("expLS","LS exp background", mass, bgLS);
	RooExponential *expUS = new RooExponential("expUS","US Exp background", mass, bgUS);
	RooCBShape *gauss1 = new RooCBShape("gauss1", "Gaussian signal 1", mass, sig1_mu, sig1_sigma, sig1_alpha, sig1_n);
	RooCBShape *gauss2 = new RooCBShape("gauss2", "Gaussian signal 2", mass, sig2_mu_bound, sig2_sigma, sig2_alpha, sig2_n);
	RooCBShape *gauss3 = new RooCBShape("gauss3", "Gaussian signal 3", mass, sig3_mu_bound, sig3_sigma, sig3_alpha, sig3_n);

	/*
	RooCrystalBall *gauss1 = new RooCrystalBall("gauss1", "Gaussian signal 1", mass, sig1_mu, sig1_sigma, sig1_alphaL, sig1_alphaR, sig1_nL, sig1_nR);
	RooCrystalBall *gauss2 = new RooCrystalBall("gauss2", "Gaussian signal 2", mass, sig2_mu, sig2_sigma, sig2_alphaL, sig2_alphaR, sig2_nL, sig2_nR);
	RooCrystalBall *gauss3 = new RooCrystalBall("gauss3", "Gaussian signal 3", mass, sig3_mu, sig3_sigma, sig3_alphaL, sig3_alphaR, sig3_nL, sig3_nR);
	*/

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	RooAddPdf *total_sim = new RooAddPdf("total_sim", "US total", RooArgList(*expLS, *expUS, *gauss1, *gauss2, *gauss3), RooArgList(bgLS_frac, bgUS_frac, sig1_frac, sig2_frac), kTRUE);

	// Simultaneous model
	RooSimultaneous * simPdf = new RooSimultaneous("simPdf", "simPdf", sample);

	// Add PDFs to the simultaneous model
	simPdf->addPdf(*expLS, "LS");
	simPdf->addPdf(*total_sim, "US");

	// Fit the model to the combined dataset
	RooFitResult *result = simPdf->fitTo(*ds_comb,
										Extended(true),
										//Minimizer("Minuit2", "scan"),
										Minimizer("Minuit2", "migrad"),
										MaxCalls(100000),
										Strategy(2),
										EvalErrorWall(1),
										RecoverFromUndefinedRegions(1000),
										Save());

	// Compute yields
	RooGenericPdf yield1("yield1", "yield from signal fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)", RooArgSet(bgLS_frac, bgUS_frac));
	RooGenericPdf yield2("yield2", "yield from signal 2 fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)*(1.0-sig1_frac)", RooArgSet(bgLS_frac, bgUS_frac, sig1_frac));
	RooGenericPdf yield3("yield3", "yield from signal 3 fit", "(1.0-bgLS_frac)*(1.0-bgUS_frac)*(1.0-sig1_frac)*(1.0-sig2_frac)", RooArgSet(bgLS_frac, bgUS_frac, sig1_frac, sig2_frac));
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
	ds_LS->plotOn(frame, MarkerColor(kRed));
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
	out << "sig1 yield: " << sig1_yield.getVal() << endl;
	out << "sig2 yield: " << sig2_yield.getVal() << endl;
	out << "sig3 yield: " << sig3_yield.getVal() << endl;
	out << "combined: " << (sig1_yield.getVal() + sig2_yield.getVal() + sig3_yield.getVal()) << endl;
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
	out.close();

	
	TLatex *latex1 = new TLatex();
	latex1->SetTextSize(0.05);
	latex1->SetNDC();
	latex1->DrawLatex(0.55, 0.85, Form("#frac{#chi^{2}}{N^{US}_{dof}} = %.2f, #frac{#chi^{2}}{N^{LS}_{dof}} = %.2f", chi2ndf_US, chi2ndf_LS));
	latex1->DrawLatex(0.55, 0.70, Form("#Upsilon(1S+2S+3S) = %.0f #pm %.0f", (yield1.getVal()*US_entries), (yield1.getPropagatedError(*result)*US_entries)));

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