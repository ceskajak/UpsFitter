#include "setup.hxx"

/////////////////////////////////////////////////////////////
void UpsFitter()
{
	// Load the data file
	f_data = new TFile(infile);

	// Load datasets for the specified bin
	TTree * tree = (TTree *)f_data->Get("tree");
		
	RooDataSet * ds_US = new RooDataSet("ds_US", "ds_US", RooArgSet(mass, charge), Import(*tree), Cut("charge==0"));
	ds_US = (RooDataSet *)ds_US->reduce(RooArgSet(mass));
	ds_US->Print();
	RooDataSet * ds_LS = new RooDataSet("ds_LS", "ds_LS", RooArgSet(mass, charge), Import(*tree), Cut("charge!=0"));
	ds_LS = (RooDataSet *)ds_LS->reduce(RooArgSet(mass));
	ds_LS->Print();

	// This calls the fitting function
	fit_bin(ds_US, ds_LS);

	// Clean up
	delete ds_US;
	delete ds_LS;
	f_data->Close();
	cout << "Fit completed!" << endl;
}

/////////////////////////////////////////////////////////////
// TODO: rename the parameters sensibly and clean this up
void fit_bin(RooDataSet * ds_US, RooDataSet * ds_LS)
{
	// Save number of entries
	double US_entries = ds_US->sumEntries();
	double LS_entries = ds_LS->sumEntries();

	// Define a sample for simultaneous fitting
	RooCategory sample("sample", "sample");
	sample.defineType("US");
	sample.defineType("LS");

	// Define a combined dataset 
	RooDataSet * ds_comb = new RooDataSet("ds_comb", "Combined dataset", RooArgList(mass), Index(sample), Import("US", *ds_US), Import("LS", *ds_LS));

	// Define parameters
	RooRealVar bgComb("bgComb","Comb", 0, -1, 1);
	RooRealVar bgRes("bgRes","Res", 0, -1, 1);
	RooRealVar sig_mu("sig_mu","#mu", 10, 5, 15);
	RooRealVar sig_sigma("sig_sigma","#sigma", 1, 0, 5);
	RooRealVar sig_alpha("sig_alpha","#alpha", 1, 0, 10);
	RooRealVar sig_n("sig_n","#n", 1, 0, 10);

	// Define recursive fractions for the components
	RooRealVar bgComb_frac("bgComb_frac","bgComb_frac", 0.5, 0., 1.);
	RooRealVar bgRes_frac("bgRes_frac","bgRes_frac", 0.5, 0., 1.);

	// Ensure the combinatorial background describes the LS data
	bgComb_frac.setVal(ds_LS->sumEntries()/ds_US->sumEntries());
	bgComb_frac.setConstant(true);

	// Create PDFs
	RooExponential *expComb = new RooExponential("expComb","Comb exp background", mass, bgComb);
	RooExponential *expRes = new RooExponential("expRes","Res exp background", mass, bgRes);
	RooCBShape *sig = new RooCBShape("gauss1", "CB signal 1", mass, sig_mu, sig_sigma, sig_alpha, sig_n);

	// Add PDFs together with recursive fractions (forces all ratios to be in [0,1])
	RooAddPdf *total_sim = new RooAddPdf("total_sim", "US total", RooArgList(*expComb, *expRes, *sig), RooArgList(bgComb_frac, bgRes_frac), kTRUE);

	// Simultaneous model
	RooSimultaneous * simPdf = new RooSimultaneous("simPdf", "simPdf", sample);

	// Add PDFs to the simultaneous model
	simPdf->addPdf(*expComb, "LS");
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
	total_sim->plotOn(frame, Components(*expComb), LineColor(kRed), LineStyle(kDashed), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kRed), FillStyle(3001));
	total_sim->plotOn(frame, Components(*expComb), LineColor(kRed), LineStyle(kDashed));
	total_sim->plotOn(frame, Components(*expRes), LineColor(kGray), LineStyle(kDashed), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kGray), FillStyle(3001));
	total_sim->plotOn(frame, Components(*expRes), LineColor(kGray), LineStyle(kDashed));
	total_sim->plotOn(frame, Components(*sig), LineColor(kMagenta), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kMagenta), FillStyle(3001));
	total_sim->plotOn(frame, Components(*sig), LineColor(kMagenta), LineStyle(kSolid));
	total_sim->plotOn(frame, Components(RooArgSet(*expComb, *expRes)), LineColor(kOrange), LineStyle(kSolid), VisualizeError(*result, sigma, useLinear), DrawOption("F"), FillColor(kOrange), FillStyle(3001));
	total_sim->plotOn(frame, Components(RooArgSet(*expComb, *expRes)), LineColor(kOrange), LineStyle(kSolid));
	total_sim->plotOn(frame);
	frame->Draw("F");
	//--------------------------------------------

	// Compute yields
	RooGenericPdf yield("yield", "yield from signal fit", "(1.0-bgComb_frac)*(1.0-bgRes_frac)", RooArgSet(bgComb_frac, bgRes_frac));
	cout << "US entries: " << US_entries << ", LS entries: " << LS_entries << endl;
	cout << "yield frac: " << yield.getVal() << ", yield: " << (yield.getVal()*US_entries) << endl;
	cout << "error frac: " << yield.getPropagatedError(*result) << ", error: " << (yield.getPropagatedError(*result)*US_entries) << endl;

	// Calculate chi2 after extracting the number of parameters - used to calculate Ndof
	Int_t nparam_US = total_sim->getParameters(*ds_US)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_US = " << nparam_US << endl;
	Double_t chi2ndf_US = frame->chiSquare(nparam_US);
	cout << "chi2/ndf (US) = " << chi2ndf_US << endl;

	Int_t nparam_LS = expComb->getParameters(*ds_LS)->selectByAttrib("Constant", kFALSE)->getSize();
	cout << "nparam_LS = " << nparam_LS << endl;
	Double_t chi2ndf_LS = frame->chiSquare(nparam_LS);
	cout << "chi2/ndf (LS) = " << chi2ndf_LS << endl;

	/*
	// pulls
	TCanvas *cnv_pull = new TCanvas("cnv_pull", "cnv_pull");
	cnv_pull->cd();

	RooHist* hpull;
	hpull = frame->pullHist("ds_US", "total_Norm[mass]");

	hpull->Draw();
	*/

	/*
	gStyle->SetOptStat(0) ;
    TH2* hcorr = result->correlationHist();

	result->correlationMatrix().Print();
	result->covarianceMatrix().Print();

	
	TLatex *latex1 = new TLatex();
	latex1->SetTextSize(0.045);
	latex1->SetNDC();
	latex1->DrawLatex(0.55, 0.85, Form("#frac{#chi^{2}}{N^{US}_{dof}} = %.2f, #frac{#chi^{2}}{N^{LS}_{dof}} = %.2f", chi2ndf_US, chi2ndf_LS));
	latex1->DrawLatex(0.55, 0.75, Form("#Upsilon(1S+2S+3S) = %.0f #pm %.0f", (yield1.getVal()*US_entries), (yield1.getPropagatedError(*result)*US_entries)));
	latex1->DrawLatex(0.55, 0.65, Form("S/#sqrt{S+B} = %.2f", (yield1.getVal()*US_entries)/ sqrt((yield1.getVal()*US_entries + US_entries))));
	*/

	cout << "S = " << (yield.getVal()*US_entries) << ", S+B = " << US_entries << endl;
	cout << "S/sqrt(S+B) = " << (yield.getVal()*US_entries)/ sqrt(US_entries) << endl;

	/*
	TLatex *latex3 = new TLatex();
	latex3->SetTextSize(0.05);
	latex3->SetNDC();
	latex3->DrawLatex(0.55, 0.6, "1 < #LT TofMult #GT < 2");
	*/

	// TODO rewrite as with pT bins!!!
	cnv->cd()->SaveAs("output.png");

	// Clean up
	delete cnv;
	delete frame;
	delete result;
}
/////////////////////////////////////////////////////////////