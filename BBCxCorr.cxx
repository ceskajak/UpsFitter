// Derives the multiplicative BBCx (luminosity) correction of TofMult and compares it with the current additive one
// Usage: root -l -b -q 'BBCxCorr.cxx()'
// Input:  minimum bias profile p_TofMultxBBC (kBBC_MB), Upsilon candidate tree (f_BHT2)
// Output: output/bbcx/bbcx_corr.txt (coefficients, slopes, bin populations), output/bbcx/*.png and *.pdf

#include "headers/common.hxx"
#include "headers/setup.hxx"
#include "headers/cuts.hxx"

const TString kBBC_MB = "../../data/BBC.histos_250228.root";
const double kFIT_LO = 1.8, kFIT_HI = 4.6;	// fit window of the MB profile (x = BBCx / 1e6), below ~1.75 the profile has a spike of 1.5 mean
const double kCAL_LO = 2.3, kCAL_HI = 4.7;	// fit window of the calibration on the candidates
const double kX_REF = 3.55;	// reference BBCx: median of the Upsilon candidates (inside the data, not extrapolated to 0)
const int kNEBIN = 12;	// BBCx bins of the flatness check
const double kXMIN = 2.0, kXMAX = 5.0;

struct Slope { double a, ea, mean; };

static Slope SlopeOf(TProfile * p)
{
	TF1 f("f_lin", "pol1", kXMIN, kXMAX);
	p->Fit(&f, "QN0", "", kXMIN, kXMAX);
	return {f.GetParameter(1), f.GetParError(1), p->GetMean(2)};
}

void BBCxCorr()
{
	gSystem->mkdir("output/bbcx", true);
	gStyle->SetOptStat(0);

	// ---- coefficients from the MB profile: TofMult(x) ~ g(x) = g0 + g1 x + g2 x^2, c(x) = g(xref) / g(x)
	TFile * f_mb = TFile::Open(kBBC_MB, "READ");
	TProfile * p_mb = (TProfile *)f_mb->Get("p_TofMultxBBC");
	TProfile * p_mb_old = (TProfile *)f_mb->Get("p_TofMultxBBC_corr");
	TProfile * p_rb = (TProfile *)p_mb->Rebin(10, "p_mb_rb");
	TF1 * f_g = new TF1("f_g", "pol2", kFIT_LO, kFIT_HI);
	p_rb->Fit(f_g, "QN0", "", kFIT_LO, kFIT_HI);
	double g0 = f_g->GetParameter(0), g1 = f_g->GetParameter(1), g2 = f_g->GetParameter(2);
	double gref = f_g->Eval(kX_REF);

	// pol1 and pol3 for the shape dependence of gref
	TF1 f_g1("f_g1", "pol1", kFIT_LO, kFIT_HI), f_g3("f_g3", "pol3", kFIT_LO, kFIT_HI);
	p_rb->Fit(&f_g1, "QN0", "", kFIT_LO, kFIT_HI);
	p_rb->Fit(&f_g3, "QN0", "", kFIT_LO, kFIT_HI);

	f_IN = new TFile(f_BHT2, "READ");
	TTree * tree = (TTree *)f_IN->Get("tree");
	TString cuts = MakeCuts("BHT2") + " && m > 8.4 && m < 11.2";

	// calibration on the candidates themselves: h(x) linear fit of the daughter-subtracted TofMult profile, c(x) = h(xref) / h(x)
	TProfile * p_cal = new TProfile("p_cal", "", 24, kXMIN, kXMAX);
	tree->Draw(s_mult_sub + ":BBCx/1e6>>p_cal", "(" + cuts + ")", "goff prof");
	TF1 f_h("f_h", "pol1", kCAL_LO, kCAL_HI);
	p_cal->Fit(&f_h, "QN0", "", kCAL_LO, kCAL_HI);
	double h0 = f_h.GetParameter(0), h1 = f_h.GetParameter(1), href = f_h.Eval(kX_REF);

	// overall scale: the reference point is a convention, fixed so that the candidate mean equals that of the current correction
	// (the comparison then isolates the removal of the BBCx trend from a rescaling of all bin edges)
	TH1D h_m0("h_m0", "", 600, 0, 60), h_m1("h_m1", "", 600, 0, 60);
	tree->Draw(s_mult_subcorr + ">>h_m0", "(" + cuts + ")", "goff");
	tree->Draw(Form("%s * %.6f / (%.6f + %.6f*BBCx/1e6) >> h_m1", s_mult_sub.Data(), href, h0, h1), "(" + cuts + ")", "goff");
	double mean_cur = h_m0.GetMean(), mean_unit = h_m1.GetMean();
	href *= mean_cur / mean_unit;

	// c(x) as tree expression: MB-derived (diagnostic) and candidate-calibrated (adopted)
	TString s_cmb = Form("(%.6f / (%.6f + %.6f*BBCx/1e6 + %.6f*pow(BBCx/1e6,2)))", gref, g0, g1, g2);
	TString s_c = Form("(%.6f / (%.6f + %.6f*BBCx/1e6))", href, h0, h1);
	TString s_sub = "(TofMult - Sum$(TofMult_id == id1 || TofMult_id == id2))";
	TString s_subcorrmb = Form("(%s * %s)", s_sub.Data(), s_cmb.Data());
	TString s_subcorrm = Form("(%s * %s)", s_sub.Data(), s_c.Data());
	TString s_corrm = Form("(TofMult * %s)", s_c.Data());

	ofstream out("output/bbcx/bbcx_corr.txt");
	out << "MB profile fit " << kFIT_LO << "-" << kFIT_HI << " (rebin 10): g(x) = " << g0 << " + " << g1 << " x + " << g2 << " x^2, chi2/ndf "
		<< f_g->GetChisquare() << "/" << f_g->GetNDF() << "\n";
	out << "pol1 fit: " << f_g1.GetParameter(0) << " + " << f_g1.GetParameter(1) << " x ; pol3 g(xref)=" << f_g3.Eval(kX_REF)
		<< " pol2 g(xref)=" << gref << " pol1 g(xref)=" << f_g1.Eval(kX_REF) << "\n";
	out << "xref " << kX_REF << "  g(xref) " << gref << "\n";
	out << Form("CAL candidates h(x) = %.4f + %.4f x  (fit %.1f-%.1f)  scale A = %.4f (mean matched to current %.4f)  relative slope %.3f /unit\n", h0, h1, kCAL_LO, kCAL_HI, href, mean_cur, h1 / (h0 + h1 * kX_REF));
	out << "EXPR_c_mb " << s_cmb << "\n";
	out << "EXPR_c " << s_c << "\n";
	out << "EXPR_subcorrm " << s_subcorrm << "\n";
	out << "EXPR_corrm " << s_corrm << "\n";

	// ---- MB flatness: old additive correction (file) vs the new multiplicative one on the profile
	TH1D * h_mb_new = new TH1D("h_mb_new", "", kNEBIN, kXMIN, kXMAX);
	TProfile * p_mb_12 = new TProfile("p_mb_12", "", kNEBIN, kXMIN, kXMAX);
	TProfile * p_mb_old_12 = new TProfile("p_mb_old_12", "", kNEBIN, kXMIN, kXMAX);
	TProfile * p_mb_new_12 = new TProfile("p_mb_new_12", "", kNEBIN, kXMIN, kXMAX);
	for (int i = 1; i <= p_mb->GetNbinsX(); i++)
	{
		double n = p_mb->GetBinEntries(i);
		if (n <= 0) continue;
		double x = p_mb->GetBinCenter(i), y = p_mb->GetBinContent(i);
		double c = gref / f_g->Eval(x);
		p_mb_12->Fill(x, y, n);
		p_mb_old_12->Fill(x, p_mb_old->GetBinContent(i), n);
		p_mb_new_12->Fill(x, y * c, n);
	}
	Slope s_mb_raw = SlopeOf(p_mb_12), s_mb_old = SlopeOf(p_mb_old_12), s_mb_new = SlopeOf(p_mb_new_12);
	out << Form("MB slope (x %.1f-%.1f)  raw %.3f +- %.3f (mean %.3f) | file additive %.3f +- %.3f (mean %.3f) | multiplicative %.3f +- %.3f (mean %.3f)\n",
		kXMIN, kXMAX, s_mb_raw.a, s_mb_raw.ea, s_mb_raw.mean, s_mb_old.a, s_mb_old.ea, s_mb_old.mean, s_mb_new.a, s_mb_new.ea, s_mb_new.mean);

	// ---- Upsilon candidates
	struct Var { TString name, expr; TString label; int col; };
	vector<Var> vars = {
		{"sub", s_mult_sub, "TofMult - daughters (no BBCx correction)", kGray + 2},
		{"subcorr", s_mult_subcorr, "current: additive, ref. BBCx = 0", kBlue + 1},
		{"subcorrmb", s_subcorrmb, "multiplicative, from MB profile", kGreen + 2},
		{"subcorrm", s_subcorrm, "new: multiplicative, from candidates", kRed + 1}};

	TCanvas * c1 = new TCanvas("c1", "", 900, 650);
	TLegend * lg = new TLegend(0.40, 0.68, 0.89, 0.89);
	lg->SetBorderSize(0);
	lg->SetFillStyle(0);
	vector<TProfile *> profs;
	vector<Slope> slopes;
	double ymin = 1e9, ymax = -1e9;
	for (auto & v : vars)
	{
		TProfile * p = new TProfile("p_" + v.name, "", kNEBIN, kXMIN, kXMAX);
		tree->Draw(v.expr + ":BBCx/1e6>>p_" + v.name, "(" + cuts + ")", "goff prof");
		profs.push_back(p);
		slopes.push_back(SlopeOf(p));
		for (int i = 1; i <= kNEBIN; i++)
			if (p->GetBinEntries(i) > 30)
			{
				ymin = min(ymin, p->GetBinContent(i) - p->GetBinError(i));
				ymax = max(ymax, p->GetBinContent(i) + p->GetBinError(i));
			}
	}
	for (size_t k = 0; k < vars.size(); k++)
	{
		TProfile * p = profs[k];
		p->SetLineColor(vars[k].col);
		p->SetMarkerColor(vars[k].col);
		p->SetMarkerStyle(20 + k);
		p->SetMinimum(ymin - 0.3);
		p->SetMaximum(ymax + 0.3);
		p->SetTitle(";BBCx / 10^{6};#LTTofMult#GT (#Upsilon candidates, 8.4 < m < 11.2 GeV)");
		p->Draw(k == 0 ? "E1" : "E1 same");
		lg->AddEntry(p, Form("%s: slope %.2f #pm %.2f", vars[k].label.Data(), slopes[k].a, slopes[k].ea), "lp");
		out << Form("UPS slope %-9s %.3f +- %.3f  mean %.3f  relative slope %.3f /unit\n", vars[k].name.Data(), slopes[k].a, slopes[k].ea, slopes[k].mean, slopes[k].a / slopes[k].mean);
	}
	lg->Draw();
	c1->SaveAs("output/bbcx/ups_profiles.png");
	c1->SaveAs("output/bbcx/ups_profiles.pdf");

	// ---- MB profile figure
	TCanvas * c2 = new TCanvas("c2", "", 900, 650);
	p_mb_12->SetLineColor(kGray + 2);
	p_mb_12->SetMarkerColor(kGray + 2);
	p_mb_12->SetMarkerStyle(20);
	p_mb_old_12->SetLineColor(kBlue + 1);
	p_mb_old_12->SetMarkerColor(kBlue + 1);
	p_mb_old_12->SetMarkerStyle(21);
	p_mb_new_12->SetLineColor(kRed + 1);
	p_mb_new_12->SetMarkerColor(kRed + 1);
	p_mb_new_12->SetMarkerStyle(22);
	p_mb_12->SetMinimum(1.0);
	p_mb_12->SetMaximum(3.6);
	p_mb_12->SetTitle(";BBCx / 10^{6};#LTTofMult#GT (minimum bias)");
	p_mb_12->Draw("E1");
	p_mb_old_12->Draw("E1 same");
	p_mb_new_12->Draw("E1 same");
	TLegend * lg2 = new TLegend(0.40, 0.68, 0.89, 0.89);
	lg2->SetBorderSize(0);
	lg2->SetFillStyle(0);
	lg2->AddEntry(p_mb_12, Form("raw: slope %.3f", s_mb_raw.a), "lp");
	lg2->AddEntry(p_mb_old_12, Form("file (additive, ref. BBCx = 0): %.3f", s_mb_old.a), "lp");
	lg2->AddEntry(p_mb_new_12, Form("new (multiplicative, ref. 3.55): %.3f", s_mb_new.a), "lp");
	lg2->Draw();
	c2->SaveAs("output/bbcx/mb_profiles.png");
	c2->SaveAs("output/bbcx/mb_profiles.pdf");

	// ---- TofMult distribution and bin populations (Upsilon candidates), nominal 5-bin and 4-bin edges
	TH1D * h_cur = new TH1D("h_cur", "", 40, 0, 20);
	TH1D * h_new = new TH1D("h_new", "", 40, 0, 20);
	tree->Draw(s_mult_subcorr + ">>h_cur", "(" + cuts + ")", "goff");
	tree->Draw(s_subcorrm + ">>h_new", "(" + cuts + ")", "goff");
	out << Form("UPS TofMult mean/N: current %.3f / %.0f , new %.3f / %.0f\n", h_cur->GetMean(), h_cur->GetEntries(), h_new->GetMean(), h_new->GetEntries());
	vector<pair<double, double>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 8}, {3, 8}};
	for (auto e : edges)
	{
		TString sel = "(" + cuts + ")";
		double nc = tree->GetEntries(sel + Form(" && (%s > 3*%g) && (%s <= 3*%g)", s_mult_subcorr.Data(), e.first, s_mult_subcorr.Data(), e.second));
		double nn = tree->GetEntries(sel + Form(" && (%s > 3*%g) && (%s <= 3*%g)", s_subcorrm.Data(), e.first, s_subcorrm.Data(), e.second));
		out << Form("BINPOP %g-%g  current %.0f  new %.0f  (%+.1f%%)\n", e.first, e.second, nc, nn, 100 * (nn - nc) / max(nc, 1.));
	}

	h_cur->SetLineColor(kBlue + 1);
	h_new->SetLineColor(kRed + 1);
	h_cur->SetLineWidth(2);
	h_new->SetLineWidth(2);
	TCanvas * c3 = new TCanvas("c3", "", 900, 650);
	h_new->SetTitle(";TofMult (daughters subtracted, BBCx corrected);#Upsilon candidates");
	h_new->SetMaximum(1.2 * max(h_cur->GetMaximum(), h_new->GetMaximum()));
	h_new->Draw("HIST");
	h_cur->Draw("HIST same");
	TLegend * lg3 = new TLegend(0.55, 0.72, 0.89, 0.89);
	lg3->SetBorderSize(0);
	lg3->SetFillStyle(0);
	lg3->AddEntry(h_cur, Form("current, mean %.2f", h_cur->GetMean()), "l");
	lg3->AddEntry(h_new, Form("new, mean %.2f", h_new->GetMean()), "l");
	lg3->Draw();
	c3->SaveAs("output/bbcx/ups_tofmult.png");
	c3->SaveAs("output/bbcx/ups_tofmult.pdf");

	// ---- correction factor vs BBCx (the additive one as its equivalent factor at <TofMult> = 6.8)
	TCanvas * c4 = new TCanvas("c4", "", 900, 650);
	TF1 * f_cm = new TF1("f_cm", Form("%.6f / (%.6f + %.6f*x)", href, h0, h1), kXMIN, kXMAX);
	f_cm->SetLineColor(kRed + 1);
	f_cm->SetLineWidth(3);
	f_cm->SetTitle(";BBCx / 10^{6};TofMult correction factor");
	TF1 * f_cb = new TF1("f_cb", Form("%.6f / (%.6f + %.6f*x + %.6f*x*x)", gref, g0, g1, g2), kXMIN, kXMAX);
	f_cb->SetLineColor(kGreen + 2);
	f_cb->SetLineWidth(3);
	TF1 * f_ca = new TF1("f_ca", "1 - (x*0.811027 + x*x*(-0.2533)) / 6.8", kXMIN, kXMAX);
	f_ca->SetLineColor(kBlue + 1);
	f_ca->SetLineWidth(3);
	f_cm->SetMinimum(0.7);
	f_cm->SetMaximum(1.7);
	f_cm->Draw();
	f_cb->Draw("same");
	f_ca->Draw("same");
	TLegend * lg4 = new TLegend(0.15, 0.68, 0.65, 0.89);
	lg4->SetBorderSize(0);
	lg4->SetFillStyle(0);
	lg4->AddEntry(f_ca, "current (additive, equivalent factor at 6.8)", "l");
	lg4->AddEntry(f_cb, "multiplicative from MB profile", "l");
	lg4->AddEntry(f_cm, "multiplicative from candidates (new)", "l");
	lg4->Draw();
	c4->SaveAs("output/bbcx/factor.png");
	c4->SaveAs("output/bbcx/factor.pdf");

	out.close();
	f_IN->Close();
	f_mb->Close();
}
