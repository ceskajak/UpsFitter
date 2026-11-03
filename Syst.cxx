#include "headers/results.hxx"
#include "headers/plots.hxx"

#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"

/////////////////////////////////////////////////////////////
// Systematic uncertainty of the self-normalised yield f = N_bin/N_int per multiplicity bin (Upsilon_systematics in the vault)
//   per source: delta = max over the variations of |f_var/f_nom - 1|
//   cut sources: kept only if |rel| > Barlow sigma of the variation that gives the maximum, otherwise set to zero
//   fit-model sources: always kept, except the diagnostic source sep (unconstrained distance scale, shown but not summed)
//   total per bin: quadrature sum of the kept sources
//   count: the nominal fits read as bin counting (US - LS - fitted residual background, divided by the fitted signal fraction in the window)
// Not included (no input in UpsFitter): smearing, efficiency vs multiplicity

struct SourceDelta
{
	Int_t n_usable = 0;		// variations that entered
	Double_t rel = 0;		// signed relative deviation of the variation with the largest |rel|
	Double_t sig = 0;		// Barlow sigma of that variation
	TString var = "";
	Bool_t kept = true;
	Double_t val = 0;		// retained |rel|
};

static SourceDelta MakeSourceDelta(const Source &s, map<TString, map<TString, Result>> &res, TString bin)
{
	SourceDelta d;
	for (auto v : s.vars)
	{
		Delta x = MakeDelta(res["nom"]["int"], res["nom"][bin], res[v]["int"], res[v][bin]);
		if (!x.ok) continue;
		d.n_usable++;
		if (d.var == "" || fabs(x.rel) > fabs(d.rel))
		{
			d.rel = x.rel;
			d.sig = x.sig;
			d.var = v;
		}
	}
	d.kept = !s.diag && (!s.cut || fabs(d.rel) > d.sig);
	d.val = d.kept ? fabs(d.rel) : 0;
	return d;
}

/////////////////////////////////////////////////////////////
// Run after Summary: root -l -b -q 'Syst.cxx("<YYMMDD>", "<kMULT>")', defaults are today and subcorr
void Syst(TString kDATE = "", TString kMULT = "subcorr")
{
	if (kDATE == "") kDATE = TodayYYMMDD();
	gStyle->SetOptStat(0);

	gSystem->mkdir("output/summary", true);
	TString s_out = "output/summary/syst_" + kDATE + "_" + kMULT;

	vector<Source> sources = MakeSources();
	Int_t ns = sources.size();
	Int_t nb = vec_mult_bins.size();

	map<TString, map<TString, Result>> res;
	for (auto v : MakeVars()) res[v] = LoadVariation(v, kMULT, kDATE);
	res["count"] = LoadVariation("count", kMULT, kDATE);

	// Nominal yield fraction and its statistical error (N_bin and N_int treated as uncorrelated, conservative)
	vector<Double_t> f_nom(nb, 0), f_stat(nb, 0), tot(nb, 0);
	vector<vector<SourceDelta>> d(ns, vector<SourceDelta>(nb));
	for (int ib = 0; ib < nb; ib++)
	{
		Result &ri = res["nom"]["int"];
		Result &rb = res["nom"][vec_mult_bins[ib]];
		if (Usable(ri) && Usable(rb))
		{
			f_nom[ib] = rb.S/ri.S;
			f_stat[ib] = f_nom[ib]*sqrt(pow(rb.dS/rb.S, 2) + pow(ri.dS/ri.S, 2));
		}
		Double_t sum2 = 0;
		for (int is = 0; is < ns; is++)
		{
			d[is][ib] = MakeSourceDelta(sources[is], res, vec_mult_bins[ib]);
			sum2 += pow(d[is][ib].val, 2);
		}
		tot[ib] = sqrt(sum2);
	}

	/*
	/////////////////////////////////////////////////////////////
	TEXT TABLE
	/////////////////////////////////////////////////////////////
	*/

	ofstream out;
	out.open((s_out + ".txt").Data());
	out << "Systematics of f = N_bin/N_int, " << kDATE << ", kMULT=" << kMULT << endl;
	out << "Relative deviation in %, per source the maximum over its variations." << endl;
	out << "(x) cut source below its Barlow sigma, or diagnostic source (sep), not in the total; n/a no usable variation; ! some variations unusable (FAILED or MISSING)." << endl;
	out << "Not included: smearing, efficiency vs multiplicity." << endl << endl;

	out << Form("%-12s", "source");
	for (auto b : vec_mult_bins) out << Form("%11s", b.Data());
	out << endl;
	for (int is = 0; is < ns; is++)
	{
		out << Form("%-12s", sources[is].name.Data());
		for (int ib = 0; ib < nb; ib++)
		{
			SourceDelta &x = d[is][ib];
			TString cell;
			if (x.n_usable == 0) cell = "n/a";
			else
			{
				cell = x.kept ? Form("%.2f", 100*x.val) : Form("(%.2f)", 100*fabs(x.rel));
				if (x.n_usable < (Int_t)sources[is].vars.size()) cell += "!";
			}
			out << Form("%11s", cell.Data());
		}
		out << endl;
	}
	out << Form("%-12s", "total syst");
	for (int ib = 0; ib < nb; ib++) out << Form("%11.2f", 100*tot[ib]);
	out << endl;
	out << Form("%-12s", "stat");
	for (int ib = 0; ib < nb; ib++) out << Form("%11.2f", f_nom[ib] > 0 ? 100*f_stat[ib]/f_nom[ib] : 0.);
	out << endl << endl;

	out << "Largest variation per source and bin (variation: signed deviation in %, Barlow sigma in %)" << endl;
	for (int is = 0; is < ns; is++)
		for (int ib = 0; ib < nb; ib++)
		{
			SourceDelta &x = d[is][ib];
			if (x.n_usable == 0) continue;
			out << Form("%-12s %-8s %-22s %8.2f %8.2f", sources[is].name.Data(), vec_mult_bins[ib].Data(), x.var.Data(), 100*x.rel, 100*x.sig) << endl;
		}
	out << endl;

	out << Form("%-8s %10s %10s %10s", "bin", "f_nom", "stat", "syst") << endl;
	for (int ib = 0; ib < nb; ib++)
		out << Form("%-8s %10.5f %10.5f %10.5f", vec_mult_bins[ib].Data(), f_nom[ib], f_stat[ib], f_nom[ib]*tot[ib]) << endl;
	out.close();

	cout << "Total systematic (%), stat (%):" << endl;
	for (int ib = 0; ib < nb; ib++)
		cout << Form("  %-8s syst %6.2f  stat %6.2f", vec_mult_bins[ib].Data(), 100*tot[ib], f_nom[ib] > 0 ? 100*f_stat[ib]/f_nom[ib] : 0.) << endl;

	/*
	/////////////////////////////////////////////////////////////
	GRAPHICS
	/////////////////////////////////////////////////////////////
	*/

	TCanvas * cnv = new TCanvas("cnv", "cnv", 1800, 1000);

	// 1. Deviation per source (filled: kept, open: cut source dropped by the Barlow criterion) and total against stat
	cnv->Divide(2, 1, 0.005, 0.005);
	cnv->cd(1);
	gPad->SetLogy();
	gPad->SetBottomMargin(0.22);
	gPad->SetLeftMargin(0.14);
	vector<TString> names;
	for (auto s : sources) names.push_back(s.name);
	TH1F * h1 = MakeFrame("h_src", names, 0.01, 200, "|f_{var}/f_{nom} - 1| (%)");
	h1->GetYaxis()->SetTitleSize(0.035);
	h1->GetYaxis()->SetTitleOffset(1.4);
	h1->Draw();
	TLegend * leg = new TLegend(0.14, 0.72, 0.4, 0.9);
	for (int ib = 0; ib < nb; ib++)
	{
		TGraphErrors * gk = new TGraphErrors();
		TGraphErrors * gd = new TGraphErrors();
		for (int is = 0; is < ns; is++)
		{
			SourceDelta &x = d[is][ib];
			if (x.n_usable == 0 || x.rel == 0) continue;
			TGraphErrors * g = x.kept ? gk : gd;
			g->SetPoint(g->GetN(), is + 0.1 + 0.8*(ib + 0.5)/nb, 100*fabs(x.rel));
		}
		SetPoint(gk, vec_colors[ib], 20, 1.0);
		SetPoint(gd, vec_colors[ib], 24, 1.0);
		if (gk->GetN()) gk->Draw("P SAME");
		if (gd->GetN()) gd->Draw("P SAME");
		leg->AddEntry(gk, vec_mult_bins[ib], "p");
	}
	leg->SetNColumns(2);
	leg->Draw();

	cnv->cd(2);
	gPad->SetBottomMargin(0.22);
	gPad->SetLeftMargin(0.12);
	TH1F * h2 = MakeFrame("h_tot", vec_mult_bins, 0, 100, "relative uncertainty of f (%)");
	Double_t ymax = 5;
	for (int ib = 0; ib < nb; ib++) ymax = max(ymax, 1.2*100*max(tot[ib], f_nom[ib] > 0 ? f_stat[ib]/f_nom[ib] : 0.));
	h2->SetMaximum(ymax);
	h2->Draw();
	TGraphErrors * g_stat = new TGraphErrors();
	TGraphErrors * g_syst = new TGraphErrors();
	for (int ib = 0; ib < nb; ib++)
	{
		g_stat->SetPoint(ib, ib + 0.35, f_nom[ib] > 0 ? 100*f_stat[ib]/f_nom[ib] : 0.);
		g_syst->SetPoint(ib, ib + 0.65, 100*tot[ib]);
	}
	SetPoint(g_stat, kBlack, 20, 1.2);
	SetPoint(g_syst, kRed + 1, 21, 1.2);
	g_stat->Draw("P SAME");
	g_syst->Draw("P SAME");
	TLegend * leg2 = new TLegend(0.15, 0.8, 0.4, 0.9);
	leg2->AddEntry(g_stat, "statistical", "p");
	leg2->AddEntry(g_syst, "systematic (total)", "p");
	leg2->Draw();

	cnv->cd();
	Title("Systematic sources (filled: kept, open: cut source below its Barlow sigma or diagnostic source) and total against statistical uncertainty of N_{bin}/N_{int}");
	cnv->Print(s_out + ".pdf(");
	cnv->SaveAs(s_out + "_1_sources.png");

	// 2. f = N_bin/N_int with statistical bars and systematic boxes, one pad per bin scheme
	cnv->Clear();
	cnv->Divide(2, 1, 0.005, 0.005);
	vector<vector<TString>> schemes = {vec_scheme5, vec_scheme4};
	vector<TString> scheme_names = {"5 bins {0,1,2,3,4,8}", "4 bins {0,1,2,3,8}"};
	for (int ip = 0; ip < 2; ip++)
	{
		cnv->cd(ip + 1);
		gPad->SetBottomMargin(0.15);
		gPad->SetLeftMargin(0.14);
		Double_t fmax = 0;
		for (auto b : schemes[ip])
		{
			Int_t ib = find(vec_mult_bins.begin(), vec_mult_bins.end(), b) - vec_mult_bins.begin();
			fmax = max(fmax, f_nom[ib]*(1 + tot[ib]) + f_stat[ib]);
		}
		TH1F * h = MakeFrame(Form("h_f%d", ip), schemes[ip], 0, 1.2*fmax, "N_{bin}/N_{int}");
		h->GetXaxis()->LabelsOption("h");
		h->GetXaxis()->SetLabelSize(0.05);
		h->Draw();
		TGraphErrors * g = new TGraphErrors();
		for (int i = 0; i < (int)schemes[ip].size(); i++)
		{
			Int_t ib = find(vec_mult_bins.begin(), vec_mult_bins.end(), schemes[ip][i]) - vec_mult_bins.begin();
			TBox * box = new TBox(i + 0.3, f_nom[ib]*(1 - tot[ib]), i + 0.7, f_nom[ib]*(1 + tot[ib]));
			box->SetFillColorAlpha(kRed + 1, 0.3);
			box->SetLineColor(kRed + 1);
			box->Draw();
			g->SetPoint(i, i + 0.5, f_nom[ib]);
			g->SetPointError(i, 0, f_stat[ib]);
		}
		SetPoint(g, kBlack, 20, 1.0);
		g->Draw("P SAME");
	}
	cnv->cd();
	Title("N_{bin}/N_{int} (nominal): statistical error bars, systematic boxes (red), " + kMULT);
	cnv->Print(s_out + ".pdf)");
	cnv->SaveAs(s_out + "_2_fraction.png");

	cout << "Written " << s_out << ".txt, " << s_out << ".pdf and 2 PNG pages" << endl;
}

/////////////////////////////////////////////////////////////
