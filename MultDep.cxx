#include "headers/results.hxx"
#include "headers/plots.hxx"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TBox.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

/////////////////////////////////////////////////////////////
// Multiplicity dependence N_Ups^bin/<N_Ups> vs N_ch^bin/<N_ch> for every TofMult variant and both bin schemes
/////////////////////////////////////////////////////////////
// Same quantity as SCDF/scripts/mult_dep/plotter/final_{4,5}bin.cxx, now with systematic boxes from Syst.cxx
// Run after Syst: root -l -b -q 'MultDep.cxx("<YYMMDD[_dcaX]>")', the default is today
//   y = (S_bin/S_int) * (nMB_int/nMB_bin), x = <TofMult>_bin/3 from the minimum-bias TofMult distribution
//   bins are (3*lo, 3*hi] in TofMult, as in the fit selection; <TofMult> = 3 as in the previous macros
//   statistical error: fit errors of S_bin and S_int and the Poisson error of the MB counts, all uncorrelated (conservative, as in Syst.cxx)
//   systematic error: relative total of Syst.cxx for the same TofMult variant (only subcorr has the full systematic study), efficiency vs multiplicity and smearing not included

TString kBBC_FILE = "/mnt/reservoir/Research/data/BBC.histos_250221.root";
Double_t kMULT_MEAN = 3;
Int_t kMAX_MULT = 20;	// last TofMult bin of the MB histograms

struct MultPoint
{
	Bool_t ok = false;
	Double_t x = 0, y = 0, stat = 0, syst = 0;	// syst is absolute, 0 when there is no systematic study for the variant
	Double_t nUps = 0, nMB = 0;
};

// Relative total systematic per bin from syst_<date>_<mult>.txt (last table: bin, f_nom, stat, syst), empty when the file does not exist
static map<TString, Double_t> ReadSystRel(TString date, TString mult)
{
	map<TString, Double_t> rel;
	ifstream in(Form("output/summary/syst_%s_%s.txt", date.Data(), mult.Data()));
	string line;
	Bool_t table = false;
	while (getline(in, line))
	{
		if (line.find("f_nom") != string::npos) { table = true; continue; }
		if (!table) continue;
		char name[64];
		Double_t f, st, sy;
		if (sscanf(line.c_str(), "%63s %lf %lf %lf", name, &f, &st, &sy) == 4 && f > 0) rel[name] = sy/f;
	}
	return rel;
}

static MultPoint MakePoint(TString bin, TString mult, TString date, TH1D * h_MB, map<TString, Double_t> &syst_rel)
{
	MultPoint p;
	Result ri = ReadResult(ResultFile("int", "nom", "exp", mult, date));
	Result rb = ReadResult(ResultFile(bin, "nom", "exp", mult, date));
	if (!Usable(ri) || !Usable(rb)) return p;

	Int_t lo = TString(bin(0, bin.Index("mult"))).Atoi();
	Int_t hi = TString(bin(bin.Index("mult") + 4, bin.Length())).Atoi();

	// TofMult in (3*lo, 3*hi]: integer-valued variants fill bin m+1 for TofMult m, the corrected (continuous) ones are cut in [3*lo, 3*hi)
	Bool_t corrected = mult.Contains("corr");
	Int_t j_lo = kMULT_MEAN*lo + (corrected ? 1 : 2);
	Int_t j_hi = std::min((Int_t)(kMULT_MEAN*hi + (corrected ? 0 : 1)), kMAX_MULT);

	Double_t n = 0, sum = 0;
	for (int j = j_lo; j <= j_hi; j++)
	{
		n += h_MB->GetBinContent(j);
		sum += h_MB->GetBinContent(j)*(h_MB->GetBinLowEdge(j) + (corrected ? 0.5*h_MB->GetBinWidth(j) : 0));
	}
	Double_t n_int = h_MB->Integral();
	if (n <= 0) return p;

	p.x = sum/n/kMULT_MEAN;
	p.y = (rb.S/ri.S)*(n_int/n);
	p.stat = p.y*sqrt(pow(rb.dS/rb.S, 2) + pow(ri.dS/ri.S, 2) + 1.0/n + 1.0/n_int);
	p.syst = syst_rel.count(bin) ? p.y*syst_rel[bin] : 0;
	p.nUps = rb.S;
	p.nMB = n;
	p.ok = true;
	return p;
}

void MultDep(TString kDATE = "")
{
	if (kDATE == "") kDATE = TodayYYMMDD();
	gStyle->SetOptStat(0);

	gSystem->mkdir("output/multdep", true);
	TString s_out = "output/multdep/multdep_" + kDATE;

	TFile * f_BBC = TFile::Open(kBBC_FILE, "READ");
	if (!f_BBC || f_BBC->IsZombie())
	{
		cout << "Cannot open " << kBBC_FILE << endl;
		return;
	}
	TH1D * h_raw = (TH1D *)f_BBC->Get("h_TofMult");
	TH1D * h_corr = (TH1D *)f_BBC->Get("h_TofMult_corr");

	// TofMult variants: label, MB histogram follows the BBCx correction (the subtraction of the pair tracks has no MB counterpart)
	vector<TString> variants = {"subcorr", "sub", "corr", "none"};
	vector<TString> overlay_variants = {"subcorr", "corr"};	// variants drawn together on the overlay page
	auto IsOverlay = [&](TString v) { return std::find(overlay_variants.begin(), overlay_variants.end(), v) != overlay_variants.end(); };
	map<TString, TString> label = {
		{"subcorr", "pair tracks subtracted + BBCx corrected (nominal)"},
		{"sub", "pair tracks subtracted"},
		{"corr", "BBCx corrected"},
		{"none", "raw TofMult"}
	};
	Int_t colors[] = {kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1};

	vector<vector<TString>> schemes = {vec_scheme4, vec_scheme5};
	vector<TString> scheme_names = {"4bin", "5bin"};

	// Run 11 reference (Phys. Rev. D 112, 032004), 4 bins, Upsilon(nS) of HEPData ins2877008 Figure 10a (no x uncertainty published)
	TGraphAsymmErrors * gr_Run11 = new TGraphAsymmErrors(4);
	Double_t run11_x[] = {0.53835, 1.3283, 2.2807, 3.3673};
	Double_t run11_y[] = {0.41, 1.39, 2.67, 4.6};
	Double_t run11_dy[] = {0.064, 0.14, 0.43, 1.8};	// statistical
	Double_t run11_dsyst[] = {0.041, 0.085, 0.31, 1.3};	// systematic
	for (int i = 0; i < 4; i++)
	{
		gr_Run11->SetPoint(i, run11_x[i], run11_y[i]);
		gr_Run11->SetPointError(i, 0, 0, run11_dy[i], run11_dy[i]);
	}
	gr_Run11->SetMarkerStyle(kFullSquare);
	gr_Run11->SetMarkerSize(1.6);
	gr_Run11->SetMarkerColor(kBlack);
	gr_Run11->SetLineColor(kBlack);
	gr_Run11->SetLineWidth(2);

	ofstream out(s_out + ".txt");
	out << "Multiplicity dependence N_Ups^bin/<N_Ups> vs N_ch^bin/<N_ch>, " << kDATE << endl;
	out << "x = <TofMult>_bin/" << kMULT_MEAN << " from " << kBBC_FILE << ", y = (S_bin/S_int)*(nMB_int/nMB_bin)" << endl;
	out << "syst = total of Syst.cxx for the same variant (0: no systematic study), efficiency vs multiplicity and smearing not included" << endl << endl;

	TCanvas * cnv = new TCanvas("cnv", "cnv", 1600, 1200);
	cnv->SetMargin(0.1, 0.02, 0.1, 0.02);
	Int_t n_pages = 0;
	TString pdf = s_out + ".pdf";

	for (int is = 0; is < 2; is++)
	{
		vector<TString> &bins = schemes[is];
		map<TString, vector<MultPoint>> pts;
		for (auto v : variants)
		{
			map<TString, Double_t> syst_rel = ReadSystRel(kDATE, v);
			TH1D * h_MB = v.Contains("corr") ? h_corr : h_raw;
			for (auto b : bins) pts[v].push_back(MakePoint(b, v, kDATE, h_MB, syst_rel));
		}

		out << "Scheme " << scheme_names[is] << endl;
		out << Form("%-8s %-8s %8s %9s %9s %9s %9s %9s", "variant", "bin", "x", "y", "stat", "syst", "S", "nMB") << endl;
		for (auto v : variants)
			for (int ib = 0; ib < (int)bins.size(); ib++)
			{
				MultPoint &p = pts[v][ib];
				if (p.ok) out << Form("%-8s %-8s %8.4f %9.4f %9.4f %9.4f %9.1f %9.0f", v.Data(), bins[ib].Data(), p.x, p.y, p.stat, p.syst, p.nUps, p.nMB) << endl;
				else out << Form("%-8s %-8s no usable fit", v.Data(), bins[ib].Data()) << endl;
			}
		out << endl;

		// One page per variant (stat bars, syst boxes, Run 11) and one overlay page of the overlay variants
		for (int iv = 0; iv <= (int)variants.size(); iv++)
		{
			Bool_t overlay = (iv == (int)variants.size());
			cnv->Clear();

			Double_t ymax = 0;
			for (int jv = 0; jv < (int)variants.size(); jv++)
				if (overlay ? IsOverlay(variants[jv]) : jv == iv)
					for (auto &p : pts[variants[jv]])
						if (p.ok) ymax = max(ymax, p.y + sqrt(p.stat*p.stat + p.syst*p.syst));
			for (int i = 0; i < 4; i++) ymax = max(ymax, run11_y[i] + sqrt(run11_dy[i]*run11_dy[i] + run11_dsyst[i]*run11_dsyst[i]));

			TH1F * frame = new TH1F(Form("fr_%d_%d", is, iv), "", 1, 0, 5.5);
			frame->SetMinimum(0);
			frame->SetMaximum(1.5*ymax);
			frame->GetXaxis()->SetTitle("N_{ch}/<N_{ch}>");
			frame->GetYaxis()->SetTitle("N_{#Upsilon}^{bin}/<N_{#Upsilon}>");
			frame->Draw();

			TF1 * f_linear = new TF1(Form("f_lin_%d_%d", is, iv), "x", 0, 5.5);
			f_linear->SetLineColor(kGray + 2);
			f_linear->SetLineStyle(kDashed);
			f_linear->Draw("same");

			TLegend * leg = new TLegend(0.13, 0.62, 0.55, 0.96);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.026);

			for (int i = 0; i < 4; i++)
			{
				TBox * box11 = new TBox(run11_x[i] - 0.07, run11_y[i] - run11_dsyst[i], run11_x[i] + 0.07, run11_y[i] + run11_dsyst[i]);
				box11->SetFillColorAlpha(kGray + 1, 0.35);
				box11->SetLineColor(kGray + 1);
				box11->Draw();
			}
			gr_Run11->Draw("P SAME");
			Int_t n_drawn = 0;
			leg->AddEntry(gr_Run11, "Run11 - Phys. Rev. D 112, 032004 (HEPData, bars stat, boxes syst)", "pe");

			for (int jv = 0; jv < (int)variants.size(); jv++)
			{
				if (overlay ? !IsOverlay(variants[jv]) : jv != iv) continue;
				TString v = variants[jv];
				Double_t shift = overlay ? 0.04*(n_drawn++) : 0;
				TGraphErrors * g = new TGraphErrors();
				Bool_t has_syst = false;
				for (int ib = 0; ib < (int)bins.size(); ib++)
				{
					MultPoint &p = pts[v][ib];
					if (!p.ok) continue;
					Int_t col = overlay ? colors[jv] : kRed + 1;
					if (p.syst > 0)
					{
						has_syst = true;
						TBox * box = new TBox(p.x + shift - 0.07, p.y - p.syst, p.x + shift + 0.07, p.y + p.syst);
						box->SetFillColorAlpha(col, 0.25);
						box->SetLineColor(col);
						box->Draw();
					}
					g->SetPoint(g->GetN(), p.x + shift, p.y);
					g->SetPointError(g->GetN() - 1, 0, p.stat);
					if (!overlay)
					{
						TLatex * t = new TLatex(p.x + 0.1, p.y + 0.05*ymax, Form("%.2f #pm %.2f #pm %.2f", p.y, p.stat, p.syst));
						t->SetTextSize(0.025);
						t->SetTextColor(kBlue + 2);
						t->Draw();
					}
				}
				SetPoint(g, overlay ? colors[jv] : kBlue + 1, kFullCircle, 1.6);
				g->SetLineWidth(2);
				g->Draw("P SAME");
				leg->AddEntry(g, Form("Run17 - %s%s", label[v].Data(), overlay ? Form(" (x + %.2f)", shift) : ""), "pe");
				if (has_syst) leg->AddEntry((TObject *)0, "  boxes: systematic (Syst.cxx), bars: statistical", "");
			}
			leg->AddEntry(f_linear, "x = y", "l");
			leg->Draw();

			TLatex * title = new TLatex();
			title->SetNDC();
			title->SetTextSize(0.025);
			title->DrawLatex(0.13, 0.975, Form("%s, %s", scheme_names[is].Data(), (overlay ? TString("TofMult variants compared") : label[variants[iv]]).Data()));

			TString name = s_out + "_" + scheme_names[is] + "_" + (overlay ? TString("subcorr_corr") : variants[iv]);
			cnv->SaveAs(name + ".png");
			cnv->Print(pdf + (n_pages == 0 ? "(" : ""));
			n_pages++;
		}
	}
	// Close the pdf
	cnv->Clear();
	cnv->Print(pdf + ")");

	out.close();
	f_BBC->Close();
	cout << "Written " << s_out << ".txt, " << pdf << " and " << n_pages << " PNG pages" << endl;
}
