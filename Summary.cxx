#include "headers/results.hxx"
#include "headers/plots.hxx"

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMarker.h"

/////////////////////////////////////////////////////////////

TString s_out;		// output path without extension
Int_t n_page = 0;

// Adds the canvas as a page of the PDF and as a PNG
static void SavePage(TCanvas * cnv, TString name)
{
	n_page++;
	cnv->Print(s_out + ".pdf" + (n_page == 1 ? "(" : ""));
	cnv->SaveAs(Form("%s_%d_%s.png", s_out.Data(), n_page, name.Data()));
}

static Int_t FlagColor(TString flag)
{
	if (flag == "OK") return kGreen + 1;
	if (flag == "WARN") return kOrange;
	if (flag == "FAILED") return kRed;
	return kGray;
}

/////////////////////////////////////////////////////////////
// Run after the jobs finished: ./Summary.cxx("<YYMMDD>", "<kMULT>"), defaults are today and subcorr
void Summary(TString kDATE = "", TString kMULT = "subcorr")
{
	if (kDATE == "") kDATE = TodayYYMMDD();
	gStyle->SetOptStat(0);

	gSystem->mkdir("output/summary", true);
	s_out = "output/summary/summary_" + kDATE + "_" + kMULT;

	vector<TString> vars = MakeVars();
	vector<Source> sources = MakeSources();
	Int_t nv = vars.size();
	Int_t nb = vec_sys_bins.size();

	map<TString, map<TString, Result>> res;
	for (auto v : vars) res[v] = LoadVariation(v, kMULT, kDATE);

	/*
	/////////////////////////////////////////////////////////////
	TEXT TABLE
	/////////////////////////////////////////////////////////////
	*/

	ofstream out;
	out.open((s_out + ".txt").Data());
	Int_t n_ok = 0, n_warn = 0, n_failed = 0, n_missing = 0;
	TString s_flagged = "";
	for (auto v : vars)
		for (auto b : vec_sys_bins)
		{
			Result &r = res[v][b];
			if (r.flag == "OK") n_ok++;
			else if (r.flag == "WARN") n_warn++;
			else if (r.flag == "FAILED") n_failed++;
			else n_missing++;
			if (r.flag == "FAILED" || r.flag == "MISSING")
				s_flagged += Form("  %-8s %-18s %-8s%s\n", b.Data(), v.Data(), r.flag.Data(), r.reason.Data());
		}

	TString s_head = Form("Summary %s, kMULT=%s: %d fits, OK %d, WARN %d, FAILED %d, MISSING %d\n", kDATE.Data(), kMULT.Data(), nv*nb, n_ok, n_warn, n_failed, n_missing);
	s_head += (s_flagged == "") ? "No failed or missing fits.\n" : "FAILED or MISSING (excluded from the systematics):\n" + s_flagged;
	cout << s_head;
	out << s_head << endl;

	out << Form("%-18s %-8s %-8s %12s %9s %7s %6s %6s %5s  %s", "variation", "bin", "flag", "S", "dS", "chi2", "att", "nlim", "covq", "reason") << endl;
	for (auto v : vars)
		for (auto b : vec_sys_bins)
		{
			Result &r = res[v][b];
			Int_t nlim = 0;
			for (auto p : r.pars) if (p.at_limit) nlim++;
			out << Form("%-18s %-8s %-8s %12.2f %9.2f %7.2f %6d %6d %5d %s", v.Data(), b.Data(), r.flag.Data(), r.S, r.dS, r.chi2, r.attempts, nlim, r.covq, r.reason.Data()) << endl;
		}
	out.close();

	/*
	/////////////////////////////////////////////////////////////
	GRAPHICS
	/////////////////////////////////////////////////////////////
	*/

	TCanvas * cnv = new TCanvas("cnv", "cnv", 1800, 1300);

	// 1. Status matrix
	{
		cnv->Clear();
		cnv->cd();
		gPad->SetLeftMargin(0.17);
		gPad->SetTopMargin(0.07);
		gPad->SetBottomMargin(0.06);
		TH2F * h = new TH2F("h_status", "", nb, 0, nb, nv, 0, nv);
		for (int ib = 0; ib < nb; ib++) h->GetXaxis()->SetBinLabel(ib + 1, vec_sys_bins[ib]);
		for (int iv = 0; iv < nv; iv++) h->GetYaxis()->SetBinLabel(nv - iv, vars[iv]);
		h->GetXaxis()->SetLabelSize(0.035);
		h->GetYaxis()->SetLabelSize(0.025);
		h->Draw("AXIS");

		TLatex t;
		t.SetTextAlign(22);
		t.SetTextSize(0.016);
		for (int iv = 0; iv < nv; iv++)
			for (int ib = 0; ib < nb; ib++)
			{
				Result &r = res[vars[iv]][vec_sys_bins[ib]];
				TBox * box = new TBox(ib + 0.03, nv - iv - 1 + 0.05, ib + 0.97, nv - iv - 0.05);
				box->SetFillColor(FlagColor(r.flag));
				box->Draw();
				if (r.found) t.DrawLatex(ib + 0.5, nv - iv - 0.5, Form("%.0f#pm%.0f", r.S, r.dS));
			}
		Title(Form("Fit status (green OK, orange WARN, red FAILED, grey MISSING) and S #pm #deltaS, %s, kMULT=%s", kDATE.Data(), kMULT.Data()));
		SavePage(cnv, "status");
	}

	// 2. Yield per variation, one pad per bin
	{
		cnv->Clear();
		cnv->Divide(4, 2, 0.005, 0.005);
		for (int ib = 0; ib < nb; ib++)
		{
			TString b = vec_sys_bins[ib];
			cnv->cd(ib + 1);
			gPad->SetBottomMargin(0.32);
			gPad->SetLeftMargin(0.14);

			Double_t ymin = 1e30, ymax = -1e30;
			for (auto v : vars)
			{
				Result &r = res[v][b];
				if (!Usable(r)) continue;
				ymin = min(ymin, r.S - r.dS);
				ymax = max(ymax, r.S + r.dS);
			}
			if (ymin > ymax) { ymin = 0; ymax = 1; }
			Double_t pad = 0.15*(ymax - ymin);

			TH1F * h = MakeFrame("h_yield_" + b, vars, ymin - pad, ymax + pad, "S (" + b + ")");
			h->Draw();
			Result &nom = res["nom"][b];
			if (Usable(nom))
			{
				TLine * l = new TLine(0, nom.S, nv, nom.S);
				l->SetLineStyle(2);
				l->Draw();
			}
			for (TString flag : {"OK", "WARN", "FAILED"})
			{
				TGraphErrors * g = new TGraphErrors();
				for (int iv = 0; iv < nv; iv++)
				{
					Result &r = res[vars[iv]][b];
					if (r.flag != flag || !std::isfinite(r.S) || !std::isfinite(r.dS)) continue;
					g->SetPoint(g->GetN(), iv + 0.5, r.S);
					g->SetPointError(g->GetN() - 1, 0, r.dS);
				}
				if (g->GetN() == 0) continue;
				SetPoint(g, flag == "OK" ? kBlack : (flag == "WARN" ? kOrange + 1 : kRed), 20, 0.8);
				g->Draw("P SAME");
			}
		}
		cnv->cd();
		Title("Yield S per variation (dashed: nominal; orange WARN, red FAILED)");
		SavePage(cnv, "yields");
	}

	// 3. Deviation of the yield fraction N_bin/N_int from nominal, one pad per multiplicity bin
	{
		cnv->Clear();
		cnv->Divide(3, 2, 0.005, 0.005);
		vector<TString> vars_var(vars.begin() + 1, vars.end());
		for (int ib = 0; ib < (int)vec_mult_bins.size(); ib++)
		{
			TString b = vec_mult_bins[ib];
			cnv->cd(ib + 1);
			gPad->SetBottomMargin(0.32);
			gPad->SetLeftMargin(0.17);

			vector<Delta> deltas;
			Double_t ymax = 0.5;
			for (auto v : vars_var)
			{
				Delta d = MakeDelta(res["nom"]["int"], res["nom"][b], res[v]["int"], res[v][b]);
				deltas.push_back(d);
				if (d.ok) ymax = max(ymax, 100*(fabs(d.rel) + d.sig));
			}

			TH1F * h = MakeFrame("h_delta_" + b, vars_var, -1.2*ymax, 1.2*ymax, "#Delta (N_{" + b + "}/N_{int}) / nominal (%)");
			h->Draw();
			TLine * l0 = new TLine(0, 0, vars_var.size(), 0);
			l0->SetLineStyle(2);
			l0->Draw();

			// Fit-model sources are always kept, cut sources only when |Delta| exceeds the Barlow sigma
			for (Int_t color : {(Int_t)kBlue, (Int_t)kRed, (Int_t)kGray + 1})
			{
				TGraphErrors * g = new TGraphErrors();
				for (int iv = 0; iv < (int)vars_var.size(); iv++)
				{
					if (!deltas[iv].ok) continue;
					Bool_t is_cut = vars_var[iv].BeginsWith("cut_");
					Bool_t kept = !is_cut || fabs(deltas[iv].rel) > deltas[iv].sig;
					Int_t c = !is_cut ? kBlue : (kept ? kRed : kGray + 1);
					if (c != color) continue;
					g->SetPoint(g->GetN(), iv + 0.5, 100*deltas[iv].rel);
					g->SetPointError(g->GetN() - 1, 0, 100*deltas[iv].sig);
				}
				if (g->GetN() == 0) continue;
				SetPoint(g, color, 20, 0.8);
				g->Draw("P SAME");
			}
		}
		cnv->cd();
		Title("Relative deviation of N_{bin}/N_{int} from nominal; bar: Barlow #sigma (blue fit model, red cut kept, grey cut below #sigma)");
		SavePage(cnv, "deltas");
	}

	// 4. Position of every floating parameter inside its bounds (0 lower bound, 1 upper bound), one pad per bin
	{
		cnv->Clear();
		cnv->Divide(4, 2, 0.005, 0.005);
		for (int ib = 0; ib < nb; ib++)
		{
			TString b = vec_sys_bins[ib];
			cnv->cd(ib + 1);
			gPad->SetBottomMargin(0.32);
			gPad->SetLeftMargin(0.12);

			vector<TString> names;
			for (auto v : vars)
				for (auto p : res[v][b].pars)
					if (find(names.begin(), names.end(), p.name) == names.end()) names.push_back(p.name);

			TH1F * h = MakeFrame("h_par_" + b, names, -0.1, 1.1, "(value - lo)/(hi - lo) (" + b + ")");
			h->Draw();
			for (Double_t y : {0., 1.})
			{
				TLine * l = new TLine(0, y, names.size(), y);
				l->SetLineStyle(2);
				l->Draw();
			}
			for (int iv = nv - 1; iv >= 0; iv--)
			{
				Result &r = res[vars[iv]][b];
				if (!r.found) continue;
				TGraphErrors * g = new TGraphErrors();
				for (auto p : r.pars)
				{
					Int_t ip = find(names.begin(), names.end(), p.name) - names.begin();
					g->SetPoint(g->GetN(), ip + 0.5 + (iv == 0 ? 0 : 0.3*((iv % 7)/6. - 0.5)), (p.val - p.lo)/(p.hi - p.lo));
				}
				SetPoint(g, iv == 0 ? kBlack : (r.flag == "FAILED" ? kRed : kAzure + 1), iv == 0 ? 24 : 20, iv == 0 ? 1.5 : 0.5);
				g->Draw("P SAME");
			}
		}
		cnv->cd();
		Title("Floating parameters inside their bounds: 0 lower bound, 1 upper bound (open black: nominal, blue: variations, red: FAILED fits)");
		SavePage(cnv, "parameters");
	}

	// 5. Lineshape and fraction parameters vs bin (nominal) against their embedding bounds
	{
		cnv->Clear();
		cnv->Divide(4, 3, 0.005, 0.005);
		vector<TString> names = {"var_sig1_mu", "var_sig1_sigma", "var_sig2_sigma", "var_sig3_sigma", "var_sig1_alpha", "var_sig2_alpha", "var_sig3_alpha", "var_sig1_n", "var_sig2_n", "var_sig3_n", "var_sig1_frac", "var_sig2_frac"};
		for (int in = 0; in < (int)names.size(); in++)
		{
			cnv->cd(in + 1);
			gPad->SetBottomMargin(0.15);
			gPad->SetLeftMargin(0.18);

			Double_t lo = 0, hi = 0;
			TGraphErrors * g = new TGraphErrors();
			for (int ib = 0; ib < nb; ib++)
			{
				Result &r = res["nom"][vec_sys_bins[ib]];
				for (auto p : r.pars)
				{
					if (p.name != names[in]) continue;
					lo = p.lo;
					hi = p.hi;
					g->SetPoint(g->GetN(), ib + 0.5, p.val);
					g->SetPointError(g->GetN() - 1, 0, p.err);
				}
			}
			if (g->GetN() == 0) continue;
			Double_t pad = 0.3*(hi - lo);
			TH1F * h = MakeFrame("h_shape_" + names[in], vec_sys_bins, lo - pad, hi + pad, names[in]);
			h->GetXaxis()->SetLabelSize(0.06);
			h->Draw();
			TBox * band = new TBox(0, lo, nb, hi);
			band->SetFillColorAlpha(kGreen + 1, 0.25);
			band->Draw();
			SetPoint(g, kBlack, 20, 0.9);
			g->Draw("P SAME");
		}
		cnv->cd();
		Title("Nominal lineshape and fraction parameters per bin (green: BHT2_int embedding bounds)");
		SavePage(cnv, "lineshape");
	}

	// 6. chi2/ndf and S/sqrt(US) per variation, one graph per bin
	{
		cnv->Clear();
		cnv->Divide(1, 2, 0.005, 0.005);
		TLegend * leg = new TLegend(0.89, 0.55, 0.995, 0.9);
		for (int ip = 0; ip < 2; ip++)
		{
			cnv->cd(ip + 1);
			gPad->SetBottomMargin(0.3);
			gPad->SetRightMargin(0.12);
			Double_t ymax = 0;
			for (auto v : vars)
				for (auto b : vec_sys_bins)
				{
					Result &r = res[v][b];
					if (r.found && std::isfinite(r.chi2)) ymax = max(ymax, ip == 0 ? r.chi2 : r.SoverRoot);
				}
			TH1F * h = MakeFrame(Form("h_q%d", ip), vars, 0, 1.2*ymax, ip == 0 ? "#chi^{2}/ndf (US)" : "S/#sqrt{US}");
			h->Draw();
			if (ip == 0)
			{
				TLine * l = new TLine(0, 3.0, nv, 3.0);
				l->SetLineStyle(2);
				l->Draw();
			}
			for (int ib = 0; ib < nb; ib++)
			{
				TGraphErrors * g = new TGraphErrors();
				for (int iv = 0; iv < nv; iv++)
				{
					Result &r = res[vars[iv]][vec_sys_bins[ib]];
					if (!r.found) continue;
					g->SetPoint(g->GetN(), iv + 0.5, ip == 0 ? r.chi2 : r.SoverRoot);
				}
				SetPoint(g, vec_colors[ib], 20 + ib, 0.8);
				g->Draw("P SAME");
				if (ip == 0) leg->AddEntry(g, vec_sys_bins[ib], "p");
			}
			if (ip == 0) leg->Draw();
		}
		cnv->cd();
		Title("Fit quality per variation (dashed: chi2/ndf warning limit)");
		SavePage(cnv, "quality");
	}

	// 7. int validation: sum of the bin yields of each scheme against the integrated yield
	{
		cnv->Clear();
		cnv->Divide(1, 2, 0.005, 0.005);
		vector<vector<TString>> schemes = {vec_scheme5, vec_scheme4};
		vector<TString> scheme_names = {"5 bins", "4 bins"};
		TLegend * leg = new TLegend(0.89, 0.7, 0.995, 0.9);
		for (int ip = 0; ip < 2; ip++)
		{
			cnv->cd(ip + 1);
			gPad->SetBottomMargin(0.3);
			gPad->SetRightMargin(0.12);
			vector<TGraphErrors *> graphs;
			Double_t ymax = ip == 0 ? 5 : 3;
			for (int is = 0; is < 2; is++)
			{
				TGraphErrors * g = new TGraphErrors();
				for (int iv = 0; iv < nv; iv++)
				{
					Result &ri = res[vars[iv]]["int"];
					Double_t sum = 0, var = 0;
					Bool_t ok = Usable(ri);
					for (auto b : schemes[is])
					{
						Result &r = res[vars[iv]][b];
						if (!Usable(r)) ok = false;
						sum += r.S;
						var += r.dS*r.dS;
					}
					if (!ok) continue;
					Double_t diff = sum - ri.S;
					Double_t value = (ip == 0) ? 100*diff/ri.S : diff/sqrt(var + ri.dS*ri.dS);
					g->SetPoint(g->GetN(), iv + 0.5, value);
					ymax = max(ymax, 1.2*fabs(value));
				}
				SetPoint(g, vec_colors[is + 1], 20 + 3*is, 1.0);
				graphs.push_back(g);
				if (ip == 0) leg->AddEntry(g, scheme_names[is], "p");
			}
			TH1F * h = MakeFrame(Form("h_val%d", ip), vars, -ymax, ymax, ip == 0 ? "(#SigmaN_{bins} - N_{int})/N_{int} (%)" : "pull (#SigmaN_{bins} - N_{int})/#sigma");
			h->Draw();
			TLine * l0 = new TLine(0, 0, nv, 0);
			l0->SetLineStyle(2);
			l0->Draw();
			for (auto g : graphs) g->Draw("P SAME");
			if (ip == 0) leg->Draw();
		}
		cnv->cd();
		Title("Validation: sum of the bin yields against the integrated fit (pull treats the bins and int as uncorrelated)");
		SavePage(cnv, "int_validation");
	}

	// 8. TofMult variants on the multiplicity bins (nominal), relative to the chosen kMULT
	{
		cnv->Clear();
		cnv->Divide(3, 2, 0.005, 0.005);
		map<TString, map<TString, Result>> res_mult;
		for (auto m : vec_mult_variants) res_mult[m] = LoadVariation("nom", m, kDATE);
		for (int ib = 0; ib < (int)vec_mult_bins.size(); ib++)
		{
			TString b = vec_mult_bins[ib];
			cnv->cd(ib + 1);
			gPad->SetBottomMargin(0.15);
			gPad->SetLeftMargin(0.16);
			Result &ref = res_mult[kMULT][b];

			Double_t ymax = 1.5;
			for (auto m : vec_mult_variants)
			{
				Result &r = res_mult[m][b];
				if (Usable(r) && Usable(ref)) ymax = max(ymax, 1.15*(r.S + r.dS)/ref.S);
			}
			TH1F * h = MakeFrame("h_mult_" + b, vec_mult_variants, 0, ymax, "S / S(" + kMULT + ") (" + b + ")");
			h->GetXaxis()->LabelsOption("h");
			h->GetXaxis()->SetLabelSize(0.06);
			h->Draw();
			TLine * l1 = new TLine(0, 1, vec_mult_variants.size(), 1);
			l1->SetLineStyle(2);
			l1->Draw();
			TGraphErrors * g = new TGraphErrors();
			for (int im = 0; im < (int)vec_mult_variants.size(); im++)
			{
				Result &r = res_mult[vec_mult_variants[im]][b];
				if (!Usable(r) || !Usable(ref)) continue;
				g->SetPoint(g->GetN(), im + 0.5, r.S/ref.S);
				g->SetPointError(g->GetN() - 1, 0, r.dS/ref.S);
			}
			SetPoint(g, kBlack, 20, 1.0);
			g->Draw("P SAME");
		}
		cnv->cd();
		Title("TofMult variants (none, sub = subtraction, corr = luminosity correction), yield relative to " + kMULT);
		SavePage(cnv, "tofmult");
	}

	cnv->Print(s_out + ".pdf]");
	cout << "Written " << s_out << ".txt, " << s_out << ".pdf and " << n_page << " PNG pages" << endl;
}

/////////////////////////////////////////////////////////////
