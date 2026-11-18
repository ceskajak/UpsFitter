// pT > 4 GeV/c version. Publication-style multiplicity dependence for projection: subcorr (daughters subtracted) and corr (no subtraction), Run17 as stars, with Run11.
// Boxes of both = Syst.cxx + BBCx model envelope in quadrature; envelope = max |S_model/S - 1| over the five models of the same subtraction choice.
void official_plot_261112_pt4_both(TString txt = "/mnt/reservoir/Research/repos/UpsFitter/output/multdep/multdep_261109_dca1_pt4.txt", TString out = "/home/jakub/workspace/scratch/bbcx_fit/official_261112_pt4")
{
	struct P { string var, bin; double x, y, st, sy, S; };
	map<string, vector<P>> sch;
	string line, cur;
	ifstream f(txt.Data());
	while (getline(f, line))
	{
		if (line.rfind("Scheme ", 0) == 0) { cur = line.substr(7); continue; }
		istringstream is(line);
		P p; double nmb;
		if (is >> p.var >> p.bin >> p.x >> p.y >> p.st >> p.sy >> p.S >> nmb) sch[cur].push_back(p);
	}
	map<string, vector<string>> models = {{"subcorr", {"subcorrm", "mbmul1", "mbmul2", "mbadd1", "mbadd2"}},
	                                      {"corr", {"corrm", "mbmul1ns", "mbmul2ns", "mbadd1ns", "mbadd2ns"}}};
	// Run 11 reference: HEPData ins2877008 Figure 10a, third table = Upsilon(1S), pT > 4 GeV/c (larger last bin; HEPData-ins2877008-v1-Figure_10a.csv)
	Double_t r11x[] = {0.53835, 1.3283, 2.2807, 3.3673}, r11y[] = {0.277, 1.33, 2.65, 11.1}, r11dy[] = {0.08, 0.2, 0.65, 3.4}, r11ds[] = {0.057, 0.1, 0.39, 2.7};

	// Okabe-Ito blue / vermilion, colour-blind safe and distinct on a projector
	int cSub = TColor::GetColor("#0072B2"), cCorr = TColor::GetColor("#D55E00"), cR11 = TColor::GetColor("#3A3A3A");
	const double dxCorr = 0.08;

	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(6);
	gStyle->SetTextFont(42);
	gStyle->SetLineScalePS(2);
	for (string name : {"5bin"})
	{
		auto & v = sch[name];
		auto Tot = [&](const P & p)
		{
			double env = 0;
			auto & mv = models[p.var];
			for (auto & q : v) if (q.bin == p.bin && find(mv.begin(), mv.end(), q.var) != mv.end()) env = max(env, fabs(q.S / p.S - 1));
			return sqrt(p.sy * p.sy + pow(env * p.y, 2));
		};
		double ymax = 0;
		for (auto & p : v) if (p.var == "subcorr" || p.var == "corr") ymax = max(ymax, p.y + max(p.st, Tot(p)));
		for (int i = 0; i < 4; i++) ymax = max(ymax, r11y[i] + r11dy[i]);
		ymax = max(7.0, ceil(1.1 * ymax));

		TCanvas * c = new TCanvas(Form("c_%s", name.c_str()), "", 1580, 1160);
		c->SetLeftMargin(0.13); c->SetBottomMargin(0.15); c->SetTopMargin(0.03); c->SetRightMargin(0.03);
		c->SetTicks(1, 1);
		TH1F * fr = new TH1F(Form("fr_%s", name.c_str()), "", 1, 0, 5.2);
		fr->SetMinimum(0); fr->SetMaximum(ymax);
		fr->GetXaxis()->SetTitle("N_{ch}/#LTN_{ch}#GT"); fr->GetYaxis()->SetTitle("N_{#Upsilon}/#LTN_{#Upsilon}#GT");
		for (auto ax : {fr->GetXaxis(), fr->GetYaxis()}) { ax->SetTitleSize(0.065); ax->SetLabelSize(0.058); ax->SetTitleFont(42); ax->SetLabelFont(42); ax->SetLabelOffset(0.012); }
		fr->GetXaxis()->SetTitleOffset(1.0); fr->GetYaxis()->SetTitleOffset(0.9);
		fr->GetXaxis()->CenterTitle(); fr->GetYaxis()->CenterTitle();
		fr->GetXaxis()->SetNdivisions(506); fr->GetYaxis()->SetNdivisions(507);
		fr->SetLineWidth(3);
		fr->Draw("AXIS");
		TF1 * lin = new TF1(Form("lin_%s", name.c_str()), "x", 0, 5.2); lin->SetLineStyle(2); lin->SetLineColor(kGray + 2); lin->SetLineWidth(4); lin->Draw("same");

		for (int i = 0; i < 4; i++)
		{
			TBox * b = new TBox(r11x[i] - 0.07, r11y[i] - r11ds[i], r11x[i] + 0.07, r11y[i] + r11ds[i]); b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineColor(cR11); b->SetLineWidth(2); b->Draw("l");
		}
		TGraphErrors * g11 = new TGraphErrors(4, r11x, r11y, nullptr, r11dy);
		g11->SetMarkerStyle(21); g11->SetMarkerSize(2.6); g11->SetMarkerColor(cR11); g11->SetLineColor(cR11); g11->SetLineWidth(4); g11->Draw("P SAME");

		TGraphErrors * gs = new TGraphErrors(), * gc = new TGraphErrors();
		TBox * boxS = nullptr, * boxC = nullptr;
		for (auto & p : v)
		{
			bool isSub = p.var == "subcorr", isCorr = p.var == "corr";
			if (!isSub && !isCorr) continue;
			double x = isCorr ? p.x + dxCorr : p.x, tot = Tot(p);
			int col = isSub ? cSub : cCorr;
			TBox * b = new TBox(x - 0.09, max(0.0, p.y - tot), x + 0.09, min(ymax, p.y + tot)); b->SetFillColorAlpha(col, 0.25); b->SetLineColor(col); b->SetLineWidth(3); b->Draw("l");
			TGraphErrors * g = isSub ? gs : gc;
			g->SetPoint(g->GetN(), x, p.y); g->SetPointError(g->GetN() - 1, 0, p.st);
			(isSub ? boxS : boxC) = b;
		}
		gc->SetMarkerStyle(29); gc->SetMarkerSize(4.2); gc->SetMarkerColor(cCorr); gc->SetLineColor(cCorr); gc->SetLineWidth(4); gc->Draw("P SAME");
		gs->SetMarkerStyle(29); gs->SetMarkerSize(4.2); gs->SetMarkerColor(cSub); gs->SetLineColor(cSub); gs->SetLineWidth(4); gs->Draw("P SAME");

		TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextColor(kRed + 1); t.SetTextSize(0.062);
		t.DrawLatex(0.155, 0.885, "STAR Preliminary");
		t.SetTextFont(42); t.SetTextColor(kBlack); t.SetTextSize(0.039);
		t.DrawLatex(0.155, 0.82, "p+p #sqrt{s} = 510 GeV, BHT2");
		t.DrawLatex(0.155, 0.76, "#Upsilon(1S+2S+3S) #rightarrow e^{+}e^{-}, |y| < 1, p_{T} > 4 GeV/c");

		// legend entries show marker, stat bar and syst box together (fill taken from the graph)
		gs->SetFillColorAlpha(cSub, 0.25); gc->SetFillColorAlpha(cCorr, 0.25); g11->SetFillColorAlpha(kGray + 1, 0.35);
		TLegend * lg = new TLegend(0.15, 0.40, 0.55, 0.735); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038); lg->SetMargin(0.16);
		lg->SetHeader("Run 17 preliminary");
		((TLegendEntry *)lg->GetListOfPrimitives()->First())->SetTextFont(62);
		lg->AddEntry(gs, "subtraction", "pef");
		lg->AddEntry(gc, "no subtraction (x+0.08)", "pef");
		lg->AddEntry(g11, "#splitline{Run 11, PRD 112, 032004}{#Upsilon(1S), p_{T} > 4 GeV/c}", "pef");
		lg->AddEntry(lin, "y = x", "l");
		lg->Draw();
		c->RedrawAxis();
		c->SaveAs(Form("%s_%s.png", out.Data(), name.c_str()));
		c->SaveAs(Form("%s_%s.pdf", out.Data(), name.c_str()));
	}
}
