// pT-integrated vs pair pT > 4 GeV/c, 5-bin: left subtraction, right no subtraction. Boxes = Syst.cxx + BBCx model envelope in quadrature (envelope per set from its own five models).
struct PP { string var, bin; double x, y, st, sy, S; };
static map<string, vector<PP>> ReadMD(TString txt)
{
	map<string, vector<PP>> sch;
	string line, cur;
	ifstream f(txt.Data());
	while (getline(f, line))
	{
		if (line.rfind("Scheme ", 0) == 0) { cur = line.substr(7); continue; }
		istringstream is(line);
		PP p; double nmb;
		if (is >> p.var >> p.bin >> p.x >> p.y >> p.st >> p.sy >> p.S >> nmb) sch[cur].push_back(p);
	}
	return sch;
}

void official_plot_261112_ptcomp(TString out = "/home/jakub/workspace/scratch/bbcx_fit/official_261112_ptcomp_5bin")
{
	TString dir = "/mnt/reservoir/Research/repos/UpsFitter/output/multdep/";
	auto all = ReadMD(dir + "multdep_261109_dca1.txt")["5bin"], hi = ReadMD(dir + "multdep_261109_dca1_pt4.txt")["5bin"];
	map<string, vector<string>> models = {{"subcorr", {"subcorrm", "mbmul1", "mbmul2", "mbadd1", "mbadd2"}},
	                                      {"corr", {"corrm", "mbmul1ns", "mbmul2ns", "mbadd1ns", "mbadd2ns"}}};
	Double_t r11x[] = {0.53835, 1.3283, 2.2807, 3.3673}, r11y[] = {0.41, 1.39, 2.67, 4.6}, r11dy[] = {0.064, 0.14, 0.43, 1.8}, r11ds[] = {0.041, 0.085, 0.31, 1.3};
	Double_t h11y[] = {0.277, 1.33, 2.65, 11.1}, h11dy[] = {0.08, 0.2, 0.65, 3.4}, h11ds[] = {0.057, 0.1, 0.39, 2.7}; // HEPData ins2877008 Fig. 10a, Upsilon(1S) pT > 4
	int cSub = TColor::GetColor("#0072B2"), cCorr = TColor::GetColor("#D55E00"), cR11 = TColor::GetColor("#3A3A3A");
	const double dxHi = 0.1;

	gStyle->SetOptStat(0); gStyle->SetEndErrorSize(6); gStyle->SetTextFont(42); gStyle->SetLineScalePS(2);
	TCanvas * c = new TCanvas("c", "", 2400, 1160);
	c->Divide(2, 1, 0, 0);
	vector<string> vars = {"subcorr", "corr"};
	vector<TGraphErrors *> keep;
	for (int k = 0; k < 2; k++)
	{
		string var = vars[k];
		int col = k == 0 ? cSub : cCorr;
		TPad * pd = (TPad *)c->cd(k + 1);
		pd->SetLeftMargin(k == 0 ? 0.17 : 0.03); pd->SetRightMargin(k == 0 ? 0.01 : 0.03); pd->SetBottomMargin(0.15); pd->SetTopMargin(0.03); pd->SetTicks(1, 1);
		TH1F * fr = new TH1F(Form("fr%d", k), "", 1, 0, 5.4);
		fr->SetMinimum(0); fr->SetMaximum(16);
		fr->GetXaxis()->SetTitle("N_{ch}/#LTN_{ch}#GT"); fr->GetYaxis()->SetTitle("N_{#Upsilon}/#LTN_{#Upsilon}#GT");
		for (auto ax : {fr->GetXaxis(), fr->GetYaxis()}) { ax->SetTitleSize(0.065); ax->SetLabelSize(0.058); ax->SetTitleFont(42); ax->SetLabelFont(42); ax->SetLabelOffset(0.012); ax->CenterTitle(); }
		fr->GetXaxis()->SetTitleOffset(1.0); fr->GetYaxis()->SetTitleOffset(1.2);
		fr->GetXaxis()->SetNdivisions(506); fr->GetYaxis()->SetNdivisions(507);
		if (k == 1) fr->GetYaxis()->SetLabelSize(0);
		fr->Draw("AXIS");
		TF1 * lin = new TF1(Form("lin%d", k), "x", 0, 5.4); lin->SetLineStyle(2); lin->SetLineColor(kGray + 2); lin->SetLineWidth(4); lin->Draw("same");
		for (int i = 0; i < 4; i++) { TBox * b = new TBox(r11x[i] - 0.07, r11y[i] - r11ds[i], r11x[i] + 0.07, r11y[i] + r11ds[i]); b->SetFillColorAlpha(kGray + 1, 0.3); b->SetLineColor(cR11); b->SetLineWidth(2); b->Draw("l"); }
		for (int i = 0; i < 4; i++) { TBox * b = new TBox(r11x[i] - 0.07, h11y[i] - h11ds[i], r11x[i] + 0.07, h11y[i] + h11ds[i]); b->SetFillStyle(0); b->SetLineColor(cR11); b->SetLineStyle(7); b->SetLineWidth(2); b->Draw("l"); }
		TGraphErrors * gh11 = new TGraphErrors(4, r11x, h11y, nullptr, h11dy);
		gh11->SetMarkerStyle(25); gh11->SetMarkerSize(2.4); gh11->SetMarkerColor(cR11); gh11->SetLineColor(cR11); gh11->SetLineWidth(4); gh11->SetLineStyle(7); gh11->Draw("P SAME");
		TGraphErrors * g11 = new TGraphErrors(4, r11x, r11y, nullptr, r11dy);
		g11->SetMarkerStyle(21); g11->SetMarkerSize(2.4); g11->SetMarkerColor(cR11); g11->SetLineColor(cR11); g11->SetLineWidth(4); g11->Draw("P SAME");

		TGraphErrors * ga = new TGraphErrors(), * gh = new TGraphErrors();
		for (int s = 0; s < 2; s++)
		{
			auto & v = s == 0 ? all : hi;
			for (auto & p : v)
			{
				if (p.var != var) continue;
				double env = 0;
				for (auto & q : v) if (q.bin == p.bin && find(models[var].begin(), models[var].end(), q.var) != models[var].end()) env = max(env, fabs(q.S / p.S - 1));
				double tot = sqrt(p.sy * p.sy + pow(env * p.y, 2)), x = p.x + (s ? dxHi : 0);
				TBox * b = new TBox(x - 0.07, max(0.0, p.y - tot), x + 0.07, min(16.0, p.y + tot));
				if (s == 0) { b->SetFillColorAlpha(col, 0.22); b->SetLineColor(col); b->SetLineStyle(1); }
				else { b->SetFillStyle(0); b->SetLineColor(col); b->SetLineStyle(7); }
				b->SetLineWidth(3); b->Draw("l");
				TGraphErrors * g = s == 0 ? ga : gh;
				g->SetPoint(g->GetN(), x, p.y); g->SetPointError(g->GetN() - 1, 0, p.st);
			}
		}
		gh->SetMarkerStyle(30); gh->SetMarkerSize(4.2); gh->SetMarkerColor(col); gh->SetLineColor(col); gh->SetLineWidth(4); gh->SetLineStyle(7); gh->Draw("P SAME");
		ga->SetMarkerStyle(29); ga->SetMarkerSize(4.2); ga->SetMarkerColor(col); ga->SetLineColor(col); ga->SetLineWidth(4); ga->Draw("P SAME");

		TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextColor(col); t.SetTextSize(0.062);
		double x0 = k == 0 ? 0.20 : 0.07;
		t.DrawLatex(x0, 0.89, k == 0 ? "Subtraction" : "No subtraction");
		t.SetTextFont(62); t.SetTextColor(kRed + 1); t.SetTextSize(0.05); t.DrawLatex(x0, 0.825, "STAR Preliminary");
		t.SetTextFont(42); t.SetTextColor(kBlack); t.SetTextSize(0.034);
		t.DrawLatex(x0, 0.765, "p+p #sqrt{s} = 510 GeV, BHT2, |y| < 1");
		gh->SetFillStyle(0); ga->SetFillColorAlpha(col, 0.22); g11->SetFillColorAlpha(kGray + 1, 0.3);
		TLegend * lg = new TLegend(x0 - 0.01, 0.25, x0 + 0.40, 0.74); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.037); lg->SetMargin(0.16);
		lg->AddEntry(ga, "p_{T}-integrated", "pef");
		lg->AddEntry(gh, "p_{T} > 4 GeV/c (x+0.1)", "pe");
		lg->AddEntry(g11, "#splitline{Run 11, PRD 112, 032004}{#Upsilon(nS), p_{T}-integrated}", "pef");
		lg->AddEntry(gh11, "#splitline{Run 11, PRD 112, 032004}{#Upsilon(1S), p_{T} > 4 GeV/c}", "pe");
		lg->AddEntry(lin, "y = x", "l");
		lg->Draw();
		pd->RedrawAxis();
	}
	c->SaveAs(out + ".png");
	c->SaveAs(out + ".pdf");
}
