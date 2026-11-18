// Multiplicity dependence (subcorr) with the BBCx model envelope added in quadrature to the systematic boxes
void syst_with_bbcx(TString txt = "/mnt/reservoir/Research/repos/UpsFitter/output/multdep/multdep_261105_dca1.txt", TString out = "/home/jakub/workspace/scratch/bbcx_fit/syst_with_bbcx")
{
	struct P { string var, bin; double x, y, st, sy, S; };
	map<string, vector<P>> sch;
	string line, cur;
	ifstream f(txt.Data());
	while (getline(f, line))
	{
		if (line.rfind("Scheme ", 0) == 0) { cur = line.substr(7); continue; }
		istringstream is(line);
		P p;
		double nmb;
		if (is >> p.var >> p.bin >> p.x >> p.y >> p.st >> p.sy >> p.S >> nmb) sch[cur].push_back(p);
	}
	vector<string> models = {"subcorrm", "mbmul1", "mbmul2", "mbadd1", "mbadd2"};
	Double_t r11x[] = {0.53835, 1.3283, 2.2807, 3.3673}, r11y[] = {0.41, 1.39, 2.67, 4.6}, r11dy[] = {0.064, 0.14, 0.43, 1.8}, r11ds[] = {0.041, 0.085, 0.31, 1.3};
	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c", "c", 1800, 800);
	c->Divide(2, 1);
	int ip = 0;
	for (string name : {"4bin", "5bin"})
	{
		c->cd(++ip);
		gPad->SetLeftMargin(0.1);
		auto & v = sch[name];
		vector<P> nom;
		for (auto & p : v) if (p.var == "subcorr") nom.push_back(p);
		TH1F * fr = new TH1F(Form("fr%d", ip), "", 1, 0, 5);
		fr->SetMinimum(0); fr->SetMaximum(6.5);
		fr->GetXaxis()->SetTitle("N_{ch}/<N_{ch}>"); fr->GetYaxis()->SetTitle("N_{#Upsilon}^{bin}/<N_{#Upsilon}>");
		fr->Draw();
		TF1 * lin = new TF1(Form("lin%d", ip), "x", 0, 5); lin->SetLineStyle(2); lin->SetLineColor(kGray + 2); lin->Draw("same");
		for (int i = 0; i < 4; i++)
		{
			TBox * b = new TBox(r11x[i] - 0.07, r11y[i] - r11ds[i], r11x[i] + 0.07, r11y[i] + r11ds[i]); b->SetFillColorAlpha(kGray + 1, 0.35); b->SetLineColor(kGray + 1); b->Draw();
		}
		TGraphErrors * g11 = new TGraphErrors(4, r11x, r11y, nullptr, r11dy); g11->SetMarkerStyle(21); g11->Draw("P SAME");
		TGraphErrors * g = new TGraphErrors();
		TBox * lb1 = nullptr, * lb2 = nullptr;
		for (auto & p : nom)
		{
			double env = 0;
			for (auto & q : v)
				if (q.bin == p.bin && find(models.begin(), models.end(), q.var) != models.end()) env = max(env, fabs(q.S / p.S - 1));
			double tot = sqrt(p.sy * p.sy + pow(env * p.y, 2));
			TBox * b2 = new TBox(p.x - 0.11, p.y - tot, p.x + 0.11, p.y + tot); b2->SetFillColorAlpha(kRed + 1, 0.2); b2->SetLineColor(kRed + 1); b2->Draw();
			TBox * b1 = new TBox(p.x - 0.07, p.y - p.sy, p.x + 0.07, p.y + p.sy); b1->SetFillColorAlpha(kBlue + 1, 0.35); b1->SetLineColor(kBlue + 1); b1->Draw();
			if (!lb1) { lb1 = b1; lb2 = b2; }
			g->SetPoint(g->GetN(), p.x, p.y); g->SetPointError(g->GetN() - 1, 0, p.st);
			bool left = p.x > 3.3; TLatex * t = new TLatex(p.x + (left ? -0.16 : 0.16), p.y - 0.02 * 6.5 - (p.x > 4 ? 1.2 : 0), Form("syst %.2f #rightarrow %.2f (env %.0f%%)", p.sy, tot, 100 * env)); t->SetTextSize(0.024); t->SetTextAlign(left ? 31 : 11); t->Draw();
		}
		g->SetMarkerStyle(20); g->SetMarkerColor(kBlue + 1); g->SetLineColor(kBlue + 1); g->Draw("P SAME");
		TLegend * lg = new TLegend(0.12, 0.72, 0.65, 0.92); lg->SetBorderSize(0); lg->SetTextSize(0.03);
		lg->AddEntry(g11, "Run11 (stat bars, syst boxes)", "pe");
		lg->AddEntry(g, Form("Run17 %s subcorr, stat bars", name.c_str()), "pe");
		lg->AddEntry(lb1, "syst (Syst.cxx)", "f");
		lg->AddEntry(lb2, "syst + BBCx model envelope", "f");
		lg->Draw();
	}
	c->SaveAs(out + ".png");
	c->SaveAs(out + ".pdf");
}
