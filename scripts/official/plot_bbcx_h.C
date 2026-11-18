void plot_bbcx_h()
{
	const double xref = 3.55, lo = 1.8, hi = 5.0;
	TFile * f = TFile::Open("/mnt/reservoir/Research/data/BBC.histos_250228.root", "READ");
	TProfile * p = (TProfile *)f->Get("p_TofMultxBBC");
	p->SetStats(0);

	TF1 * h1 = new TF1("h1", "pol1", lo, hi);
	TF1 * h2 = new TF1("h2", "pol2", lo, hi);
	p->Fit(h1, "0QR");
	p->Fit(h2, "0QR");
	TF1 * hu = new TF1("hu", "[0] + (0.811027*x - 0.2533*x*x)", lo, hi);
	hu->SetParameter(0, 1.9);
	p->Fit(hu, "0QR");
	printf("pol1: %.4f %.4f  chi2/ndf %.1f/%d  h(xref)=%.4f\n", h1->GetParameter(0), h1->GetParameter(1), h1->GetChisquare(), h1->GetNDF(), h1->Eval(xref));
	printf("pol2: %.4f %.4f %.5f  chi2/ndf %.1f/%d  h(xref)=%.4f\n", h2->GetParameter(0), h2->GetParameter(1), h2->GetParameter(2), h2->GetChisquare(), h2->GetNDF(), h2->Eval(xref));
	printf("UpsFitter shape, offset %.4f  chi2/ndf %.1f/%d  h(xref)=%.4f\n", hu->GetParameter(0), hu->GetChisquare(), hu->GetNDF(), hu->Eval(xref));
	printf("relative slope at xref: pol1 %.3f  pol2 %.3f per 1e6\n", h1->Derivative(xref) / h1->Eval(xref), h2->Derivative(xref) / h2->Eval(xref));

	TCanvas * c = new TCanvas("c", "", 900, 800);
	TPad * pt = new TPad("pt", "", 0, 0.32, 1, 1);
	TPad * pb = new TPad("pb", "", 0, 0, 1, 0.32);
	pt->SetBottomMargin(0.02);
	pb->SetTopMargin(0.02);
	pb->SetBottomMargin(0.30);
	pt->Draw();
	pb->Draw();

	pt->cd();
	p->SetTitle(";;<TofMult>");
	p->SetMarkerStyle(20);
	p->SetMarkerSize(0.7);
	p->GetXaxis()->SetRangeUser(1.5, 5.0);
	p->GetXaxis()->SetLabelSize(0);
	p->GetYaxis()->SetRangeUser(0.0, 3.0);
	p->Draw("E");
	h1->SetLineColor(kBlue + 1);
	h2->SetLineColor(kGreen + 2);
	hu->SetLineColor(kRed + 1);
	TF1 * fs[3] = {h1, h2, hu};
	for (TF1 * g : fs) { g->SetLineWidth(3); g->SetRange(1.5, 5.0); g->Draw("same"); }
	TLegend * lg = new TLegend(0.40, 0.62, 0.89, 0.89);
	lg->SetBorderSize(0);
	lg->AddEntry(p, "BBC MB data", "lp");
	lg->AddEntry(h1, Form("pol1 fit, chi2/ndf %.1f", h1->GetChisquare() / h1->GetNDF()), "l");
	lg->AddEntry(h2, Form("pol2 fit, chi2/ndf %.1f", h2->GetChisquare() / h2->GetNDF()), "l");
	lg->AddEntry(hu, Form("UpsFitter shape, offset only, chi2/ndf %.0f", hu->GetChisquare() / hu->GetNDF()), "l");
	lg->Draw();

	pb->cd();
	TH1D * r1 = (TH1D *)p->ProjectionX("r1");
	TH1D * r2 = (TH1D *)p->ProjectionX("r2");
	TH1D * ru = (TH1D *)p->ProjectionX("ru");
	TF1 * hs[3] = {h1, h2, hu};
	TH1D * rs[3] = {r1, r2, ru};
	for (int k = 0; k < 3; k++)
		for (int i = 1; i <= p->GetNbinsX(); i++)
			if (p->GetBinEntries(i) > 0)
			{
				double x = p->GetBinCenter(i), g = hs[k]->Eval(x);
				rs[k]->SetBinContent(i, p->GetBinContent(i) / g);
				rs[k]->SetBinError(i, p->GetBinError(i) / g);
			}
			else { rs[k]->SetBinContent(i, 0); rs[k]->SetBinError(i, 0); }
	Color_t col[3] = {kBlue + 1, kGreen + 2, kRed + 1};
	for (int k = 0; k < 3; k++)
	{
		rs[k]->SetMarkerStyle(20);
		rs[k]->SetMarkerSize(0.5);
		rs[k]->SetMarkerColor(col[k]);
		rs[k]->SetLineColor(col[k]);
	}
	r1->SetTitle(";BBCx [10^{6}];data / curve");
	r1->GetXaxis()->SetRangeUser(1.8, 5.0);
	r1->GetYaxis()->SetRangeUser(0.7, 1.3);
	r1->GetXaxis()->SetTitleSize(0.11);
	r1->GetXaxis()->SetLabelSize(0.09);
	r1->GetYaxis()->SetTitleSize(0.10);
	r1->GetYaxis()->SetTitleOffset(0.45);
	r1->GetYaxis()->SetLabelSize(0.08);
	r1->Draw("E");
	r2->Draw("E same");
	ru->Draw("E same");
	TLine * l = new TLine(1.8, 1, 5.0, 1);
	l->SetLineStyle(2);
	l->Draw();

	c->SaveAs("bbcx_h.png");
	c->SaveAs("bbcx_h.pdf");
}
