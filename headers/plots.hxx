#pragma once

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TStyle.h"

/*
/////////////////////////////////////////////////////////////
PLOT HELPERS (Summary.cxx, Syst.cxx)
/////////////////////////////////////////////////////////////
*/

vector<Int_t> vec_colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1, kCyan + 2};

// Empty frame with one labelled bin per entry of labels
static TH1F * MakeFrame(TString name, vector<TString> labels, Double_t ymin, Double_t ymax, TString ytitle)
{
	TH1F * h = new TH1F(name, "", labels.size(), 0, labels.size());
	for (int i = 0; i < (int)labels.size(); i++) h->GetXaxis()->SetBinLabel(i + 1, labels[i]);
	h->GetXaxis()->LabelsOption("v");
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetTitle(ytitle);
	h->GetYaxis()->SetTitleSize(0.045);
	h->GetYaxis()->SetLabelSize(0.04);
	h->SetMinimum(ymin);
	h->SetMaximum(ymax);
	h->SetStats(0);
	return h;
}

static void SetPoint(TGraphErrors * g, Int_t color, Int_t style, Double_t size)
{
	g->SetMarkerColor(color);
	g->SetLineColor(color);
	g->SetMarkerStyle(style);
	g->SetMarkerSize(size);
}

static void Title(TString text)
{
	TLatex t;
	t.SetNDC();
	t.SetTextFont(42);
	t.SetTextSize(0.022);
	t.DrawLatex(0.01, 0.985, text);
}

/////////////////////////////////////////////////////////////
