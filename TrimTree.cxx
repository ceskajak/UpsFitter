#include "headers/common.hxx"
#include "headers/setup.hxx"
#include "headers/cuts.hxx"
#include "headers/variations.hxx"

/////////////////////////////////////////////////////////////
// Trims the full tree into one that passes the most relaxed version of every cut (all other cuts at nominal)
// No multiplicity cut is applied, so every bin, TofMult variant and cut variation can run on the trimmed tree
// Rerun whenever a nominal cut value changes (e.g. DCA)
void TrimTree()
{
	TFile * f_full = new TFile(f_BHT2_full, "READ");
	TTree * tree = (TTree *)f_full->Get("tree");
	cout << "Full tree entries: " << tree->GetEntries() << endl;

	// Nominal cut strings, kept to quantify how many events each loosened cut adds
	TString s_cuts_nom = MakeCuts("BHT2");
	vector<pair<TString, TString>> vec_nom = {{"nSigmaE", MakeCutNSigmaE()},
												{"ETOWoECLU", MakeCutETOWoECLU()},
												{"EoP", MakeCutEoP()},
												{"Dphi", MakeCutDphi()},
												{"DZ", MakeCutDZ()},
												{"nHitsFit", MakeCutNHitsFit()},
												{"nHitsRatio", MakeCutNHitsRatio()}};

	SetVariation("loosest");
	TString s_cuts_loose = MakeCuts("BHT2");
	cout << "Loosest cuts: " << s_cuts_loose << endl;

	TFile * f_trim = new TFile(f_BHT2_trim, "RECREATE");
	TTree * tree_trim = tree->CopyTree(s_cuts_loose);
	tree_trim->Write();

	cout << "Trimmed tree entries: " << tree_trim->GetEntries() << endl;
	cout << "  of which pass the nominal cuts: " << tree_trim->GetEntries(s_cuts_nom) << endl;
	for (auto c : vec_nom)
		cout << "  added by loosening " << c.first << " alone: " << tree_trim->GetEntries("!(" + c.second + ")") << endl;

	f_trim->Close();
	f_full->Close();
	cout << "Trimmed tree written to " << f_BHT2_trim << endl;
}

/////////////////////////////////////////////////////////////
