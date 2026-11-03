#pragma once

#include "common.hxx"
#include "setup.hxx"
#include "cuts.hxx"

/////////////////////////////////////////////////////////////

Double_t kRANGE_SHIFT = 0.10;	// relative shift of the fit range edges
Double_t kCUT_SHIFT = 0.10;		// relative shift of the selection cuts

// Cuts that can be varied, name as used in the variation string (cut_<name>_tight / cut_<name>_loose)
// DCA, nSigmaE_hi and ADC have no loose variation
vector<TString> vec_cuts = {"nHitsFit", "nHitsRatio", "DCA", "nSigmaE_lo", "nSigmaE_hi", "ETOWoECLU", "EoP_lo", "EoP_hi", "Dphi", "DZ", "ADC"};
vector<TString> vec_cuts_noloose = {"DCA", "nSigmaE_hi", "ADC"};

/////////////////////////////////////////////////////////////

// Shifts a cut value by kCUT_SHIFT, tight = fewer events: lower bounds go up, upper bounds go down
// Integer cuts are rounded to the tight side (tight rounds up for lower bounds, loose rounds down)
static float ShiftCut(float val, bool lower, bool tight, bool integer = false)
{
	bool up = (lower == tight);
	float shifted = val + (up ? 1 : -1) * fabs(val) * kCUT_SHIFT;
	if (integer) shifted = up ? ceil(shifted - 1e-4) : floor(shifted + 1e-4);
	return shifted;
}

/////////////////////////////////////////////////////////////

// Applies the named variation to the run settings (call before the cuts and datasets are built)
// Returns false for an unknown or not yet implemented name, so a run never silently falls back to nominal
static bool SetVariation(TString kVAR)
{
	if (kVAR == "nom")
		return true;

	if (kVAR.BeginsWith("range_"))
	{
		if (kVAR == "range_lo_up") kMASS_MIN *= (1 + kRANGE_SHIFT);
		else if (kVAR == "range_lo_dn") kMASS_MIN *= (1 - kRANGE_SHIFT);
		else if (kVAR == "range_hi_up") kMASS_MAX *= (1 + kRANGE_SHIFT);
		else if (kVAR == "range_hi_dn") kMASS_MAX *= (1 - kRANGE_SHIFT);
		else return false;

		mass.setRange(kMASS_MIN, kMASS_MAX);
		mass.setRange("full", kMASS_MIN, kMASS_MAX);
		return true;
	}

	// Lineshape and ratio: the parameters are fixed at the BHT2_int bounds inside fit_bin
	if (kVAR == "lineshape_up" || kVAR == "lineshape_dn")
	{
		kSHAPE_FIX = (kVAR == "lineshape_up") ? 1 : -1;
		return true;
	}
	if (kVAR == "ratio_up" || kVAR == "ratio_dn")
	{
		kRATIO_FIX = (kVAR == "ratio_up") ? 1 : -1;
		return true;
	}

	// Handled inside fit_bin: meanfree floats the 1S mean in every bin (nominal: fixed at the integrated-bin value), sepscale floats the 2S and 3S distance scale
	if (kVAR == "meanfree")
		return true;

	// Diagnostic, not a systematic source: upper bound of the 1S mean, e.g. meanmax_10.0
	if (kVAR.BeginsWith("meanmax_"))
	{
		kMEAN_MAX = TString(kVAR(8, kVAR.Length() - 8)).Atof();
		return kMEAN_MAX > 0;
	}

	// Option, not a systematic source: 2S and 3S mass distances to the 1S mean scaled by one floating parameter
	if (kVAR == "sepscale")
	{
		kSEP_SCALE = true;
		return true;
	}

	// 2S and 3S distance scale floating with a Gaussian constraint (mean and width from the sepscale fit of the integrated bin, read in fit_bin)
	if (kVAR == "sepcon")
	{
		kSEP_SCALE = true;
		kSEP_CON = true;
		return true;
	}

	// All cuts at their loosest variation (used for the tree trimmer)
	if (kVAR == "loosest")
	{
		for (auto c : vec_cuts)
			if (find(vec_cuts_noloose.begin(), vec_cuts_noloose.end(), c) == vec_cuts_noloose.end())
				SetVariation("cut_" + c + "_loose");
		return true;
	}

	if (kVAR.BeginsWith("cut_"))
	{
		bool tight = kVAR.EndsWith("_tight");
		bool loose = kVAR.EndsWith("_loose");
		if (!tight && !loose) return false;

		TString name = kVAR(4, kVAR.Length() - 4 - 6);
		if (find(vec_cuts.begin(), vec_cuts.end(), name) == vec_cuts.end()) return false;
		if (loose && find(vec_cuts_noloose.begin(), vec_cuts_noloose.end(), name) != vec_cuts_noloose.end()) return false;

		if (name == "nHitsFit") kNHITSFIT_LO = ShiftCut(kNHITSFIT_LO, true, tight, true);
		else if (name == "nHitsRatio") kNHITSRATIO_LO = ShiftCut(kNHITSRATIO_LO, true, tight);
		else if (name == "DCA") kDCA_HI = ShiftCut(kDCA_HI, false, tight);
		else if (name == "nSigmaE_lo") kSIGMAE_LO = ShiftCut(kSIGMAE_LO, true, tight);
		else if (name == "nSigmaE_hi") kSIGMAE_HI = ShiftCut(kSIGMAE_HI, false, tight);
		else if (name == "ETOWoECLU") kETOWOECLU_LO = ShiftCut(kETOWOECLU_LO, true, tight);
		else if (name == "EoP_lo") kEOP_LO = ShiftCut(kEOP_LO, true, tight);
		else if (name == "EoP_hi") kEOP_HI = ShiftCut(kEOP_HI, false, tight);
		else if (name == "Dphi") kDphi_HI = ShiftCut(kDphi_HI, false, tight);
		else if (name == "DZ") kDZ_HI = ShiftCut(kDZ_HI, false, tight);
		else if (name == "ADC")
		{
			kADC_TRIG_BHT1 = ShiftCut(kADC_TRIG_BHT1, true, tight, true);
			kADC_TRIG_BHT2 = ShiftCut(kADC_TRIG_BHT2, true, tight, true);
		}
		return true;
	}

	return false;
}

/////////////////////////////////////////////////////////////
