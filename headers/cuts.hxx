#pragma once

#include "common.hxx"

/////////////////////////////////////////////////////////////

const float kY_LO = -1;			    // lower bound for quarkonium rapidity
const float kY_HI = 1;			    // upper bound for quarkonium rapidity

const float kSIGMAE_LO = -1.2;	    // lower bound for nSigmaElectron (same for both electrons)

const float kETOWOECLU_LO = 0.5;	// lower bound for E_TOW/E_CLU (same for both electrons)

const float kEOP_LO = 0.7;	        // lower bound for E_CLU/p (same for both electrons)
const float kEOP_HI = 1.45;         // upper bound for E_CLU/p (same for both electrons)

const float kPAR_PT_LO = 1;	        // lower bound for electron1 pT

const int kADC_TRIG_BHT1 = 18;	    // DSM-ADC for trigger electron (BHT1 dataset)
const int kADC_TRIG_BHT2 = 25;	    // DSM-ADC for trigger electron (BHT2 dataset)
const int kAED_PAR = -1;	        // DSM-ADC for partner electron. off due to greater than

const float kDZ_HI = 5;			    // ZDist
const float kDphi_HI = 0.05;	    // PhiDist

// Systematic variations
const float kNHITSFIT_LO = 20;		// nHitsFit for electron

const float kNHITSRATIO_LO = 0.52;	// nHitsFit/nHitsMax for electron

const float kDCA_HI = 1;			// DCA for electron

const float kVZDIFF_HI = 9999;		// maximum difference between VzTPC and VzVPD

/////////////////////////////////////////////////////////////

TString s_cut_VzDiff = Form("(abs(VzTPC - VzVPD) <= %f)", kVZDIFF_HI);

TString s_cut_BEMC_valid = Form("E_TOW1 > 0 && E_CLU1 > 0 && E_TOW2 > 0 && E_CLU2 > 0");

TString s_cut_y = Form("(y > %f && y < %f)", kY_LO, kY_HI);
TString s_cut_nSigmaE = Form("(nSigmaElectron1 >= %f && nSigmaElectron2 >= %f)", kSIGMAE_LO, kSIGMAE_LO);
TString s_cut_ETOWoECLU = Form("(E_TOW1/E_CLU1 >= %f && E_TOW2/E_CLU2 >= %f)", kETOWOECLU_LO, kETOWOECLU_LO);
TString s_cut_EoP = Form("(E_CLU1/p1 > %f && E_CLU1/p1 < %f && E_CLU2/p2 > %f && E_CLU2/p2 < %f)",
                               kEOP_LO, kEOP_HI, kEOP_LO, kEOP_HI);
static TString MakeCutADC(int adc)
{
    return Form("((adc1 >= %f && pt2 >= %f) || (adc2 >= %f && pt1 >= %f))", adc * 16.0, kPAR_PT_LO, adc * 16.0, kPAR_PT_LO);
}

TString s_cut_DZ = Form("(abs(bemcZDist1) <= %f && abs(bemcZDist2) <= %f)", kDZ_HI, kDZ_HI);
TString s_cut_Dphi = Form("(abs(bemcPhiDist1) <= %f && abs(bemcPhiDist2) <= %f)", kDphi_HI, kDphi_HI);

TString s_cut_nHitsFit = Form("(nHitsFit1 >= %f && nHitsFit2 >= %f)", kNHITSFIT_LO, kNHITSFIT_LO);
TString s_cut_nHitsRatio = Form("(nHitsFit1/nHitsMax1 >= %f && nHitsFit2/nHitsMax2 >= %f)", kNHITSRATIO_LO, kNHITSRATIO_LO);
TString s_cut_DCA = Form("(abs(DCA1) <= %f && abs(DCA2) <= %f)", kDCA_HI, kDCA_HI);

// TofMult variants, evaluated as tree expressions: daughter-electron subtraction and luminosity (BBCx) correction on/off
const TString s_mult_none = "TofMult";
const TString s_mult_sub = "(TofMult - Sum$(TofMult_id == id1 || TofMult_id == id2))";
const TString s_mult_corr = "(TofMult - BBCx/1000000.0 * 0.811027 - pow(BBCx/1000000.0, 2) * -0.2533)";
const TString s_mult_subcorr = "(TofMult - Sum$(TofMult_id == id1 || TofMult_id == id2) - BBCx/1000000.0 * 0.811027 - pow(BBCx/1000000.0, 2) * -0.2533)";

TString s_mult = s_mult_subcorr;	// TofMult variant used for the multiplicity bins

TString s_cut_mult;

/////////////////////////////////////////////////////////////

static TString AppendCuts(const TString &cut1, const TString &cut2)
{
    if (cut1.Length() == 0)
        return cut2;
    if (cut2.Length() == 0)
        return cut1;
    TString combined = cut1;
    combined += " && ";
    combined += cut2;
    return combined;
}

/////////////////////////////////////////////////////////////

static TString CombineCuts(std::initializer_list<TString> cuts)
{
    TString combined;
    for (const auto &cut : cuts)
    {
        if (cut.Length() > 0)
            combined = AppendCuts(combined, cut);
    }
    return combined;
}

/////////////////////////////////////////////////////////////
