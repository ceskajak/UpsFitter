#!/bin/bash

# Runs every variation needed for the full systematics (see Upsilon_systematics in the vault) for all systematics bins
# Usage: ./runSyst.sh [-j jobs] [-f]
# Variation names are defined in headers/variations.hxx (cut names: vec_cuts, no loose variant for vec_cuts_noloose)

cd "$(dirname "$0")"

CUTS="nHitsFit nHitsRatio DCA nSigmaE_lo nSigmaE_hi ETOWoECLU EoP_lo EoP_hi Dphi DZ ADC"
NOLOOSE="DCA nSigmaE_hi ADC"

MULTBINS="0mult1 1mult2 2mult3 3mult4 4mult8 3mult8"
VARS="range_lo_up range_lo_dn range_hi_up range_hi_dn lineshape_up lineshape_dn ratio_up ratio_dn meanfree sepscale"
for c in $CUTS; do
	VARS="$VARS cut_${c}_tight"
	case " $NOLOOSE " in
		*" $c "*) ;;
		*) VARS="$VARS cut_${c}_loose" ;;
	esac
done

# 1S mean: fixed at the nominal integrated-bin value in all other bins, meanfree floats it; sepscale: 2S and 3S distances scaled by a floating parameter; sepcon: the same, constrained to the sepscale fit of the integrated bin
# Exponential background: all variations; Chebyshev background: nominal only
# The nominal fits come first: every variation starts from the fitted values of its nominal bin
# The integrated bin comes first: every other bin reads its 1S mean from it
./runAll.sh -b int -v "nom" "$@"
./runAll.sh -b "$MULTBINS" -v "nom" "$@"
./runAll.sh -b int -v sepscale "$@"
./runAll.sh -v "$VARS" "$@"
# sepcon reads the sepscale fit of the integrated bin, so it runs after it is finished
./runAll.sh -v sepcon "$@"
./runAll.sh -b int -v "nom" -g cheby "$@"
./runAll.sh -b "$MULTBINS" -v "nom" -g cheby "$@"

# TofMult variants (nominal only, the integrated bin has no multiplicity cut and is identical for all of them)
for m in sub corr none; do
	./runAll.sh -b "$MULTBINS" -v "nom" -m $m "$@"
done

# Status table and plots of everything run today, failed or missing fits are listed first
root -l -b -q 'Summary.cxx("", "subcorr")' 2>&1 | grep -v "^Info in"

# Systematic per source and bin, total per bin
root -l -b -q 'Syst.cxx("", "subcorr")' 2>&1 | grep -v "^Info in"
