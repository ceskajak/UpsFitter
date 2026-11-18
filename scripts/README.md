# scripts

ROOT macros that produce the final (presentation) plots. Copied unchanged on 2026-11-18 from `~/workspace/scratch/bbcx_fit/`, where they were developed. Run with `root -l -b -q`.

## official/

All macros read the result files of `MultDep.cxx` (`output/multdep/multdep_<YYMMDD>_dca1[_pt4].txt`, one block per scheme `4bin` / `5bin`, columns `var bin x y stat syst S N_MB`). Boxes are `Syst.cxx` systematic and the BBCx model envelope in quadrature; the envelope is the maximum of `|S_model/S - 1|` over the five BBCx models of the same subtraction choice. Run17 points are stars, Run11 is the reference.

Variants: `subcorr` = daughters subtracted (models `subcorrm`, `mbmul1`, `mbmul2`, `mbadd1`, `mbadd2`); `corr` = no subtraction (models `corrm`, `mbmul1ns`, `mbmul2ns`, `mbadd1ns`, `mbadd2ns`).

| Macro | Input | Output (`out` + suffix, PNG and PDF) | Deck plot (`UpsFitter_systematics_261112/plots/official/`) |
|---|---|---|---|
| `official_plot_261112.C` | `multdep_261109_dca1.txt` | `official_261112_{subcorr,corr}_{4bin,5bin}` | `multdep_{subcorr,corr}_{4bin,5bin}` |
| `official_plot_261112_pt4.C` | `multdep_261109_dca1_pt4.txt` | `official_261112_pt4_{subcorr,corr}_{4bin,5bin}` | `multdep_pt4_{subcorr,corr}_{4bin,5bin}` |
| `official_plot_261112_pt4_both.C` | `multdep_261109_dca1_pt4.txt` | `official_261112_pt4_{4bin,5bin}` (both sets in one plot) | `multdep_full_syst_5bin_pt4` |
| `official_plot_261112_ptcomp.C` | both txt files above, 5-bin | `official_261112_ptcomp_5bin` (left subtraction, right no subtraction) | `multdep_5bin_ptcomp` |
| `official_plot_261109.C` | `multdep_261109_dca1.txt` | `official_261109_{4bin,5bin}` (both sets in one plot, own envelope each) | `multdep_full_syst_{4bin,5bin}` |
| `syst_with_bbcx.C` | `multdep_261105_dca1.txt` | `syst_with_bbcx` (subcorr, boxes with and without the envelope) | `bbcx_correction/syst_with_bbcx` |
| `plot_bbcx_h.C` | `/mnt/reservoir/Research/data/BBC.histos_250228.root` | `bbcx_h` (MB profile fits with the correction models, written to the current directory) | `bbcx_correction/mb_fits_models` |

Run, from `scripts/official/` (the `out` path defaults to `~/workspace/scratch/bbcx_fit/`, so pass `out` or edit the default):

```
root -l -b -q 'official_plot_261112.C("../../output/multdep/multdep_261109_dca1.txt","official_261112")'
```

`official_plot_261112_ptcomp.C` takes only `out` (it reads both txt files from `/mnt/reservoir/Research/repos/UpsFitter/output/multdep/`); `plot_bbcx_h.C` takes no arguments.

## Run11 reference

- pT-integrated plots: Υ(nS) values `y = 0.41, 1.39, 2.67, 4.6`, stat `0.064, 0.14, 0.43, 1.8`, syst `0.041, 0.085, 0.31, 1.3` at `x = 0.53835, 1.3283, 2.2807, 3.3673`.
- pT > 4 GeV/c plots (decided 2026-11-12): third table of `HEPData-ins2877008-v1-Figure_10a.csv` (Υ(1S), pT > 4 GeV/c): `y = 0.277, 1.33, 2.65, 11.1`, stat `0.08, 0.2, 0.65, 3.4`, syst `0.057, 0.1, 0.39, 2.7`.
- The values are typed into each macro; the CSV is not read. Legends: "Run 11, PRD 112, 032004" plus the selection.

## Not in this directory

- `MultDep.cxx`, `Syst.cxx`, `Summary.cxx`, `BBCxCorr.cxx` (repo root) produce the inputs; `BBCxCorr.cxx` writes `output/bbcx/{mb_profiles,ups_profiles,factor,ups_tofmult}`, `Syst.cxx` the `output/summary/syst_<date>_dca1_<variant>_{1_sources,2_fraction}` figures, `MultDep.cxx` the `output/multdep/` plots, `UpsFitter.cxx` the fit figures in `output/fig/`.
- Deck assembly: `/mnt/reservoir/Research/Presentations/UpsFitter_systematics_261112/scripts/collect_plots.py` (copies the figures into `plots/<topic>/`; its `DATE` is `261105` for the syst and multdep plots) and `build.py`.
- Superseded and not copied: `official_plot.C` (261105), `official_plot_261109_pt4.C`, `official_plot_261109_ptcomp.C`.
