# UpsFitter
## Overview

UpsFitter is a C++ ROOT-based analysis tool designed to fit the invariant mass distribution of Upsilon candidates from [CandReader](https://github.com/ceskajak/CandReader). It uses RooFit to perform simultaneous fits on unlike-sign (signal) and like-sign (background) pair samples across multiple multiplicity and $p_T$ bins.

## Project Structure

- **fit_simul.cxx** - Main fitting macro that loads data and orchestrates the fitting procedure
- **setup.hxx** - Configuration header with global variables, file paths, and fitting parameters
- **pars/** - Parameter directory containing initial fit values and constraints for different triggers:
  - `parsBHT1.hxx` - Parameters for BHT1 trigger
  - `parsBHT2.hxx` - Parameters for BHT2 trigger
  - Additional variants with `_pT` and `_Jpsi` suffixes for specialized analyses

## Usage

### Prerequisites
- ROOT framework with RooFit
- Input data files (specified in `setup.hxx`)

### Quick Start

1. **Initialize output directories:**
   ```bash
   ./init.sh
   ```

2. **Run fits for a single bin:**
   ```bash
   root -l -b -q 'fit_simul.cxx("BHT2", "int")'
   ```
   where `"BHT2"` is the trigger (BHT1 or BHT2) and `"int"` is the bin name (e.g., "int", "0mult3", "3mult6", etc.)

3. **Run fits for all bins in parallel:**
   ```bash
   ./runAll.sh
   ```

## Fitting Approach

The code performs simultaneous fits of the invariant mass distribution using:

- **Signal Model**: Multiple CrystalBall (typically 3) to model the different Upsilon resonance states
- **Background Model**: Exponential functions for both unlike-sign and like-sign backgrounds
- **Parameters**: Each bin contains:
  - 8 parameters for background and signal shapes
  - 4 fraction parameters to control relative contributions

## Configuration

Edit `setup.hxx` to customize:
- Input file paths
- Mass range for fitting (`min_mass`, `max_mass`)
- Bin width for histograms (`bin_width`)
- Output tag and naming conventions

Edit parameter files in `pars/` to adjust initial fit values and constraints for different triggers and kinematic regions.

## Output

Generated output is organized in the `output/` directory:
- `fig/` - Fit plots and canvas files
- `results/` - Fit results and summaries
- `corr/` - Correlation matrices
- `rootfiles/` - ROOT files with fit results

## Notes

- The code processes data split into multiple bins (typically by multiplicity or pT ranges)
- Simultaneous fitting of US and LS samples improves background constraint
- Parameter files contain detailed comments explaining the meaning of each parameter

