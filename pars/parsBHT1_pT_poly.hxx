#include "RooRealVar.h"

// in order to fit a new bin/dataset/whatever, you need to define a fitting parameter array and a recursive fraction array
// these look like this:
/*
	Double_t par[18][3] - a fitting parameter array
		like-sign p0
		like-sign p1
		like-sign p2
		unlike-sign p0
		unlike-sign p1
		unlike-sign p2
		1. Gauss mean
		2. Gauss mean
		3. Gauss mean
		1. Gauss width
		2. Gauss width
		3. Gauss width
		1. Gauss alpha
		2. Gauss alpha
		3. Gauss alpha
		1. Gauss n
		2. Gauss n
		3. Gauss n

	Double_t frac[4][3] - a recursive fraction array
		LS background fraction
		US background fraction
		1. Gauss fraction
		2. Gauss fraction
*/
// where in each row you specify: initial fit value; lower constraint; upper constraint

/////////////////////////////////////////////////////////////
// 0 <= pT < 2
Double_t parBHT1_0pT2[18][3] =
{
	{0, -1000, 1000},		// LS bg p0
	{0, -1000, 1000},		// LS bg p1
	{0, -1000, 1000},		// LS bg p2
	{0, -1000, 1000},		// US bg p0
	{0, -1000, 1000},		// US bg p1
	{0, -1000, 1000},		// US bg p2
	{9.442, 9.440, 9.444},		// sig1 mean
	{10.006, 10.008, 10.010},	// sig2 mean
	{10.334, 10.332, 10.336},	// sig3 mean
	{0.194, 0.193, 0.195},	// sig1 width
	{0.208, 0.206, 0.210},	// sig2 width
	{0.214, 0.211, 0.217},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_0pT2[4][3] =
{
	{0.7, 0.5, 0.9},	// LS bg frac
	{0.7, 0.5, 0.9},	// US bg frac
	{0.71, 0.7, 0.72},	// sig1 frac
	{0.68, 0.67, 0.69}	// sig2 frac
};