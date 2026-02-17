#include "../setup.hxx"

// in order to fit a new bin/dataset/whatever, you need to define a fitting parameter array and a recursive fraction array
// these look like this:
/*
	Double_t par[8][3] - a fitting parameter array
		like-sign background exponent
		unlike-sign background exponent
		1. Gauss mean
		2. Gauss mean
		3. Gauss mean
		1. Gauss width
		2. Gauss width
		3. Gauss width

	Double_t frac[4][3] - a recursive fraction array
		LS background fraction
		US background fraction
		1. Gauss fraction
		2. Gauss fraction
*/
// where in each row you specify: initial fit value; lower constraint; upper constraint

/////////////////////////////////////////////////////////////
// pT and TofMult integrated
Double_t parBHT2_int[14][3] =
{
	{100.0, -1000.0, 1000.0},	// LS p0
	{5.0, -100.0, 100.0},		// LS p1
	{0.1, -10.0, 10.0},		// LS p2
	{0.1, -10.0, 10.0},		// LS p3
	{100.0, -1000.0, 1000.0},	// US p0
	{1.0, -100.0, 100.0},		// US p1
	{0.1, -10.0, 10.0},		// US p2
	{0.1, -10.0, 10.0},		// US p3
	{3.1, 2.8, 3.5},		// sig1 mean
	{0.1, 0.05, 0.2},		// sig1 width
	{1.4099, 1.386, 1.4338},		// sig1 alpha_L
	{1.8, 0.1, 2.0},		// sig1 alpha_R
	{1.7572, 1.7108, 1.8036},		// sig1 n_L
	{3.1, 0.1, 5.0},		// sig1 n_R
};

Double_t fracBHT2_int[1][3] =
{
	{0.2, 0.0, 1.0},	// bg frac
};

Double_t parBHT2_0mult3[8][3];
Double_t fracBHT2_0mult3[2][3];

Double_t parBHT2_3mult6[8][3];
Double_t fracBHT2_3mult6[2][3];

Double_t parBHT2_6mult9[8][3];
Double_t fracBHT2_6mult9[2][3];

Double_t parBHT2_9mult12[8][3];
Double_t fracBHT2_9mult12[2][3];

Double_t parBHT2_12mult20[8][3];
Double_t fracBHT2_12mult20[2][3];

Double_t parBHT2_9mult20[8][3];
Double_t fracBHT2_9mult20[2][3];

Double_t parBHT2_6mult20[8][3];
Double_t fracBHT2_6mult20[2][3];