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
Double_t parBHT1_int[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.3, 9.3, 9.35},		// sig1 mean
	{9.9, 9.9, 10.1},	// sig2 mean
	{10.3, 10.2, 10.35},	// sig3 mean
	{0.2, 0.15, 0.30},	// sig1 width
	{0.32, 0.31, 0.32},	// sig2 width
	{0.336, 0.32, 0.33},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_int[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.72, 1.0},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.68, 0.6, 0.7}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 0 <= TofMult < 3
Double_t parBHT1_0mult3[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.3, 9.3, 9.35},		// sig1 mean
	{9.9, 9.9, 10.1},	// sig2 mean
	{10.3, 10.2, 10.35},	// sig3 mean
	{0.2, 0.15, 0.30},	// sig1 width
	{0.32, 0.31, 0.32},	// sig2 width
	{0.336, 0.32, 0.33},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_0mult3[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.4, 1.0},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.8, 0.75, 0.9},	// sig1 frac
	{0.68, 0.6, 0.7}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 3 <= TofMult < 6
Double_t parBHT1_3mult6[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.3, 9.3, 9.35},		// sig1 mean
	{9.9, 9.9, 10.0},	// sig2 mean
	{10.3, 10.25, 10.35},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.3, 0.31, 0.32},	// sig2 width
	{0.33, 0.32, 0.33},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_3mult6[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.61, 0.5, 1.0},	// LS bg frac
	{0.51, 0.4, 1.0},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.68, 0.6, 0.8}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 6 <= TofMult < 9
Double_t parBHT1_6mult9[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.3, 9.3, 9.35},		// sig1 mean
	{9.9, 9.9, 10.0},	// sig2 mean
	{10.3, 10.25, 10.35},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.3, 0.3, 0.31},	// sig2 width
	{0.33, 0.32, 0.33},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_6mult9[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.7, 0.67, 0.72},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.75, 0.6, 0.8},	// sig1 frac
	{0.68, 0.65, 0.8}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 9 <= TofMult < 12
Double_t parBHT1_9mult12[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.442, 9.3, 9.5},		// sig1 mean
	{10.006, 9.9, 10.1},	// sig2 mean
	{10.334, 10.2, 10.5},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.32, 0.31, 0.33},	// sig2 width
	{0.336, 0.33, 0.34},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_9mult12[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.72, 1.0},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.68, 0.6, 0.8}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 12 <= TofMult < 20
Double_t parBHT1_12mult20[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.442, 9.3, 9.5},		// sig1 mean
	{10.006, 9.9, 10.1},	// sig2 mean
	{10.334, 10.2, 10.5},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.32, 0.31, 0.33},	// sig2 width
	{0.336, 0.33, 0.34},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_12mult20[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.72, 1.0},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.68, 0.6, 0.8}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 9 <= TofMult < 20
Double_t parBHT1_9mult20[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.25, 9.25, 9.35},		// sig1 mean
	{9.9, 9.9, 9.95},	// sig2 mean
	{10.3, 10.25, 10.35},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.3, 0.31, 0.32},	// sig2 width
	{0.33, 0.32, 0.33},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_9mult20[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.6, 0.75},	// LS bg frac
	{0.63, 0.5, 0.85},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.7, 0.6, 0.8}	// sig2 frac
};

/////////////////////////////////////////////////////////////
// 6 <= TofMult < 20
Double_t parBHT1_6mult20[16][3] =
{
	{-0.2, -0.5, -0.0},		// LS bg exponent
	{0.2, 0.0, 2.0},		// LS bg exp corr
	{-0.2, -0.5, -0.0},		// US bg exponent
	{0.2, 0.0, 2.0},		// US bg exp corr
	{9.442, 9.3, 9.5},		// sig1 mean
	{10.006, 9.9, 10.1},	// sig2 mean
	{10.334, 10.2, 10.5},	// sig3 mean
	{0.2, 0.15, 0.3},	// sig1 width
	{0.32, 0.31, 0.33},	// sig2 width
	{0.336, 0.33, 0.34},	// sig3 width
	{1.71, 1.66, 1.76},		// sig1 alpha
	{1.71, 1.65, 1.77},		// sig2 alpha
	{1.77, 1.73, 1.81},		// sig3 alpha
	{1.62, 1.48, 1.76},		// sig1 n
	{1.57, 1.44, 1.70},		// sig2 n
	{1.36, 1.27, 1.45}		// sig3 n
};

Double_t fracBHT1_6mult20[6][3] =
{
	{0.5, 0.0, 1.0},	// LS corr frac
	{0.5, 0.0, 1.0},	// US corr frac
	{0.74, 0.72, 1.0},	// LS bg frac
	{0.63, 0.5, 1.0},	// US bg frac
	{0.71, 0.6, 0.8},	// sig1 frac
	{0.68, 0.6, 0.8}	// sig2 frac
};